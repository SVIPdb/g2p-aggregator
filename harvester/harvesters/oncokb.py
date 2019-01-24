#!/usr/bin/python

from pathlib import Path
import pandas as pd
import requests
from urllib import quote_plus

import cosmic_lookup_table

from lookups import evidence_label as el, evidence_direction as ed
from normalizers.gene_enricher import get_gene
from utils import str_or_none


# OncoKB harvester now pulls from the below downloadable OncoKB files
# and supplements with additional variant data pulled from their public
# API. This is because pulling from the private API gives unpredictable
# results and their is no endpoint in the public API that gives the
# same drug-gene-variant association information as was being
# pulled from the private API.
clinv = Path('../data/oncokb_allActionableVariants.txt')
biov = Path('../data/oncokb_allAnnotatedVariants.txt')


# used to get COSMIC info about genes/alterations
LOOKUP_TABLE = cosmic_lookup_table.CosmicLookup("./cosmic_lookup_table.tsv")


def harvest(genes):
    levels = requests.get('http://oncokb.org/api/v1/levels').json()

    if not genes:
        all_genes = requests.get('http://oncokb.org/api/v1/genes').json()
        genes = set(gene['hugoSymbol'] for gene in all_genes)

    # get all variants
    print 'Gathering all OncoKB variants'
    variants = requests.get('http://oncokb.org/api/v1/variants').json()

    # load clinical (aka, predictive) records
    print 'Loading OncoKB clinical TSV'
    # then use it to harvest from oncokb actionable
    v = pd.read_csv(clinv, sep='\t')
    v = v[v['Gene'].isin(genes)]
    cols = {
        'Gene': 'gene',
        'Alteration': 'variant',
        'Cancer Type': 'cancerType',
        'Level': 'level',
        'Drugs(s)': 'drug',
        'PMIDs for drug': 'drugPmids',
        'Abstracts for drug': 'drugAbstracts'
    }
    v = v.rename(columns=cols)
    v = v.fillna('')
    for idx, row in v.iterrows():
        r = requests.get(
            'http://oncokb.org/api/v1/variants/lookup?hugoSymbol={}&variant={}'.format(row['gene'], row['variant'])
        )
        for ret in r.json():
            if str(ret['name']) == v['variant'][idx]:
                v.at[idx, 'variant'] = ret


    # load biological (aka, predisposing) records
    print 'Loading OncoKB biological TSV'
    # then use it to harvest from oncokb biologic
    b = pd.read_csv(biov, sep='\t')
    b = b[b['Gene'].isin(genes)]
    cols = {
        'Gene': 'gene',
        'Alteration': 'variant',
        'Oncogenicity': 'oncogenic',
        'Mutation Effect': 'mutationEffect',
        'PMIDs for Mutation Effect': 'mutationEffectPmids',
        'Abstracts for Mutation Effect': 'mutationEffectAbstracts'
    }
    b = b.rename(columns=cols)
    b = b.fillna('')
    for idx, row in b.iterrows():
        FLAG = False
        r = requests.get(
            'http://oncokb.org/api/v1/variants/lookup?hugoSymbol={}&variant={}'.format(row['gene'], row['variant'])
        )
        for ret in r.json():
            # if we find a matching variant in the api, use that instead
            if str(ret['name']) == b['variant'][idx]:
                b.at[idx, 'variant'] = ret
                FLAG = True
        if not FLAG:
            # we didn't find a match, so
            check = row['variant'].replace('?', ' ').split()
            for var in variants:
                i = 0
                if var['gene']['hugoSymbol'] == row['gene']:
                    for bits in check:
                        if bits in var['name']:
                            i = i + 1
                    if i == len(check):
                        b.at[idx, 'variant'] = var

    for gene in genes:
        gene_data = {'gene': gene, 'oncokb': {}}
        gene_data['oncokb']['clinical'] = v[v['gene'].isin([gene])].to_dict(orient='records')
        for clinical in gene_data['oncokb']['clinical']:
            key = "LEVEL_{}".format(clinical['level'])
            if key in levels:
                clinical['level_label'] = levels[key]
            else:
                print '{} not found'.format(clinical['level'])
        gene_data['oncokb']['biological'] = b[b['gene'].isin([gene])].to_dict(orient='records')
        yield gene_data


def _enrich_feature(gene, alteration, feature):
    matches = LOOKUP_TABLE.get_entries(gene, alteration)

    if len(matches) > 0:
        # FIXME: just using the first match for now;
        # it's not clear what to do if there are multiple matches.
        match = matches[0]
        feature['chromosome'] = str_or_none(match['chrom'])
        feature['start'] = match['start']
        feature['end'] = match['end']
        feature['ref'] = match['ref']
        feature['alt'] = match['alt']
        feature['referenceName'] = str_or_none(match['build'])
    else:
        # attempt to use the gene info to get some info, e.g. the chromosome, at least
        gene_meta = get_gene(gene)
        feature['chromosome'] = gene_meta[0]['chromosome']

    return feature


def convert(gene_data):
    gene = gene_data['gene']
    oncokb = {'clinical': [], 'biological': []}

    if 'oncokb' in gene_data:
        oncokb = gene_data['oncokb']

    # this section yields a feature association for the predictive evidence
    for clinical in oncokb['clinical']:
        variant = clinical['variant']
        alteration = variant['alteration']
        gene_data = variant['gene']

        # if feature is 'Oncogenic Mutations' then merge in biological
        features = []
        if variant['name'] == 'Oncogenic Mutations':
            for biological in oncokb['biological']:
                if biological['oncogenic'] in ['Likely Oncogenic', 'Oncogenic']:
                    var = biological['variant']
                    gene_data = var['gene']
                    alteration = var['alteration']
                    feature = {
                        'geneSymbol': gene,
                        'description': '{} {}'.format(gene, var['name']),
                        'name': var['name'],
                        'entrez_id': gene_data['entrezGeneId'],
                        'refseq': str_or_none(var['gene']['curatedRefSeq']),
                        'isoform': str_or_none(var['gene']['curatedIsoform']),
                        'biomarker_type': variant['consequence']['term']
                    }
                    feature = _enrich_feature(gene, alteration, feature)
                    features.append(feature)

        feature = {
            'geneSymbol': gene,
            'name': variant['name'],
            'description': '{} {}'.format(gene, variant['name']),
            'entrez_id': gene_data['entrezGeneId'],
            'refseq': str_or_none(variant['gene']['curatedRefSeq']),
            'isoform': str_or_none(variant['gene']['curatedIsoform']),
            'biomarker_type': variant['consequence']['term']
        }
        feature = _enrich_feature(gene, alteration, feature)
        features.append(feature)

        association = {
            'description': clinical['level_label'],
            'variant_name': variant['name'],
            'environmentalContexts': []
        }
        for drug in clinical['drug'].split(', '):
            association['environmentalContexts'].append({'description': drug})
        association['phenotypes'] = [{ 'description': clinical['cancerType'] }]

        # grab all publications from abstracts or PMIDs for piblication list
        abstract = []
        if clinical['drugAbstracts'] != '':
            absts = clinical['drugAbstracts'].split('; ')
            for i in range(len(absts)):
                print absts[i]
                abstract.append({'text': absts[i], 'link': ''})
                for bit in abstract[i]['text'].split():
                    if 'http' in bit:
                        abstract[i]['link'] = bit
        if clinical['drugPmids'] != '':
            pmids = clinical['drugPmids'].split(', ')
            for i in range(len(pmids)):
                abstract.append({'text': '', 'link': 'http://www.ncbi.nlm.nih.gov/pubmed/{}'.format(pmids[i])})
        clinical['drugAbstracts'] = abstract

        association['evidence'] = [{
            "evidenceType": {
                "sourceName": "oncokb",
                "id": '{}-{}'.format(gene,
                                     clinical['cancerType'])
            },
            'type': 'Predictive',  # all clinical records regard therapy, so they're all Predictive
            'description': clinical['level'],
            'info': {
                'publications':
                    [drugAbstract['link']
                        for drugAbstract in clinical['drugAbstracts']]
            }
        }]
        # add summary fields for Display
        association['source_link'] = 'http://oncokb.org/#/gene/{}/variant/{}'.format(gene, quote_plus(variant['name']))
        association = el.evidence_label(clinical['level'], association, na=True)
        association = ed.evidence_direction(clinical['level_label'], association, na=True)

        if len(clinical['drugAbstracts']) > 0:
            association['publication_url'] = clinical['drugAbstracts'][0]['link']  # NOQA
        else:
            for drugPmid in clinical['drugPmids'].split(', '):
                association['publication_url'] = 'http://www.ncbi.nlm.nih.gov/pubmed/{}'.format(drugPmid)  # NOQA
                break

        association['drug_labels'] = ','.join([drug for drug in clinical['drug']])   # NOQA
        feature_names = ', '.join(['{}:{}'.format(f["geneSymbol"], f["name"]) for f in features])

        source_url = None
        if len(features) > 0:
            f = features[0]
            source_url = 'http://oncokb.org/#/gene/{}/variant/{}'.format(f["geneSymbol"], f["name"])

        feature_association = {
            'genes': [gene],
            'features': features,
            'feature_names': feature_names,
            'association': association,
            'source': 'oncokb',
            'source_url': source_url,
            'oncokb': {'clinical': clinical}
        }
        yield feature_association

    # this section yields a feature association for the oncogenic evidence
    for biological in oncokb['biological']:
        variant = biological['variant']
        gene_data = variant['gene']
        alteration = variant['alteration']

        feature = {
            'geneSymbol': gene,
            'name': variant['name'],
            'description': '{} {}'.format(gene.encode('utf8'), variant['name'].encode('utf8')),
            'entrez_id': gene_data['entrezGeneId'],
            'refseq': str_or_none(variant['gene']['curatedRefSeq']),
            'isoform': str_or_none(variant['gene']['curatedIsoform']),
            'biomarker_type': variant['consequence']['term']
        }

        # Look up variant and add position information.
        feature = _enrich_feature(gene, alteration, feature)

        association = {
            'variant_name': variant['name'],
            'description': variant['consequence']['description'],
            'environmentalContexts': []
        }

        if biological['oncogenic'] in ['Likely Oncogenic', 'Oncogenic']:
            association['phenotypes'] = [{'description': 'cancer'}]

        association['evidence'] = [{
            "evidenceType": {
                "sourceName": "oncokb",
                "id": '{}-{}'.format(gene,
                                     biological['oncogenic'])
            },
            'type': 'Predisposing',  # all biological records regard pathogenicity, so they're all Predisposing
            'description': biological['mutationEffect'],
            'info': {
                'publications':
                    ['http://www.ncbi.nlm.nih.gov/pubmed/{}'.format(Pmid)
                        for Pmid in biological['mutationEffectPmids'].split(', ')]
            }
        }]

        association['oncogenic'] = biological['oncogenic']
        association['evidence_label'] = None

        if len(biological['mutationEffectPmids']) > 0:
            for drugPmid in biological['mutationEffectPmids'].split(', '):
                association['publication_url'] = 'http://www.ncbi.nlm.nih.gov/pubmed/{}'.format(drugPmid)  # NOQA
                break

        feature_names = feature["geneSymbol"] + ' ' + feature["name"]

        source_url = 'http://oncokb.org/#/gene/{}/variant/{}'.format(
            feature["geneSymbol"].encode('utf8'), feature["name"].encode('utf8'))

        feature_association = {
            'genes': [gene],
            'features': [feature],
            'feature_names': feature_names,
            'association': association,
            'source': 'oncokb',
            'source_url': source_url,
            'oncokb': {'biological': biological}
        }

        yield feature_association


def harvest_and_convert(genes):
    """ get data from oncokb, convert it to ga4gh and return via yield """
    for gene_data in harvest(genes):
        # print "harvester_yield {}".format(jax_evidence.keys())
        for feature_association in convert(gene_data):
            # print "convert_yield {}".format(feature_association.keys())
            yield feature_association


def _test():
    for feature_association in harvest_and_convert(None):
        print feature_association.keys()
        break

if __name__ == '__main__':
    _test()
