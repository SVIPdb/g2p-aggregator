#!/usr/bin/python
import os
import re
import logging
import urllib2

from pathlib import Path
import pandas as pd
import requests

from tqdm import tqdm
from retrying import retry

from lookups import cosmic_lookup_table
from normalizers.gene_enricher import get_gene
from utils_ex.downloading import acquire_files
from utils_ex.formatting import unicode_or_none


# ----------------------------------------------------------------------------------------------------------------
# -- file definitions, updating
# ----------------------------------------------------------------------------------------------------------------

# OncoKB harvester now pulls from the below downloadable OncoKB files
# and supplements with additional variant data pulled from their public
# API. This is because pulling from the private API gives unpredictable
# results and their is no endpoint in the public API that gives the
# same drug-gene-variant association information as was being
# pulled from the private API.

data_paths = acquire_files({
    'oncokb_allActionableVariants.txt': {
        'path': '../data/oncokb_allActionableVariants.txt',
        'url': 'http://oncokb.org/api/v1/utils/allActionableVariants.txt'
    },
    'oncokb_allAnnotatedVariants.txt': {
        'path': '../data/oncokb_allAnnotatedVariants.txt',
        'url': 'http://oncokb.org/api/v1/utils/allAnnotatedVariants.txt'
    }
})

clinv = data_paths['oncokb_allActionableVariants.txt']
biov = data_paths['oncokb_allAnnotatedVariants.txt']


# used to get COSMIC info about genes/alterations
LOOKUP_TABLE = cosmic_lookup_table.CosmicLookup("../data/cosmic_lookup_table.tsv")

# as of dec 2019, oncokb now suggests using an API key
# NOTE (jan 30, 2020): it seems the v1 api is still accessible without it, so this is just a warning, not a requirement
ONCOKB_API_KEY = os.environ.get('ONCOKB_API_KEY')


@retry(wait_exponential_multiplier=1000, wait_exponential_max=10000, stop_max_attempt_number=5)
def get_oncokb_variant(session, gene, variant):
    return session.get(
        'http://oncokb.org/api/v1/variants/lookup?hugoSymbol={}&variant={}'.format(gene, variant)
    )


# ----------------------------------------------------------------------------------------------------------------
# -- harvesting
# ----------------------------------------------------------------------------------------------------------------

def harvest(genes):
    # create a session to send authorized requests
    session = requests.Session()
    if ONCOKB_API_KEY:
        session.headers.update({
            'Authorization': 'Bearer %s' % ONCOKB_API_KEY
        })
    else:
        logging.info("ONCOKB_API_KEY not set; currently this is fine, but an API key will eventually be required")

    levels = session.get('http://oncokb.org/api/v1/levels').json()

    # always get all the gene metadata, since we're going to use this later to reconstruct the variant object's gene key
    # index it by hugoSymbol, which i've verified to be unique
    all_genes = {x['hugoSymbol']: x for x in session.get('http://oncokb.org/api/v1/genes').json()}

    if not genes:
        genes = set(all_genes.keys())

    # get all variants
    print 'Gathering all OncoKB variants'
    variants = session.get('http://oncokb.org/api/v1/variants').json()

    # load clinical (aka, predictive) records
    print 'Loading OncoKB clinical TSV'

    # then use it to harvest from oncokb actionable
    v = pd.read_csv(clinv, sep='\t')
    v = v[v['Hugo Symbol'].isin(genes)]
    cols = {
        'Hugo Symbol': 'gene',
        'Alteration': 'variant',
        'Cancer Type': 'cancerType',
        'Level': 'level',
        'Drugs(s)': 'drug',
        'PMIDs for drug': 'drugPmids',
        'Abstracts for drug': 'drugAbstracts'
    }
    v = v.rename(columns=cols)
    v = v.fillna('')

    # FIXME: the API returns a dict for 'variant' (consisting of mostly redundant info), whereas the file
    #  oncokb_allActionableVariants.txt simply contains the protein change. ideally we should normalize
    #  the contents of variant somewhere, since downstream stuff is expecting a dict...
    # FIXME: i'm also unclear on why we need to issue special requests per variant when we get all the variants
    #  above, in the call to http://oncokb.org/api/v1/variants...

    for idx, row in tqdm(v.iterrows(), total=len(v.index), desc="oncokb clinical lookups"):
        r = get_oncokb_variant(session, row['gene'], row['variant'])
        matched = False

        # attempt to find a match in the responses for this gene
        for ret in r.json():
            if unicode(ret['name']) == v['variant'][idx]:
                v.at[idx, 'variant'] = ret
                matched = True

        if matched is None:
            # we didn't get a match, so we need to synthesize the structure we'd normally get from oncokb
            # first, parse the variant's name, then use the gene data from before to populate the gene key
            m = re.match(r'^(?P<ref>[A-Z]+)(?P<pos>[0-9]+)(_(?P<pos2>[0-9]+))?(?P<alt>[A-Z]+)?$', v['variant'][idx])

            if not m:
                # we can't do much if anything if we couldn't parse it, so move along
                continue

            groups = m.groupdict()

            v.at[idx, 'variant'] = {
                u'variantResidues': groups.get('alt'),  # will be None in the case of a deletion
                u'proteinStart': groups['pos'],  # 601
                u'name': v['variant'][idx],  # K601
                u'proteinEnd': groups.get('pos2', groups['pos']),  # 601
                u'refResidues': groups['ref'],
                u'alteration': v['variant'][idx],  # K601
                # FIXME: should we return a consequence at all?
                u'consequence': {
                    u'term': u'NA', u'description': u'NA', u'isGenerallyTruncating': False
                },
                u'gene': all_genes[v['gene'][idx]]
            }


    # load biological (aka, predisposing) records
    print 'Loading OncoKB biological TSV'

    # then use it to harvest from oncokb biologic
    b = pd.read_csv(biov, sep='\t')
    b = b[b['Hugo Symbol'].isin(genes)]
    cols = {
        'Hugo Symbol': 'gene',
        'Alteration': 'variant',
        'Oncogenicity': 'oncogenic',
        'Mutation Effect': 'mutationEffect',
        'PMIDs for Mutation Effect': 'mutationEffectPmids',
        'Abstracts for Mutation Effect': 'mutationEffectAbstracts'
    }
    b = b.rename(columns=cols)
    b = b.fillna('')

    for idx, row in tqdm(b.iterrows(), total=len(b.index), desc="oncokb biological lookups"):
        FLAG = False
        r = get_oncokb_variant(session, row['gene'], row['variant'])

        for ret in r.json():
            # if we find a matching variant in the api, use that instead
            if unicode(ret['name']) == b['variant'][idx]:
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
        feature['chromosome'] = unicode_or_none(match['chrom'])
        feature['start'] = match['start']
        feature['end'] = match['end']
        feature['ref'] = match['ref']
        feature['alt'] = match['alt']
        feature['referenceName'] = unicode_or_none(match['build'])
    else:
        # attempt to use the gene info to get some info, e.g. the chromosome, at least
        gene_meta = get_gene(gene)
        feature['chromosome'] = gene_meta[0]['chromosome']

    return feature


def create_feature(variant):
    gene_data = variant['gene']
    gene = gene_data['hugoSymbol']

    feature = {
        'geneSymbol': gene,
        'description': u'{} {}'.format(gene, variant['name']),
        'name': variant['name'].replace(u'\u2013', '-'),
        'entrez_id': gene_data['entrezGeneId'],
        'refseq': unicode_or_none(variant['gene']['curatedRefSeq']),
        'isoform': unicode_or_none(variant['gene']['curatedIsoform']),
        'biomarker_type': variant['consequence']['term']
    }
    feature['source_link'] = u'http://oncokb.org/#/gene/{}/variant/{}'.format(feature["geneSymbol"], feature["name"])

    return _enrich_feature(gene, variant['alteration'], feature)


def convert(gene_data):
    gene = gene_data['gene']
    oncokb = {'clinical': [], 'biological': []}

    if 'oncokb' in gene_data:
        oncokb = gene_data['oncokb']

    # this section yields a feature association for the predictive evidence
    for clinical in oncokb['clinical']:
        variant = clinical['variant']

        # we may not have been able to create a variant structure earlier; if that's the case, variant
        # will just be a string and there isn't much point in passing it along, so we'll drop it with a warning
        if type(variant) is not dict:
            logging.warn("OncoKB convert(): received non-dict variant %s, dropping it" % variant)
            continue

        # if feature is 'Oncogenic Mutations' then merge in biological
        features = []
        if variant['name'] == 'Oncogenic Mutations':
            for biological in oncokb['biological']:
                if biological['oncogenic'] in ['Likely Oncogenic', 'Oncogenic']:
                    features.append(create_feature(biological['variant']))

        # create the main feature for this clinical entry, too, and append it
        feature = create_feature(variant)
        features.append(feature)

        clinical_level = clinical['level'].strip("LEVEL_")

        association = {
            'description': clinical['level_label'],
            'variant_name': variant['name'],
            'environmentalContexts': [],

            'evidence_type': 'Predictive',
            'evidence_direction': None,  # FIXME: does oncokb provide this? do we presume it supports it?
            'clinical_significance': 'Resistance' if clinical_level.startswith("R") else 'Sensitivity/Response',
            'evidence_level': clinical_level,
        }
        for drug in clinical['drug'].split(', '):
            association['environmentalContexts'].append({'description': drug})
        association['phenotypes'] = [{ 'description': clinical['cancerType'] }]

        # grab all publications from abstracts or PMIDs for piblication list
        abstract = []
        if clinical['drugAbstracts'] != '':
            absts = clinical['drugAbstracts'].split('; ')
            for i in range(len(absts)):
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

        association['source_link'] = feature['source_link']

        # establish evidence label, level, clinical_significance here
        # FIXME: we disabled these because they're borked, but is that the right choice?
        # association = el.evidence_label(clinical['level'], association, na=True)
        # association['clinical_significance'] = ed.evidence_direction(clinical['level_label'], na=True)

        if len(clinical['drugAbstracts']) > 0:
            association['publication_url'] = clinical['drugAbstracts'][0]['link']  # NOQA
        else:
            for drugPmid in clinical['drugPmids'].split(', '):
                association['publication_url'] = u'http://www.ncbi.nlm.nih.gov/pubmed/{}'.format(drugPmid)  # NOQA
                break

        association['drug_labels'] = clinical['drug']  # NOQA
        feature_names = u', '.join(['{}:{}'.format(f["geneSymbol"], f["name"]) for f in features])

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

        try:
            feature = create_feature(variant)
        except ValueError as ex:
            # this probably means the gene is missing, which means we can't really do anything...
            logging.warn(str(ex))
            continue

        association = {
            'variant_name': variant['name'],
            'description': variant['consequence']['description'],
            'environmentalContexts': [],

            'evidence': [{
                "evidenceType": {
                    "sourceName": "oncokb",
                    "id": '{}-{}'.format(gene,
                                         biological['oncogenic'])
                },
                'info': {
                    'publications':
                        ['http://www.ncbi.nlm.nih.gov/pubmed/{}'.format(Pmid)
                            for Pmid in biological['mutationEffectPmids'].split(', ')]
                }
            }],
            'oncogenic': biological['oncogenic'],
            'source_link': feature['source_link'],

            # FIXME: why do we only assign a phenotype of 'cancer' if it's oncogenic?
            #  just made the executive decision to disable the check and make predisposing entries always about cancer;
            #  verify if that's ok later. old code below for posterity.
            # if biological['oncogenic'] in ['Likely Oncogenic', 'Oncogenic']:
            #   association['phenotypes'] = [{'description': 'cancer'}]
            'phenotypes': [{'description': 'cancer'}],

            'evidence_type': 'Predisposing',  # all biological records regard pathogenicity, so they're all Predisposing
            'evidence_direction': None,  # FIXME: does oncokb provide this? do we presume it supports it?
            'clinical_significance': '%s, %s' % (biological['oncogenic'], biological['mutationEffect']),
            'evidence_level': None,
        }

        if len(biological['mutationEffectPmids']) > 0:
            for drugPmid in biological['mutationEffectPmids'].split(', '):
                association['publication_url'] = 'http://www.ncbi.nlm.nih.gov/pubmed/{}'.format(drugPmid)  # NOQA
                break

        # FIXME: is there a reason that this is a duplicate of source_link, but is generated in a different way
        source_url = 'http://oncokb.org/#/gene/{}/variant/{}'.format(
            feature["geneSymbol"].encode('utf8'), feature["name"].encode('utf8'))

        feature_association = {
            'genes': [gene],
            'features': [feature],
            'feature_names': feature["geneSymbol"] + ' ' + feature["name"],
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
