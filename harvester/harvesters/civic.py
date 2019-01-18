#!/usr/bin/python

import requests
import copy
from lookups import evidence_label as el, evidence_direction as ed
import logging
from tqdm import tqdm

from utils import str_or_none


def harvest(genes):
    """ given an array of gene symbols, harvest them from civic"""
    # harvest all genes
    if not genes:
        r = requests.get('https://civic.genome.wustl.edu/api/genes?count=99999')  # NOQA
        for record in r.json()['records']:
            variants = record['variants']
            gene = record['name']
            variants_details = []
            for variant in tqdm(variants, desc="fetching civic variant data for %s" % gene):
                r = requests.get('https://civic.genome.wustl.edu/api/variants/{}'.format(variant['id']))   # NOQA
                variants_details.append(r.json())
            gene_data = {'gene': gene, 'civic': {'variants': variants_details}}
            yield gene_data
    else:
        # harvest some genes
        for gene in set(genes):
            r = requests.get('https://civic.genome.wustl.edu/api/genes/{}?identifier_type=entrez_symbol'.format(gene))  # NOQA
            if r.status_code != 200 or len(r.json()['variants']) == 0:
                # print "{} Found no variants in civic".format(gene)
                gene_data = {'gene': gene, 'civic': {}}
            else:
                variants = r.json()['variants']
                variants_details = []
                for variant in tqdm(variants, desc="fetching civic variant data for %s" % gene):
                    r = requests.get('https://civic.genome.wustl.edu/api/variants/{}'.format(variant['id']))   # NOQA
                    variants_details.append(r.json())
                gene_data = {'gene': gene,
                             'civic': {'variants': variants_details}}
            yield gene_data


def convert(gene_data):
    """ given gene data from civic, convert it to ga4gh """
    try:
        variants = gene_data['civic']['variants']
        for variant in variants:
            feature = {
                'geneSymbol': variant['entrez_name'],
                'entrez_id': variant['entrez_id'],
                'start': variant['coordinates']['start'],
                'end': variant['coordinates']['stop'],
                'referenceName': str_or_none(variant['coordinates']['reference_build']),
                'refseq': None,
                'isoform': str_or_none(variant['coordinates']['representative_transcript']),
                'chromosome': str_or_none(variant['coordinates']['chromosome']),
                'ref': str_or_none(variant['coordinates']['reference_bases']),
                'alt': str_or_none(variant['coordinates']['variant_bases']),
                'name': variant['name'],
                'description': '{} {}'.format(variant['entrez_name'], variant['name'])
            }

            # debugging
            # if variant['name'] == 'N842S':
            #     print "****: %s" % str(variant)

            if 'variant_types' in variant and len(variant['variant_types']) > 0:
                feature['biomarker_type'] = variant['variant_types'][0]['display_name']

            for evidence_item in variant['evidence_items']:
                association = {}

                for part in variant['name'].split():
                    if '-' not in part and not part == variant['entrez_name']:
                        association['variant_name'] = part

                association['description'] = evidence_item['description']
                association['environmentalContexts'] = []
                environmentalContexts = association['environmentalContexts']

                for drug in evidence_item['drugs']:
                    environmentalContexts.append({
                        'term': drug['name'],
                        'description': drug['name'],
                        'id': drug['pubchem_id']
                    })

                association['phenotypes'] = [{
                    'description': evidence_item['disease']['name'],
                    'id': evidence_item['disease']['url']
                }]

                association['evidence'] = [{
                    "evidenceType": {
                        "sourceName": "CIVIC",
                        "id": '{}'.format(evidence_item['id'])
                    },
                    'type': evidence_item.get('evidence_type'),
                    'description': evidence_item['clinical_significance'],
                    'info': {
                        'publications': [
                            evidence_item['source']['source_url']
                        ]
                    }
                }]

                # add summary fields for Display
                association = el.evidence_label(
                    evidence_item['evidence_level'], association, na=True
                )
                association = ed.evidence_direction(
                    evidence_item['clinical_significance'], association
                )

                association[
                    'source_link'] = 'https://civic.genome.wustl.edu/events/genes/{}/summary/variants/{}/summary'.format(
                    variant['gene_id'], variant['id'])  # NOQA
                association['publication_url'] = evidence_item['source']['source_url'],  # NOQA
                if len(evidence_item['drugs']) > 0:
                    association['drug_labels'] = ','.join([drug['name'] for drug in evidence_item['drugs']])  # NOQA
                # create snapshot of original data
                v = copy.deepcopy(variant)
                del v['evidence_items']
                v['evidence_items'] = [evidence_item]
                source_url = "https://civicdb.org/events/genes/{}/summary/variants/{}/summary/evidence/{}/summary#evidence".format(
                    variant['gene_id'], variant['id'], evidence_item['id'])  # NOQA

                feat_assoc = {
                    'genes': [gene_data['gene']],
                    'features': [feature],
                    'feature_names': evidence_item['name'],
                    'association': association,
                    'source': 'civic',
                    'source_url': source_url,
                    'civic': v
                }

                yield feat_assoc

    except Exception as e:
        logging.error(gene_data['gene'], exc_info=1, ex=e)


def harvest_and_convert(genes):
    """ get data from civic, convert it to ga4gh and return via yield """
    for gene_data in harvest(genes):
        # print "harvester_yield {}".format(gene_data.keys())
        for feat_assoc in convert(gene_data):
            # print "convert_yield {}".format(feature_association.keys())
            yield feat_assoc


# main
if __name__ == '__main__':
    for feature_association in harvest_and_convert(["MDM2"]):
        logging.info(feature_association)
