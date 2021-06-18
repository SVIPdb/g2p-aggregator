#!/usr/bin/python

import requests
import copy
import logging
from tqdm import tqdm

from normalizers.gene_enricher import get_gene
from utils_ex.formatting import unicode_or_none
from utils_ex.iterables import matched
from lookups.accession_mapping import NoMatchError, ensembl_txac_to_refseq

# CIVIC_API_URL = "civic.genome.wustl.edu"
CIVIC_API_URL = "civicdb.org"
# CIVIC_API_URL = "civic.nexus.ethz.ch"


def harvest(genes):
    """ given an array of gene symbols, harvest them from civic"""

    # harvest all genes
    if not genes:
        r = requests.get('https://%s/api/genes?count=99999' % CIVIC_API_URL)
        for record in r.json()['records']:
            variants = record['variants']
            gene = record['name']
            variants_details = []
            for variant in tqdm(variants, desc="fetching civic variant data for %s" % gene):
                r = requests.get(
                    'https://{}/api/variants/{}'.format(CIVIC_API_URL, variant['id']))
                variants_details.append(r.json())
            gene_data = {'gene': gene, 'civic': {'variants': variants_details}}
            yield gene_data

    else:
        # harvest some genes
        for gene in set(genes):
            r = requests.get(
                'https://{}/api/genes/{}?identifier_type=entrez_symbol'.format(CIVIC_API_URL, gene))

            if r.status_code != 200 or len(r.json()['variants']) == 0:
                logging.info(
                    "Found no variants in CIViC for gene {}".format(gene))
                gene_data = {'gene': gene, 'civic': {'variants': []}}
            else:
                variants = r.json()['variants']
                variants_details = []

                for variant in tqdm(variants, desc="fetching civic variant data for %s" % gene):
                    v = requests.get(
                        'https://{}/api/variants/{}'.format(CIVIC_API_URL, variant['id']))
                    variants_details.append(v.json())

                gene_data = {'gene': gene, 'civic': {
                    'variants': variants_details}}

            yield gene_data


def _extract_name(variant):
    """

    :param variant:
    :return:
    """
    for part in variant['name'].split():
        if '-' not in part and not part == variant['entrez_name']:
            return part


def accession(hgvs_str):
    if not hgvs_str:
        return None
    return hgvs_str[:(hgvs_str.index(':'))]


def convert(gene_data):
    """ given gene data from civic, convert it to ga4gh """
    variants = gene_data['civic']['variants']

    # we can get some coarse location info, e.g. the chromosome, from the gene symbol itself
    # we'll retrieve that as a failover in case the civic entry is missing that info
    try:
        gene_meta = get_gene(gene_data['gene'])[0]
    except ValueError as ex:
        # this probably means the gene is missing, which means we can't really do anything...
        logging.warn(str(ex))
        return

    for variant in variants:
        # parse out hgvs strings from 'hgvs_expression' and assign each to a type
        hgvs_exprs = variant['hgvs_expressions']
        hgvs_types = {
            'hgvs_g': matched(hgvs_exprs, lambda x: x.startswith("NC_")),
            'hgvs_c': matched(hgvs_exprs, lambda x: x.startswith("NM_")),
            'hgvs_p': matched(hgvs_exprs, lambda x: x.startswith("NP_")),
            'hgvs_ensembl_c': matched(hgvs_exprs, lambda x: x.startswith("ENST")),
        }

        feature = {
            'geneSymbol': variant['entrez_name'],
            'entrez_id': variant['entrez_id'],
            'start': variant['coordinates']['start'],
            'end': variant['coordinates']['stop'],
            'referenceName': unicode_or_none(variant['coordinates']['reference_build']),
            'refseq': accession(hgvs_types['hgvs_c']),
            'isoform': unicode_or_none(variant['coordinates']['representative_transcript']),
            'chromosome': unicode_or_none(variant['coordinates']['chromosome']),
            'ref': unicode_or_none(variant['coordinates']['reference_bases']),
            'alt': unicode_or_none(variant['coordinates']['variant_bases']),
            'name': variant['name'],
            'description': u'{} {}'.format(variant['entrez_name'], variant['name']),
        }

        # also insert the hgvs strings to potentially save work for the normalizers downstream
        feature.update(hgvs_types)

        # if our feature is lacking information we can infer from the gene metadata, fill that in
        if not feature['chromosome'] and gene_meta['chromosome']:
            feature['chromosome'] = gene_meta['chromosome']

        if 'variant_types' in variant and len(variant['variant_types']) > 0:
            feature['biomarker_type'] = variant['variant_types'][0]['display_name']

        # if the referenceName (aka the assembly) is missing, we might be able to infer it from the ensembl version
        if feature['referenceName'] is None and variant['coordinates']['ensembl_version'] == 75:
            feature['referenceName'] = 'GRCh37'

        # if the refseq is still missing but isoform is specified, attempt to convert that into an NCBI refseq
        if not feature['refseq'] and feature['isoform']:
            try:
                feature['refseq'] = ensembl_txac_to_refseq(feature['isoform'])
            except NoMatchError as ex:
                logging.warn(ex)
                feature['refseq'] = None

        for evidence_item in variant['evidence_items']:
            # FIXME: maybe we should skip entries where the evidence item was rejected; see evidence_item['status']
            # example of erroneous submission:
            # https://civicdb.org/events/genes/5/summary/variants/842/summary/evidence/1941/talk/comments
            if evidence_item['status'] == 'rejected':
                logging.warn("Skipping evidence item %s because it has status %s" % (
                    evidence_item['id'], evidence_item['status']))
                continue

            evidence_url = "https://{}/events/genes/{}/summary/variants/{}/summary/evidence/{}/summary#evidence".format(
                CIVIC_API_URL,
                variant['gene_id'], variant['id'], evidence_item['id']
            )

            association = {
                'variant_name': _extract_name(variant),
                'description': evidence_item['description'],
                'environmentalContexts': [
                    {
                        'term': drug.get('name'),
                        'description': drug.get('name'),
                        'id': drug.get('pubchem_id')
                    }
                    for drug in evidence_item['drugs']
                ],
                'drug_interaction_type': evidence_item['drug_interaction_type'],
                'phenotypes': [{
                    'description': evidence_item['disease']['name'] if 'disease' in evidence_item and evidence_item['disease'] is not None and 'name' in evidence_item['disease'] else 'N/A',
                    'id': evidence_item['disease']['url'] if 'disease' in evidence_item and evidence_item['disease'] is not None and 'url' in evidence_item['disease'] else 'N/A'
                }],
                'evidence': [{
                    "evidenceType": {
                        "sourceName": "CIVIC",
                        "id": '{}'.format(evidence_item['id'])
                    },
                    'info': {
                        'publications': [
                            evidence_item['source']['source_url']
                        ]
                    }
                }],

                'evidence_type': evidence_item['evidence_type'],
                'evidence_direction': evidence_item['evidence_direction'],
                'clinical_significance': evidence_item['clinical_significance'],
                'evidence_level': evidence_item['evidence_level'],

                'source_link': evidence_url,
                'publication_url': (evidence_item['source']['source_url'],),

                'drug_labels': u', '.join([drug['name'] for drug in evidence_item['drugs']]) if len(evidence_item['drugs']) > 0 else None
            }

            # create snapshot of original data
            v = copy.deepcopy(variant)
            del v['evidence_items']
            v['evidence_items'] = [evidence_item]

            variant_url = 'https://{}/events/genes/{}/summary/variants/{}/summary'.format(
                CIVIC_API_URL,
                variant['gene_id'], variant['id']
            )

            feat_assoc = {
                'genes': [gene_data['gene']],
                'features': [feature],
                'feature_names': evidence_item['name'],
                'association': association,
                'source': 'civic',
                'source_url': variant_url,
                'civic': v
            }

            yield feat_assoc


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
