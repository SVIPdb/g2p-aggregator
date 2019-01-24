#!/usr/bin/python

import requests
import copy
from lookups import evidence_label as el, evidence_direction as ed
import logging
from tqdm import tqdm


# TODO: implement clinvar harvester (likely from a file instead of an API; the below is a copy of civic)

def harvest(genes):
    """ given an array of gene symbols, harvest them from clinvar"""
    # harvest all genes
    if not genes:
        r = requests.get('https://clinvar_url_here/api/genes?count=99999')
        for record in r.json()['records']:
            variants = record['variants']
            gene = record['name']
            variants_details = []
            for variant in tqdm(variants, desc="fetching clinvar variant data for %s" % gene):
                r = requests.get('https://clinvar_url_here/api/variants/{}'.format(variant['id']))
                variants_details.append(r.json())
            gene_data = {'gene': gene, 'clinvar': {'variants': variants_details}}
            yield gene_data
    else:
        # harvest some genes
        for gene in set(genes):
            r = requests.get('https://clinvar_url_here/api/genes/{}?identifier_type=entrez_symbol'.format(gene))
            if r.status_code != 200 or len(r.json()['variants']) == 0:
                # print "{} Found no variants in clinvar".format(gene)
                gene_data = {'gene': gene, 'clinvar': {}}
            else:
                variants = r.json()['variants']
                variants_details = []
                for variant in tqdm(variants, desc="fetching clinvar variant data for %s" % gene):
                    r = requests.get('https://clinvar_url_here/api/variants/{}'.format(variant['id']))
                    variants_details.append(r.json())
                gene_data = {'gene': gene,
                             'clinvar': {'variants': variants_details}}
            yield gene_data


def convert(gene_data):
    """ given gene data from clinvar, convert it to ga4gh """
    try:
        variants = gene_data['clinvar']['variants']
        for variant in variants:
            # feature = {
            #     'geneSymbol': variant['entrez_name'],
            #     'entrez_id': variant['entrez_id'],
            #     'start': variant['coordinates']['start'],
            #     'end': variant['coordinates']['stop'],
            #     'referenceName': str(variant['coordinates']['reference_build']),
            #     'refseq': None,
            #     'isoform': str(variant['coordinates']['representative_transcript']),
            #     'chromosome': str(variant['coordinates']['chromosome']),
            #     'ref': str(variant['coordinates']['reference_bases']),
            #     'alt': str(variant['coordinates']['variant_bases']),
            #     'name': variant['name'],
            #     'description': '{} {}'.format(variant['entrez_name'], variant['name'])
            # }
            #
            # if 'variant_types' in variant and len(variant['variant_types']) > 0:
            #     feature['biomarker_type'] = variant['variant_types'][0]['display_name']

            for evidence_item in variant['evidence_items']:
                pass
                # association = {}
                #
                # for part in variant['name'].split():
                #     if '-' not in part and not part == variant['entrez_name']:
                #         association['variant_name'] = part
                #
                # association['description'] = evidence_item['description']
                # association['environmentalContexts'] = []
                # environmentalContexts = association['environmentalContexts']

                # source_url = "https://clinvar_url_here/gene/%s/variant/%s/evidence/%s"
                #     variant['gene_id'], variant['id'], evidence_item['id'])
                #
                # feat_assoc = {
                #     'genes': [gene_data['gene']],
                #     'features': [feature],
                #     'feature_names': evidence_item['name'],
                #     'association': association,
                #     'source': 'clinvar',
                #     'source_url': source_url,
                #     'clinvar': payload
                # }
                #
                # yield feat_assoc

    except Exception as e:
        logging.error(gene_data['gene'], exc_info=1, ex=e)


def harvest_and_convert(genes):
    """ get data from clinvar, convert it to ga4gh and return via yield """
    for gene_data in harvest(genes):
        # print "harvester_yield {}".format(gene_data.keys())
        for feat_assoc in convert(gene_data):
            # print "convert_yield {}".format(feature_association.keys())
            yield feat_assoc


# main
if __name__ == '__main__':
    for feature_association in harvest_and_convert(["MDM2"]):
        logging.info(feature_association)
