#!/usr/bin/python
import csv
import itertools
import json
import operator
import os
import re

import logging
from collections import defaultdict

import tqdm

# hgvs stuff
import hgvs.normalizer
import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.assemblymapper


# globals that we use for parsing variant representations
from cosmic_lookup_table import parse_hgvc_c
from lookups.accession_mapping import ensembl_txac_to_refseq
from utils_ex.formatting import capitalize_words
from utils_ex.tqdm_logger import std_out_err_redirect_tqdm

coord_matcher = re.compile(r'(?P<chrom>[^:]+):(?P<start>[0-9]+)-(?P<end>[0-9]+)')

# these shared assembly mappers will allow us to convert HGVS g. variants to c. and p. later on
hdp = hgvs.dataproviders.uta.connect()
hgnorm = hgvs.normalizer.Normalizer(hdp)
hgvsparser = hgvs.parser.Parser()
am = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh37', normalize=True)

GENE_COUNT_CACHE_FILE = './_cached_cosmic_genes.json'
COSMIC_MUTANTS_SORTED_FILE = './CosmicMutantExport_sorted.tsv'


# ==========================================================================================
# === harvester implementation
# ==========================================================================================

def harvest(genes=None):
    with open(COSMIC_MUTANTS_SORTED_FILE) as cme_fp:
        cme_reader = csv.DictReader(cme_fp, dialect='excel-tab')

        it = itertools.groupby(cme_reader, operator.itemgetter('Gene name'))
        if genes is not None:
            it = itertools.ifilter(lambda x: x[0] in genes, it)

        for gene_symbol, samples in it:
            sample = samples.next()  # peek at the first sample in the list...

            if not sample:
                # no idea how this would happen, but if there aren't any samples for this gene, skip the gene
                continue

            full_samples = itertools.chain([sample], samples)  # ...but put it back when we iterate later

            # all genes have the same accession number, so we can do it on the gene level here
            refseq_ac = ensembl_txac_to_refseq(sample['Accession Number'])

            yield {
                'gene_symbol': gene_symbol,
                'refseq_ac': refseq_ac,
                'samples': full_samples
            }


# noinspection PyTypeChecker
def convert(gene_data, tq):
    # try:

    # NOTE: we used to grab the first variant in the Mutation AA group, but later noticed that not every
    # sample under a 'Mutation AA' entry had the same fields for, say, "Mutation CDS" (which in retrospect
    # makes perfect sense.) anyway, we went back to processing each sample individually instead of trying to
    # group them under the same AA change, which is incidentally faster due to not having to sort the samples anymore.

    for sample in gene_data['samples']:
        if (
            not sample['Mutation AA'] or not sample['Mutation CDS']
            or '?' in sample['Mutation AA'] or '?' in sample['Mutation CDS']
            # or '_' in sample['Mutation CDS']
        ):
            # afaik, we can't parse variants that have a question mark in their cDNA representation
            # logging.warn("Skipping COSMIC variant %s (AA: %s) b/c of missing cDNA rep: %s" % (
            #     sample['Mutation ID'], sample['Mutation AA'], mapped_cds))
            tq.update()  # move the counter ahead anyway
            continue

        # the following is shared by all samples under this variant, although we emit an evidence item
        # for each sample
        match = coord_matcher.match(sample['Mutation genome position'])
        coords = match.groupdict() if match else {}

        # this is way faster and gets us basically the same thing as using the hgvs library, since we don't even
        # have the right accession anyway...
        pos_info = parse_hgvc_c(sample['Mutation CDS'])

        feature = {
            'geneSymbol': gene_data['gene_symbol'],
            'entrez_id': sample['HGNC ID'],
            'start': coords.get('start'),
            'end': coords.get('stop'),
            'referenceName': sample['GRCh'],
            'refseq': gene_data['refseq_ac'],
            'isoform': sample['Accession Number'],
            'chromosome': coords.get('chrom'),
            'ref': pos_info['ref'] if pos_info else None,
            'alt': pos_info['alt'] if pos_info else None,
            'name': sample['Mutation AA'][2:],
            'description': "%s %s" % (gene_data['gene_symbol'], sample['Mutation AA'][2:]),
            'biomarker_type': sample['Mutation Description']
        }

        association = {
            'variant_name': feature['name'],
            'description': 'n/a',  # FIXME: what's a good description for this?

            'extras': json.dumps({
                'fathmm_prediction': sample['FATHMM prediction'],
                'fathmm_score': sample['FATHMM score'],
            }),

            # for COSMIC, the tissue in which the sample was found
            'environmentalContexts': [{
                'description': sample["Primary site"],
                'type': 'tissue'
            }],

            'evidence': [{
                "evidenceType": {
                    "sourceName": "cosmic",
                    "id": sample['ID_sample']
                },
                # all biological records regard pathogenicity, so they're all Predisposing
                'type': 'Predisposing',
                'description': "%s (%s)" % (
                    capitalize_words(sample['FATHMM prediction']),
                    sample['FATHMM score']
                ),
                'info': {
                    'publications': ['http://www.ncbi.nlm.nih.gov/pubmed/%s' % sample['Pubmed_PMID']]
                }
            }],

            'source_link': 'https://cancer.sanger.ac.uk/cosmic/sample/overview?id=%s' % sample['ID_sample'],

            'phenotypes': [{
                'description': sample['Primary histology'].replace("_", " ").replace("neoplasm", "cancer")
            }],

            # FIXME: figure out what these should be for COSMIC
            'evidence_type': None,
            'evidence_direction': None,
            'clinical_significance': None,
            'evidence_level': None,
        }

        # FIXME: this should really be named 'variant URL' or somesuch
        # remove the COSM prefix on the entry's ID when constructing the URL
        source_url = 'https://cancer.sanger.ac.uk/cosmic/mutation/overview?id=%s' % sample['Mutation ID'][4:],

        feat_assoc = {
            'genes': [gene_data['gene_symbol']],  # there's only ever one gene/variant
            'features': [feature],
            'feature_names': feature["geneSymbol"] + ' ' + feature["name"],
            'association': association,
            'source': 'cosmic',
            'source_url': source_url,
            'cosmic': sample
        }

        tq.update()

        yield feat_assoc

    # except Exception as e:
    #     raise e
    #     # logging.exception(gene_data['gene_symbol'])


def harvest_and_convert(genes):
    gene_counts = _get_gene_counts()
    sample_count = sum(gene_counts[g] for g in genes) if genes else gene_counts['_total']

    with std_out_err_redirect_tqdm() as orig_stdout:
        with tqdm.tqdm(total=sample_count, desc="harvesting %s" % genes, file=orig_stdout, dynamic_ncols=True) as tq:
            for gene_data in harvest(genes):
                for feat_assoc in convert(gene_data, tq):
                    yield feat_assoc


def _get_gene_counts():
    # returns a dict of sample counts per gene and a special key, _total, which is the total number of samples
    # (the dict is computed from COSMIC_MUTANTS_SORTED_FILE the first time, then returned from a cache afterward)

    if os.path.isfile(GENE_COUNT_CACHE_FILE) and \
       os.path.getmtime(GENE_COUNT_CACHE_FILE) >= os.path.getmtime(COSMIC_MUTANTS_SORTED_FILE):
        # a cachefile with a more recent modification date exists, so send that back
        with open(GENE_COUNT_CACHE_FILE, 'r') as fp:
            gene_counts = json.load(fp)
    else:
        # we need to re-count the variants from the source file and recreate/refresh the cache
        gene_counts = defaultdict(int)
        total_samples = 0

        with open(COSMIC_MUTANTS_SORTED_FILE) as cme_fp:
            cme_reader = csv.DictReader(cme_fp, dialect='excel-tab')
            # 6699780 comes from running wc -l CosmicMutantExport_sorted.tsv, because it takes a minute+ to compute
            # technically, there are 6699779 samples, but we still need to process the last line...
            for row in tqdm.tqdm(cme_reader, desc='building per-gene variant count cache', total=6699780):
                gene_counts[row['Gene name']] += 1
                total_samples += 1

        gene_counts['_total'] = total_samples

        # write out the gene counts to the cache
        with open(GENE_COUNT_CACHE_FILE, 'w') as fp:
            json.dump(gene_counts, fp)

    return gene_counts


if __name__ == '__main__':
    for feature_association in harvest_and_convert(["MDM2"]):
        logging.info(feature_association)
