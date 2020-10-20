#!/usr/bin/python
import csv
import itertools
import json
import operator
import os
import re

import logging
import sqlite3
from collections import defaultdict

import tqdm

# hgvs stuff
import hgvs.normalizer
import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.assemblymapper


# globals that we use for parsing variant representations
from lookups.accession_mapping import ensembl_txac_to_refseq
from utils_ex.formatting import capitalize_words
from utils_ex.tqdm_logger import std_out_err_redirect_tqdm

coord_matcher = re.compile(r'(?P<chrom>[^:]+):(?P<start>[0-9]+)-(?P<end>[0-9]+)')

# these shared assembly mappers will allow us to convert HGVS g. variants to c. and p. later on
hdp = hgvs.dataproviders.uta.connect()
hgnorm = hgvs.normalizer.Normalizer(hdp)
hgvsparser = hgvs.parser.Parser()
am = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh37', normalize=True)

COSMIC_MUTANTS_DB_FILE = '/data/CosmicMutantExport.sqlite'


def parse_cosmic_cds_change(in_cds):
    """
    Parses the COSMIC Mutation CDS string to extract the 'ref' and 'alt' sections. Returns a dict with 'ref' and 'alt'
    keys.

    We can't use a regular hgvs-parsing library since this is non-standard nomenclature.

    :param in_cds: the coding-sequence change to parse.
    :return: A dict with 'ref' and 'alt' keys.

    >>> parse_cosmic_cds_change("c.570C>T")
    {'alt': 'T', 'ref': 'C'}
    >>> parse_cosmic_cds_change("c.570CC>TT")
    {'alt': 'TT', 'ref': 'CC'}
    >>> parse_cosmic_cds_change("c.393delC")
    {'alt': None, 'ref': 'C'}
    >>> parse_cosmic_cds_change("c.64_66delTTG")
    {'alt': None, 'ref': 'TTG'}
    """
    patterns = [
        r"(?P<ref>[ACTG]+)\>(?P<alt>[ACTG]+)$",
        r"del(?P<ref>[ACTG]+)$",
        r"ins(?P<alt>[ACTG]+)$",
    ]

    for pattern in patterns:
        m = re.search(pattern, in_cds)
        if m:
            groups = m.groupdict()
            return {
                'ref': groups.get('ref'),
                'alt': groups.get('alt')
            }

    return None


# ==========================================================================================
# === harvester implementation
# ==========================================================================================

def _get_rows(cursor):
    while True:
        results = cursor.fetchmany()
        if not results:
            break
        for result in results:
            yield result

def harvest(genes=None, conn=None):
    # if the connection hasn't been opened, open it; otherwise, reuse the existing one
    if not conn:
        conn = sqlite3.connect(COSMIC_MUTANTS_DB_FILE)

    cursor = conn.cursor()

    if genes:
        cursor.execute('select * from variants where "Gene name" in (%s) order by "Gene name"' % ','.join('?' * len(genes)), genes)
    else:
        cursor.execute('select * from variants order by "Gene name"')

    cur_gene = None
    refseq_ac = None

    for sample in _get_rows(cursor):
        if not cur_gene or sample['Gene name'] != cur_gene:
            refseq_ac = ensembl_txac_to_refseq(sample['Accession Number'])

        yield sample, refseq_ac


# noinspection PyTypeChecker
def convert(sample, refseq_ac, tq):
    # NOTE: we used to grab the first variant in the Mutation AA group, but later noticed that not every
    # sample under a 'Mutation AA' entry had the same fields for, say, "Mutation CDS" (which in retrospect
    # makes perfect sense.) anyway, we went back to processing each sample individually instead of trying to
    # group them under the same AA change, which is incidentally faster due to not having to sort the samples anymore.

    if (
        not sample['Mutation AA'] or not sample['Mutation CDS']
        or '?' in sample['Mutation AA'] or '?' in sample['Mutation CDS']
        # or '_' in sample['Mutation CDS']
    ):
        # afaik, we can't parse variants that have a question mark in their cDNA representation
        # logging.warn("Skipping COSMIC variant %s (AA: %s) b/c of missing cDNA rep: %s" % (
        #     sample['Mutation ID'], sample['Mutation AA'], mapped_cds))
        tq.update()  # move the counter ahead anyway
        return None

    # the following is shared by all samples under this variant, although we emit an evidence item
    # for each sample
    match = coord_matcher.match(sample['Mutation genome position'])
    coords = match.groupdict() if match else {}

    # cosmic uses a weird format for non-single-nucleotide variants, e.g. it uses 98_99CC>TT when it should be
    #  98_99delinsTT. similarly, it shows the deleted bases in, e.g., c.253_254delGG; the spec states that the
    #  deleted bases shouldn't be included. since we just need to extract ref and alt here, we can do some ad-hoc
    #  parsing for each case to get what we need

    # var = hgvsparser.parse_c_variant("%s:%s" % (sample['Accession Number'], sample['Mutation CDS']))
    # pos_info = {
    #     'ref': var.posedit.edit.ref,
    #     'alt': var.posedit.edit.alt
    # }

    pos_info = parse_cosmic_cds_change(sample['Mutation CDS'])

    # TODO: we should use the accession field plus the CDS change to find the genomic position,
    #  which is clumsily accomplished later via the location normalizer (just with the protein change, which
    #  frankly is probably wrong...)

    feature = {
        'geneSymbol': sample['Gene name'],
        'entrez_id': sample['HGNC ID'],
        'start': coords.get('start'),
        'end': coords.get('stop'),
        'referenceName': sample['GRCh'],
        'refseq': refseq_ac,
        'isoform': sample['Accession Number'],
        'chromosome': coords.get('chrom'),
        'ref': pos_info['ref'] if pos_info else None,
        'alt': pos_info['alt'] if pos_info else None,
        'name': sample['Mutation AA'][2:],
        'description': "%s %s" % (sample['Gene name'], sample['Mutation AA'][2:]),
        'biomarker_type': sample['Mutation Description'],
        'strand': sample['Mutation strand'],

        'somatic_status': sample['Mutation somatic status']
    }

    association = {
        'variant_name': feature['name'],
        'description': 'n/a',  # FIXME: what's a good description for this?

        'extras': {
            'fathmm_prediction': sample['FATHMM prediction'],
            'fathmm_score': sample['FATHMM score'],
        },

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
    source_url = 'https://cancer.sanger.ac.uk/cosmic/mutation/overview?id=%s' % sample['MUTATION_ID'],

    feat_assoc = {
        'genes': [sample['Gene name']],  # there's only ever one gene/variant
        'features': [feature],
        'feature_names': feature["geneSymbol"] + ' ' + feature["name"],
        'association': association,
        'source': 'cosmic',
        'source_url': source_url,
        'cosmic': sample
    }

    tq.update()

    return feat_assoc

def harvest_and_convert(genes):
    conn = sqlite3.connect(COSMIC_MUTANTS_DB_FILE)
    curs = conn.cursor()

    if genes:
        curs.execute('select count(*) from variants where "Gene name" in (%s)' % ','.join('?' * len(genes)), genes)
    else:
        curs.execute("select count(*) from variants")

    sample_count = curs.fetchone()[0]

    with std_out_err_redirect_tqdm() as orig_stdout:
        with tqdm.tqdm(total=sample_count, desc="harvesting %s" % (genes,), file=orig_stdout, dynamic_ncols=True) as tq:
            for sample, refseq_ac in harvest(genes):
                feat_assoc = convert(sample, refseq_ac, tq=tq)
                if feat_assoc:
                    yield feat_assoc


if __name__ == '__main__':
    for feature_association in harvest_and_convert(["MDM2"]):
        logging.info(feature_association)
