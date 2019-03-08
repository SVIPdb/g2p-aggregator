#!/usr/bin/python
import csv
import itertools
import json
import operator
import re

import logging
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


# ==========================================================================================
# === harvester implementation
# ==========================================================================================

def harvest(genes=None):
    with open("./CosmicMutantExport_sorted.tsv") as cme_fp:
        cme_reader = csv.DictReader(cme_fp, dialect='excel-tab')

        it = itertools.groupby(cme_reader, operator.itemgetter('Gene name'))
        if genes is not None:
            it = itertools.ifilter(lambda x: x[0] in genes, it)

        for gene_symbol, variants in it:
            variant = variants.next()  # peek at the first variant in the list...
            full_variants = itertools.chain([variant], variants)  # ...but put it back when we iterate later

            # all genes have the same accession number, so we can do it on the gene level here
            refseq_ac = ensembl_txac_to_refseq(variant['Accession Number'])

            yield {
                'gene_symbol': gene_symbol,
                'refseq_ac': refseq_ac,
                'variants': full_variants
            }


# noinspection PyTypeChecker
def convert(gene_data, tq):
    # try:

    for variant_aa, samples in itertools.groupby(gene_data['variants'], operator.itemgetter('Mutation AA')):
        # we know all samples for the same variant share these fields, so we'll just grab the first one
        # (we later have to put it back on the generator, which is accomplished via chain() later on)
        variant = samples.next()
        mapped_cds = variant['Mutation CDS']

        if not variant:
            raise Exception("No variants found for AA %s!" % variant_aa)

        if '?' in mapped_cds:
            # afaik, we can't parse variants that have a question mark in their cDNA representation
            logging.warn("Skipping COSMIC variant %s (AA: %s) b/c of missing cDNA rep: %s" % (
                variant['Mutation ID'], variant_aa, mapped_cds))
            tq.update(len(list(samples)) + 1)
            continue

        # the following is shared by all samples under this variant, although we emit an evidence item
        # for each sample
        match = coord_matcher.match(variant['Mutation genome position'])
        coords = match.groupdict() if match else {}

        # # for some reason the incorrect notation for a delins that specifies the source rep is common in COSMIC
        # # (e.g., c.1798_1799GT>AA, which should instead be c.1798_1799GTdelinsAA)
        # # FIXME: should we be attempting to map these to their appropriate hgvs rep, as below, or should we ignore them?
        # bad_delins_match = re.match(r'c.(?P<range>[0-9]+_[0-9]+)[A-Z]+>(?P<alt>[A-Z]+)', mapped_cds)
        # if bad_delins_match:
        #     mapped_cds = 'c.%(range)sdelins%(alt)s' % bad_delins_match.groupdict()
        # in order to produce ref, alt, we'll need to parse the HGVS cDNA field
        # seqvar = None
        # try:
        #     seqvar = hgvsparser.parse_c_variant("%s:%s" % (gene_data['refseq_ac'], mapped_cds))
        # except:
        #     logging.exception("Couldn't parse hgvs cDNA string, skipping")
        #     tq.update(len(list(samples)) + 1)
        #     continue

        # this is way faster and gets us basically the same thing as using the hgvs library, since we don't even have
        # the right accession anyway...
        pos_info = parse_hgvc_c(variant['Mutation CDS'])

        feature = {
            'geneSymbol': gene_data['gene_symbol'],
            'entrez_id': variant['HGNC ID'],
            'start': coords.get('start'),
            'end': coords.get('stop'),
            'referenceName': variant['GRCh'],
            'refseq': gene_data['refseq_ac'],
            'isoform': variant['Accession Number'],
            'chromosome': coords.get('chrom'),
            'ref': pos_info['ref'] if pos_info else None,
            'alt': pos_info['alt'] if pos_info else None,
            'name': variant['Mutation AA'][2:],
            'description': "%s %s" % (gene_data['gene_symbol'], variant['Mutation AA'][2:]),
            'biomarker_type': variant['Mutation Description']
        }

        for sample in itertools.chain([variant], samples):
            # TODO: figure out the evidence label and level for all cosmic variants

            association = {
                'variant_name': feature['name'],
                'description': 'n/a',  # FIXME: what's a good description for this?

                'extras': json.dumps({
                    'fathmm_prediction': variant['FATHMM prediction'],
                    'fathmm_score': variant['FATHMM score'],
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
            source_url = 'https://cancer.sanger.ac.uk/cosmic/mutation/overview?id=%s' % variant['Mutation ID'][4:],

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
    # 6699780 comes from running wc -l CosmicMutantExport_sorted.tsv, because it takes a minute+ to compute otherwise
    if genes:
        print("Counting samples to insert...")
        sample_count = 0
        with open("./CosmicMutantExport_sorted.tsv") as cme_fp:
            cme_reader = csv.DictReader(cme_fp, dialect='excel-tab')
            for row in tqdm.tqdm(cme_reader, total=6699780):
                if not genes or row['Gene name'] in genes:
                    sample_count += 1
        print("done!")
    else:
        sample_count = 6699780

    with std_out_err_redirect_tqdm() as orig_stdout:
        with tqdm.tqdm(total=sample_count, desc="harvesting %s" % genes, file=orig_stdout, dynamic_ncols=True) as tq:
            for gene_data in harvest(genes):
                for feat_assoc in convert(gene_data, tq):
                    yield feat_assoc


# main
if __name__ == '__main__':
    for feature_association in harvest_and_convert(["MDM2"]):
        logging.info(feature_association)
