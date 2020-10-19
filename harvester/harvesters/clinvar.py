#!/usr/bin/python
import functools
import json
import os
from itertools import ifilter
from operator import itemgetter

import requests
import copy
from lookups import evidence_label as el, evidence_direction as ed
import logging
from tqdm import tqdm

from lxml import etree


# --------------------------------------------------------------------------------------------------------------
# --- support methods
# --------------------------------------------------------------------------------------------------------------
from utils_ex.instrumentation import OpTimeLogger

three_to_one = {
    '*': 'X',
    'Ter': 'X',  # this appears to be an NCBI-specific thing
    'Ala': 'A',
    'Arg': 'R',
    'Asn': 'N',
    'Asp': 'D',
    'Asx': 'B',
    'Cys': 'C',
    'Gln': 'Q',
    'Glu': 'E',
    'Glx': 'Z',
    'Gly': 'G',
    'His': 'H',
    'Ile': 'I',
    'Leu': 'L',
    'Lys': 'K',
    'Met': 'M',
    'Phe': 'F',
    'Pro': 'P',
    'Sec': 'U',
    'Ser': 'S',
    'Thr': 'T',
    'Trp': 'W',
    'Tyr': 'Y',
    'Val': 'V',
    'Xaa': 'Xa'
}


def remap_prots(three_prot_change):
    """
    >>> remap_prots("Val600Glu")
    'V600E'
    """
    return reduce(lambda acc, x: acc.replace(x[0], x[1]), three_to_one.items(), three_prot_change)


def first(l, default=None):
    return next(iter(l), default)


DISEASE_DB_URLS = {
    'OMIM': 'https://www.omim.org/entry/%s',
    'Genetic Alliance': 'https://www.diseaseinfosearch.org/%s',
    'OrphaNet': 'https://www.orpha.net/consor/cgi-bin/OC_Exp.php?Expert=%s',
    'SNOMED CT': 'http://purl.bioontology.org/ontology/SNOMEDCT/%s',
    'Human Phenotype Ontology': 'https://hpo.jax.org/app/browse/term/HP:%s'
}


def fast_iter(context):
    """
    http://lxml.de/parsing.html#modifying-the-tree
    Based on Liza Daly's fast_iter
    http://www.ibm.com/developerworks/xml/library/x-hiperfparse/
    See also http://effbot.org/zone/element-iterparse.htm
    """
    for event, elem in context:
        yield event, elem
        # It's safe to call clear() here because no descendants will be
        # accessed
        elem.clear()
        # Also eliminate now-empty references from the root node to elem
        for ancestor in elem.xpath('ancestor-or-self::*'):
            while ancestor.getprevious() is not None:
                del ancestor.getparent()[0]
    del context


def is_valid_record(elem, gene_set):
    if gene_set:
        these_genes = set(elem.xpath("InterpretedRecord/*/GeneList/Gene/@Symbol", smart_strings=False))
        if gene_set.isdisjoint(these_genes):
            return False

    # if it's a non-human sample or not current, we should probably bail, too
    if elem.find('RecordStatus').text != 'current' or elem.find('Species').text != 'Homo sapiens':
        return False

    # finally, if it doesn't have a good protein change, also bail
    if not len(elem.xpath('InterpretedRecord/*/HGVSlist/HGVS[@Type="coding"]/ProteinExpression')):
        return False

    return True


# --------------------------------------------------------------------------------------------------------------
# --- harvester implementation
# --------------------------------------------------------------------------------------------------------------

# if true, uses the cached variant counts if available
USE_COUNTS_CACHE = True

# the XML file to use for clinvar data, retrieved by ../data/clinvar/get_clinvar_latest.sh
CLINVAR_XML_FILE = "/data/clinvar/ClinVarVariationRelease_00-latest.xml"


def harvest(genes):
    """ given an array of gene symbols, harvest them from clinvar"""
    gene_set = set(genes) if genes else None

    is_valid_record_w_genes = functools.partial(is_valid_record, gene_set=gene_set)

    with open(CLINVAR_XML_FILE) as fp:
        cache_file_path = os.path.join(os.path.dirname(CLINVAR_XML_FILE), ".var_counts_cache.json")

        # verify that the cache is valid, i.e. it refers to the correct file, gene set, etc.
        cache_valid = False
        if  os.path.exists(cache_file_path):
            with open(cache_file_path, "r") as cache_fp:
                cached_figures = json.load(cache_fp)

                cache_valid = (
                    all(k in cached_figures for k in ('filename', 'gene_set', 'total_variants', 'matching_variants',)) and
                    cached_figures['filename'] == CLINVAR_XML_FILE and
                    set(cached_figures['gene_set']) == set(gene_set)
                )

        if cache_valid and USE_COUNTS_CACHE:
            total_variants = cached_figures['total_variants']
            matching_variants = cached_figures['matching_variants']
        else:
            # first, count off the elements so we can make a progress bar
            total_variants = 0
            matching_variants = 0

            with OpTimeLogger(name="Clinvar Variant Counting"):
                print("Calculating elements in ClinVar input set...")
                ctx = etree.iterparse(fp, events=("end",), tag="VariationArchive")
                for action, elem in fast_iter(ctx):
                    # first, if we have a gene list, check if this one's in it and skip it if it's not
                    total_variants += 1
                    if is_valid_record_w_genes(elem):
                        matching_variants += 1

            # write the results to the cache, too
            with open(cache_file_path, "w") as cache_fp:
                json.dump({
                    "filename": CLINVAR_XML_FILE,
                    "gene_set": gene_set,
                    "total_variants": total_variants,
                    "matching_variants": matching_variants
                }, cache_fp)

            # and reset the file pointer so we can re-read it all
            fp.seek(0)

        print("Matching vs. total variants: %d/%d" % (matching_variants, total_variants))

        ctx = etree.iterparse(fp, events=("end",), tag="VariationArchive")

        for action, elem in tqdm(fast_iter(ctx), desc="ClinVar", total=total_variants):
            if not is_valid_record_w_genes(elem):
                continue

            # basically, just extract each VariationArchive, which represents a single variant.

            # (we can't aggregate these under a single gene like we do with other harvesters b/c we're not
            # guaranteed to get variants in gene order and sorting/aggregating them in a 7GB XML file
            # would be expensive anyway, and obviate the advantages of SAX-parsing the monster...)
            yield elem


def extract_phenotype(y):
    preferred_term = y.xpath('Name/ElementValue[@Type="Preferred"]')[0]
    xref = preferred_term.getparent().find('XRef')

    return {
        # ontology link, e.g. http://purl.obolibrary.org/obo/DOID_4329
        'source': DISEASE_DB_URLS[xref.attrib['DB']] % xref.attrib['ID'] if xref is not None and xref.attrib['DB'] in DISEASE_DB_URLS else None,
        # disease name, e.g. 'Erdheim-Chester disease'
        'term': preferred_term.text,
        # some kind of ontology reference, e.g. 'DOID:4329'
        'id': "%s:%s" % (xref.attrib['DB'], xref.attrib['ID']) if xref is not None else None,
        # e.g., 'histiocytosis'
        'family': None,
        # a readable version of the term, e.g. 'Erdheim-Chester disease'
        'description': preferred_term.text
    }


def convert(root):
    """
    Transforms a single VariationArchive entry from ClinVarVariation into a feat_assoc object that we can import into our silo.
    :param root: the VariationArchive lxml element
    :return: a feat_assoc object
    """

    # general values
    source_url = 'https://www.ncbi.nlm.nih.gov/clinvar/variation/%d' % int(root.attrib['VariationID'])

    # -- features-related values --
    genes = root.xpath('InterpretedRecord/*/GeneList/Gene')
    gene_symbols = [x.attrib['Symbol'] for x in genes]

    try:
        grch37_pos = root.xpath('//SimpleAllele/Location/SequenceLocation[@Assembly="GRCh37"]')[0]
    except IndexError:
        # this variant has no GRCh37 values; we can't include it, so we have to skip it
        logging.warn("ClinVar entry has no GRCh37 SequenceLocation, skipping...")
        return

    # FIXME: when running over the full clinvar set, this expression throws an error...maybe because there's no hg19 expression?
    #  ideally we'd have better error reporting for all these deep, fragile references into the structure.
    # FIXME: we should also decide what the appropriate behavior is if something critical is missing: do we ignore the
    #  the variant, or attempt to infer these values from other fields? i'm going with ignore for hgvs_g, but maybe not for _c and _p
    #  since it could be in a non-coding or inter-gene region, or it could be a variant that involves multiple genes,
    #  i.e a fusion. hrm...
    try:
        hgvs_g = root.xpath('InterpretedRecord/*/HGVSlist/HGVS[@Assembly="GRCh37" and @Type="genomic, top-level"]/NucleotideExpression/Expression/text()')[0]
    except IndexError:
        logging.warn("Unable to locate HGVS.g field in variant@ %s (gene %s), skipping..." % (source_url, ", ".join(gene_symbols)))
        return

    # find the hgvs.c representation that uses the NCBI transcript accessions (prefixed by NM_), not LRG
    # (since pyhgvs can't handle projections w/it afaik)
    hgvs_c_first_node = first(root.xpath('InterpretedRecord/*/HGVSlist/HGVS[@Type="coding"]/NucleotideExpression[starts-with(@sequenceAccessionVersion, "NM_")]'))
    # looking for <HGVS type="protein"> gives UniProt accession numberrs for some reason
    hgvs_p_first_node = first(root.xpath('InterpretedRecord/*/HGVSlist/HGVS[@Type="coding"]/ProteinExpression'))

    prot_change = remap_prots(hgvs_p_first_node.attrib['change'][2:]) if hgvs_p_first_node is not None else None

    # FIXME: consider skipping the entry if we can't get the protein change,
    #  since it's likely intronic or outside the gene

    # FIXME: currently some variants (looking at you, EGFR) have multiple associated genes
    #  for now i'm going to ignore them, but we'll revisit them later
    if len(gene_symbols) != 1:
        logging.warn("Non-singular gene list specified (%s), skipping..." % ", ".join(gene_symbols))
        return

    # each variant will produce multiple assocations for the same set of features
    basis = {
        'source': 'clinvar',
        'source_url': source_url,
        # should be the URL of the variant

        'genes': gene_symbols,  # first element required

        'features': [
            {
                'geneSymbol': gene_symbols[0],
                'name': prot_change if prot_change else hgvs_g.split(':')[1],  # typically <PROTEIN-CHANGE>, e.g. "V600E"
                'sequence_ontology': {  # optional, injected by biomarker_normalizer, relies on biomarker_type
                    'soid': 'required',
                    'so_name': 'required',
                    'hierarchy': None
                },

                'description': "%s %s" % (gene_symbols[0], prot_change),
                # typically <GENE-SYMBOL> <NAME>, e.g. "BRAF V600E"
                'referenceName': 'GRCh37',  # the assembly, in all cases so far GRCh37
                'refseq': hgvs_c_first_node.attrib['sequenceAccessionVersion'] if hgvs_c_first_node is not None else None,
                'isoform': None,  # presumbly an ensembl accession, starts with ENST
                'biomarker_type': root.xpath('InterpretedRecord/SimpleAllele/VariantType/text()')[0],
                'chromosome': grch37_pos.attrib['Chr'],
                'start': grch37_pos.attrib['start'],
                'end': grch37_pos.attrib['stop'],
                'ref': grch37_pos.attrib.get('referenceAllele', grch37_pos.attrib.get('referenceAlleleVCF')),
                'alt': grch37_pos.attrib.get('alternateAllele', grch37_pos.attrib.get('alternateAlleleVCF')),
                'hgvs_g': hgvs_g,
                'hgvs_c': hgvs_c_first_node.find('Expression').text if hgvs_c_first_node is not None and len(hgvs_c_first_node) else None,
                'hgvs_p': hgvs_p_first_node.find('Expression').text if hgvs_p_first_node is not None and len(hgvs_p_first_node) else None,
                'dbsnp_ids': root.xpath('//XRef[@Type="rs" and @DB="dbSNP"]/@ID'),
                'myvariant_hg19': None,  # injected by normalizer
                'mv_info': None,  # injected by normalizer
                'crawl_status': None  # used to report crawling issues, injected as defaultdict by utils_ex.instrumentation
            }
        ]
    }

    for rec in extract_rcvs(basis, gene_symbols, prot_change, root):
        yield rec


def extract_rcvs(basis, gene_symbols, prot_change, root):
    for rcv_accession in root.xpath('InterpretedRecord/RCVList/RCVAccession'):
        feat_assoc = basis.copy()

        condition = first(rcv_accession.xpath("InterpretedConditionList/InterpretedCondition"))

        # noinspection PyTypeChecker
        feat_assoc['association'] = {
            'description': "%s %s" % (gene_symbols[0], prot_change),
            'drug_labels': None,  # a comma-delimited list of drug names (none for clinvar)
            'drug_interaction_type': None,  # only supplied by CIViC afaict, and only ever 'Substitutes' or null
            'variant_name': prot_change,
            # 'source_link': source_url,  # should be the URL of the evidence item(s) # FIXME: validate this
            # 'source_link': 'https://www.ncbi.nlm.nih.gov/clinvar/?term="%s"' % x['ID'],
            'source_link': 'https://www.ncbi.nlm.nih.gov/clinvar/%(Accession)s.%(Version)s/' % rcv_accession.attrib,
            'evidence_type': 'Predisposing',  # all clinvar records are about clinical significance
            'evidence_direction': 'Supports',
            'clinical_significance': rcv_accession.attrib['Interpretation'],
            'evidence_level': rcv_accession.attrib['ReviewStatus'],
            # FIXME: what's the evidence level for clinvar? i can't find the stars in the XML
            'crawl_status': {},
            'extras': {
                'num_submissions': rcv_accession.attrib['SubmissionCount']
            },

            'phenotypes': [
                {
                    # ontology link, e.g. http://purl.obolibrary.org/obo/DOID_4329
                    'source': (
                        DISEASE_DB_URLS[condition.attrib['DB']] % condition.attrib['ID']
                        if condition is not None and 'DB' in condition.attrib and condition.attrib['DB'] in DISEASE_DB_URLS
                        else None
                    ),
                    # disease name, e.g. 'Erdheim-Chester disease'
                    'term': condition.text,
                    # some kind of ontology reference, e.g. 'DOID:4329'
                    'id': (
                        "%s:%s" % (condition.attrib['DB'], condition.attrib['ID'])
                        if condition is not None and 'DB' in condition.attrib
                        else None
                    ),
                    # e.g., 'histiocytosis'
                    'family': None,
                    # a readable version of the term, e.g. 'Erdheim-Chester disease'
                    'description': condition.text
                }
            ],

            'evidence': [],

            'environmentalContexts': []
        }

        yield feat_assoc


def extract_scvs(basis, gene_symbols, prot_change, root):
    for scv in root.xpath('//ClinicalAssertionList/ClinicalAssertion'):
        feat_assoc = basis.copy()

        acc = first(scv.xpath('ClinVarAccession'))
        condition = first(scv.xpath("TraitSet/Trait[@Type='Disease']/Name/ElementValue[@Type='Preferred']"))

        # noinspection PyTypeChecker
        feat_assoc['association'] = {
            'description': "%s %s" % (gene_symbols[0], prot_change),
            'drug_labels': None,  # a comma-delimited list of drug names (none for clinvar)
            'drug_interaction_type': None,  # only supplied by CIViC afaict, and only ever 'Substitutes' or null
            'variant_name': prot_change,
            # 'source_link': source_url,  # should be the URL of the evidence item(s) # FIXME: validate this
            'source_link': 'https://www.ncbi.nlm.nih.gov/clinvar/?term="%(Accession)s"' % acc.attrib,
            # 'source_link': 'https://www.ncbi.nlm.nih.gov/clinvar/%(Accession)s.%(Version)s/' % scv.attrib,
            'evidence_type': 'Predisposing',  # all clinvar records are about clinical significance
            'evidence_direction': 'Supports',
            'clinical_significance': first(scv.xpath('Interpretation/Description/text()')),
            'evidence_level': first(scv.xpath('ReviewStatus/text()')),
            # FIXME: what's the evidence level for clinvar? i can't find the stars in the XML
            'crawl_status': {},
            'extras': None,

            'phenotypes': [
                {
                    # ontology link, e.g. http://purl.obolibrary.org/obo/DOID_4329
                    'source': (
                        DISEASE_DB_URLS[condition.attrib['DB']] % condition.attrib['ID']
                        if condition is not None and condition.attrib['DB'] in DISEASE_DB_URLS
                        else None
                    ),
                    # disease name, e.g. 'Erdheim-Chester disease'
                    'term': condition.text,
                    # some kind of ontology reference, e.g. 'DOID:4329'
                    'id': "%s:%s" % (condition.attrib['DB'], condition.attrib['ID']) if condition is not None else None,
                    # e.g., 'histiocytosis'
                    'family': None,
                    # a readable version of the term, e.g. 'Erdheim-Chester disease'
                    'description': condition.text
                }
            ],

            'evidence': [
                {
                    'info': {  # optional
                        'publications': [
                            'http://www.ncbi.nlm.nih.gov/pubmed/%s' % y[0]
                            for y in scv.xpath('Citation/ID[@Source="PubMed"]/text()')
                        ]
                    },
                    'evidenceType': {  # optional
                        'id': "%s-%s" % (
                            gene_symbols[0],
                            scv.xpath('ConditionList/TraitSet[@Type="Disease"]/Trait[@Type="Disease"]/Name/ElementValue[@Type="Preferred"]/text()')[0]
                        ),  # e.g., 'BRAF-Erdheim-Chester Disease' (gene name, disease name)
                        'sourceName': 'clinvar'
                    }
                }
            ],

            'environmentalContexts': []
        }

        yield feat_assoc


def harvest_and_convert(genes):
    """ get data from clinvar, convert it to ga4gh and return via yield """
    for entry in harvest(genes):
        # print "harvester_yield {}".format(gene_data.keys())
        for feat_assoc in convert(entry):
            # print "Yielded; %s" % feat_assoc
            yield feat_assoc


if __name__ == '__main__':
    for feature_association in harvest_and_convert(["MDM2"]):
        logging.info(feature_association)
