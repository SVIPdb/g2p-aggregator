#!/usr/bin/python
# -*- coding: utf-8 -*-
import re
import json
from collections import defaultdict
from itertools import chain

from utils_ex.downloading import acquire_files


# load gene names
# ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/json/non_alt_loci_set.json
# other data sources:
# - ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
#   * contains a table mapping NCBI gene IDs to their canonical names + list of synyonyms

# ----------------------------------------------------------------------------------------------------------------
# -- file updating
# ----------------------------------------------------------------------------------------------------------------

data_paths = acquire_files({
    'non_alt_loci_set.json': {
        'path': '/data/non_alt_loci_set.json',
        'url': 'ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/json/non_alt_loci_set.json',
        'compressed': False
    },
    'Homo_sapiens.gene_info': {
        'path': '/data/ncbi_gene/Homo_sapiens.gene_info',
        'url': 'ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz',
        'compressed': True
    }
})


# ----------------------------------------------------------------------------------------------------------------
# -- gene list creation
# ----------------------------------------------------------------------------------------------------------------

GENES = {}
ALIASES = defaultdict(list)

# trim payload, we only need symbol and ensembl
with data_paths['non_alt_loci_set.json'].open() as fp:
    data = json.load(fp)

# also grab synonyms from the clinvar dataset
# we index it by clinvar gene id (aka entrez id) to eliminate collisions
clinvar_genes = {}
with data_paths['Homo_sapiens.gene_info'].open() as fp:
    header = fp.readline()[1:].strip().split('\t')
    for line in fp:
        fields = dict(zip(header, line.strip().split('\t')))
        idx = fields['GeneID']
        if idx not in clinvar_genes:
            clinvar_genes[idx] = fields
        else:
            if str(clinvar_genes[idx]) != str(fields):
                # ensure we're not overwriting anything
                raise Exception("GeneID %s already in clinvar_genes, but other data varies" % fields['GeneID'])


# FIXME: commented out b/c i'm unsure if we'll need to use the uniprot mapping data later on
# currently we were only getting the uniprot id for a gene name, but that's in non_alt_loci_set...

# # read in uniprot ID mapping reference file so we can annotate the gene info with it
# csv.field_size_limit(sys.maxsize)
# uniprot_fp = csv.reader(open('/data/HUMAN_9606_idmapping_selected.tab', 'r'), dialect='excel-tab')
# uniprot_header = uniprot_fp.next()
#
# # maps entrez id to uniprot accession ID, but only for entries w/a PIR column
# # (this is because there are multiple entries w/the same entrez id, most of which are unreviewed;
# # it seems that reviewed entries tend to have a valid PIR column and unreviewed ones don't)
# uniprot_map = dict((r[2], r[0]) for r in uniprot_fp if r[11])


chromosome_extractor = re.compile(r'^([0-9XY]+)')

for doc in data['response']['docs']:
    entrez_id = doc.get('entrez_id', None)
    gene = {
        'symbol': doc['symbol'],
        'ensembl_gene_id': doc.get('ensembl_gene_id', None),
        'entrez_id': entrez_id,
        'location': doc.get('location', None),
        'uniprot_ids': doc.get('uniprot_ids', None),
        'aliases': [],
        'prev_symbols': []
    }
    GENES[doc['symbol']] = [gene]

    if gene['location']:
        result = chromosome_extractor.match(gene['location'])
        if result:
            gene['chromosome'] = result.group(0)

    if gene['ensembl_gene_id']:
        ALIASES[gene['ensembl_gene_id']].append(gene)

    if gene['entrez_id']:
        ALIASES[gene['entrez_id']].append(gene)

    for alias in doc.get('alias_symbol', []):
        ALIASES[alias].append(gene)
        gene['aliases'].append(alias)

    for prev in doc.get('prev_symbol', []):
        ALIASES[prev].append(gene)
        gene['prev_symbols'].append(prev)

    # add in the aliases from clinvar as well, if present
    try:
        synonyms = clinvar_genes[entrez_id]['Synonyms'].split('|')
        for syn in synonyms:
            # ignore placeholder null synonyms
            if syn == '-':
                continue

            ALIASES[syn].append(gene)
            # merge together the existing aliases and any synonyms that aren't listed as previous symbols,
            # and ensure that there are no duplicates via coercing to a set(), then back to a list
            gene['aliases'] = list(set(chain(
                gene['aliases'],
                (syn for syn in synonyms if syn not in gene['prev_symbols'])
            )))
    except KeyError:
        pass

del data


def get_gene(identifier):
    """ return gene for identifier """
    for store in [GENES, ALIASES]:
        genes = store.get(identifier, None)
        if genes and len(genes) == 1:
            return genes

    raise ValueError('gene reference does not exist or refers to multiple genes')


def normalize_feature_association(feature_association):
    """ add gene_identifiers array to feature_association """
    feature_association['gene_identifiers'] = [
        get_gene(gene_symbol)[0] for gene_symbol in feature_association['genes']
    ]
