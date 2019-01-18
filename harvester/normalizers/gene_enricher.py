#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import json
import csv
from collections import defaultdict

# load gene names
# ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/json/non_alt_loci_set.json
GENES = {}
ALIASES = defaultdict(list)

# trim payload, we only need symbol and ensembl
data = json.load(open('../data/non_alt_loci_set.json'))

# FIXME: commented out b/c i'm unsure if we'll need to use the uniprot mapping data later on
# currently we were only getting the uniprot id for a gene name, but that's in non_alt_loci_set...

# # read in uniprot ID mapping reference file so we can annotate the gene info with it
# csv.field_size_limit(sys.maxsize)
# uniprot_fp = csv.reader(open('../data/HUMAN_9606_idmapping_selected.tab', 'r'), dialect='excel-tab')
# uniprot_header = uniprot_fp.next()
#
# # maps entrez id to uniprot accession ID, but only for entries w/a PIR column
# # (this is because there are multiple entries w/the same entrez id, most of which are unreviewed;
# # it seems that reviewed entries tend to have a valid PIR column and unreviewed ones don't)
# uniprot_map = dict((r[2], r[0]) for r in uniprot_fp if r[11])

for doc in data['response']['docs']:
    entrez_id = doc.get('entrez_id', None)
    gene = {
        'symbol': doc['symbol'],
        'ensembl_gene_id': doc.get('ensembl_gene_id', None),
        'entrez_id': entrez_id,
        'location': doc.get('location', None),
        'uniprot_ids': doc.get('uniprot_ids', None)
    }
    GENES[doc['symbol']] = [gene]

    if gene['ensembl_gene_id']:
        ALIASES[gene['ensembl_gene_id']].append(gene)

    if gene['entrez_id']:
        ALIASES[gene['entrez_id']].append(gene)

    for alias in doc.get('alias_symbol', []):
        ALIASES[alias].append(gene)

    for prev in doc.get('prev_symbol', []):
        ALIASES[prev].append(gene)

del data


def get_gene(identifier):
    """ return gene for identifier """
    for store in [GENES, ALIASES]:
        genes = store.get(identifier, None)
        if genes and len(genes) == 1:
            return genes
        else:
            raise ValueError('gene reference does not exist or refers to multiple genes')


def normalize_feature_association(feature_association):
    """ add gene_identifiers array to feature_association """
    gene_identifiers = []
    for gene_symbol in feature_association['genes']:
        try:
            gene = get_gene(gene_symbol)
        except:
            gene = None
        if (gene):
            gene_identifiers.extend(gene)
    feature_association['gene_identifiers'] = gene_identifiers
