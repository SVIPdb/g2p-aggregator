import sys
import re
import pandas
import json
import copy
import logging

import cosmic_lookup_table
import evidence_label as el
import evidence_direction as ed
import mutation_type as mut

LOOKUP_TABLE = None


""" https://www.cancergenomeinterpreter.org/biomarkers """


def _get_evidence(gene_ids, path='../data'):
    """ load tsv, yield object where values_cols are lists """
    df = pandas.read_table(path + '/cgi_biomarkers_per_variant.tsv')
    # change nan to blank string
    df = df.fillna('')
    # if no gene list return all
    if gene_ids:
        df = df.loc[df['Gene'].isin(gene_ids)]

    groupby_cols = [
        'Alteration',
        'Alteration type',
        'Assay type',
        'Association',
        'Biomarker',
        'Curator',
        'Drug',
        'Drug family',
        'Drug full name',
        'Drug status',
        'Evidence level',
        'Gene',
        'Metastatic Tumor Type',
        'Primary Tumor type',
        'Source',
        'Targeting'
        ]
    values_cols = [
        'gDNA',
        'cDNA',
        'individual_mutation',
        'info',
        'region',
        'strand',
        'transcript'
        ]
    grouped = df.groupby(groupby_cols)

    for name, group in grouped:
        evidence = {}
        for p in groupby_cols:
            evidence[p] = group[p].values[0]
        for p in values_cols:
            evidence[p] = list(group[p].values)
        yield evidence


def convert(evidence):
    """
    ['Primary Tumor type', 'Drug family', 'Alteration type', 'Targeting',
    'Assay type', 'Evidence level', 'Biomarker', 'Drug', 'Alteration',
    'Source', 'Curator', 'Comments', 'Drug status', 'Drug full name',
    'TCGI included', 'Curation date', 'Gene', 'Metastatic Tumor Type',
    'Association']
    {'Primary Tumor type': 'GIST', 'Drug family': '[HSP90 inhibitor]',
     'Alteration type': 'MUT', 'Targeting': nan, 'Assay type': nan,
     'Evidence level': 'Pre-clinical',
     'Biomarker': 'KIT mutation in exon 9 or 17',
     'Drug': '[]',
     'Alteration': 'KIT:788-828,449-514',
     'Source': 'PMID:21737509', 'Curator': 'RDientsmann',
     'Comments': nan, 'Drug status': nan,
     'Drug full name': 'HSP90 inhibitors',
     'TCGI included': True, 'Curation date': '01/16',
     'Gene': 'KIT', 'Metastatic Tumor Type': nan,
     'Association': 'Responsive'}
    """

    def split_gDNA(gDNA):
        ''' Split gDNA field of the form 'chr9:g.133747588G>C' and return dictionary. '''

        # TODO: handle non-SNPs like chr1:g.43815009_43815010delGGinsTT
        try:
            chrom, remainder = gDNA.split(':g.')
            if chrom.startswith('chr'):
                chrom = chrom[3:]
            start = re.search(r'(\d+)', remainder).group()
            ref, alt = remainder[len(start):].split(">")
            return {
                'chromosome': str(chrom),
                'start': int(start),
                'ref': ref,
                'alt': alt
            }
        except Exception as e:
            return {}

    # need gobal LOOKUP_TABLE variable
    global LOOKUP_TABLE

    # Create document for insertion.
    genes = re.split('\W+', evidence['Gene'])

    alteration_types = re.split('\W+', evidence['Alteration type'])
    genes = filter(len, genes)
    alteration_types = filter(len, alteration_types)
    features = []

    gDNA = evidence['gDNA']
    indiv_mut = evidence['individual_mutation']
    if len(gDNA) != 0 or len(indiv_mut) != 0:
        if not LOOKUP_TABLE:
            LOOKUP_TABLE = cosmic_lookup_table.CosmicLookup(
                "./cosmic_lookup_table.tsv")
        for idx, emut in enumerate(indiv_mut):
            # get genomic locus from COSMIC; if mutation not in COSMIC,
            # get locus info from given gDNA instead
            feature = {}
            if len(emut) > 0:
                gene, prot = emut.split(':')
                matches = LOOKUP_TABLE.get_entries(gene, prot)
                if len(matches) > 0:
                    # FIXME: just using the first match for now;
                    # it's not clear what to do if there are multiple matches.
                    match = matches[0]
                    feature['chromosome'] = str(match['chrom'])
                    feature['start'] = match['start']
                    feature['end'] = match['end']
                    feature['ref'] = match['ref']
                    feature['alt'] = match['alt']
                    feature['referenceName'] = str(match['build'])
            else:
                if len(gDNA[idx]) == 0:
                    continue
                feature = split_gDNA(gDNA[idx])
            alteration_type = alteration_types[0]
            if len(alteration_types) > idx:
                alteration_type = alteration_types[idx]
            if not gene:
                gene = genes[0]
                if len(genes) > idx:
                    gene = genes[idx]
            feature['biomarker_type'] = mut.norm_biomarker(alteration_type, emut)
            feature['name'] = emut
            feature['description'] = emut.replace(':', '  ')
            feature['referenceName'] = 'GRCh37'
            feature['geneSymbol'] = gene
            features.append(feature)

    if len(features) == 0:
        description_parts = re.split(' +|:|__', evidence['Biomarker'].strip())
        features.append({
            'description': ' '.join(description_parts),
            'name': ' '.join(description_parts),
            'geneSymbol': genes[0],
            'biomarker_type': mut.norm_biomarker(evidence['Alteration type'], evidence['Biomarker'])
        })

    association = {}

    association['description'] = '{} {} {}'.format(' '.join(genes),
                                                   evidence['Drug full name'],
                                                   evidence['Association'])
    association['environmentalContexts'] = []
    association['environmentalContexts'].append({
        'description': evidence['Drug full name']})
    association['phenotypes'] = [{ 'description' : evidence['Primary Tumor type']} ]
    if not evidence['Metastatic Tumor Type'] == '':
        association['phenotypes'].append({ 'description' : evidence['Metastatic Tumor Type'] })

    pubs = []
    for p in evidence['Source'].split(';'):
        t = None
        if ':' in p:
            t, id = p.split(':')
        if t == 'PMID':
            pubs.append('http://www.ncbi.nlm.nih.gov/pubmed/{}'.format(id))
        else:
            pubs.append('https://www.google.com/#q={}'.format(p))

    association['evidence'] = [{
        "evidenceType": {
            "sourceName": "cgi"
        },
        'description': evidence['Association'],
        'info': {
            'publications': pubs
        }
    }]
    # add summary fields for Display

    association = el.evidence_label(evidence['Evidence level'], association)
    association = ed.evidence_direction(evidence['Association'], association)

    if "oncogenic" in evidence['Biomarker']:
        association['oncogenic'] = evidence['Biomarker']

    association['publication_url'] = pubs[0]
    association['drug_labels'] = evidence['Drug full name']
    feature_association = {'genes': genes,
                           'features': features,
                           'feature_names': evidence['Biomarker'],
                           'association': association,
                           'source': 'cgi',
                           'source_url': 'https://www.cancergenomeinterpreter.org/biomarkers',
                           'cgi': evidence}

    yield feature_association


def harvest(genes=None, drugs=None):
    """ get data from cgi """
    for evidence in _get_evidence(genes):
        yield evidence


def harvest_and_convert(genes=None, drugs=None):
    """ get data from cgi, convert it to ga4gh and return via yield """
    for evidence in harvest(genes, drugs):
        # print "harvester_yield {}".format(evidence.keys())
        # print evidence
        for feature_association in convert(evidence):
            # print "cgi convert_yield {}".format(feature_association.keys())
            yield feature_association


def _test():
    for feature_association in harvest_and_convert(['KIT']):
        logging.info(feature_association.keys())

if __name__ == '__main__':
    _test()
