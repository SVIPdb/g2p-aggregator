import csv

import requests
import tqdm

ensembl_refseq_map = {}


def ensembl_txac_to_refseq_ucsc():
    # uses UCSS genome browser export file to map ensembl transcript accessions to refseq transcripts
    with open("./ensembl_refseq_slim.tsv") as fp:
        eref_reader = csv.DictReader(fp, dialect='excel-tab')
        for row in tqdm.tqdm(eref_reader, desc="parsing ensembl-to-refseq table"):
            if not row['hg19.kgXref.refseq']:
                continue
            # FIXME: if there are multiple ensembl transcripts, only the last will be retained.
            #  we should somehow figure out which is most likely to be the right one...
            ensembl_refseq_map[row['#hg19.knownToEnsembl.value']] = {
                'refseq': row['hg19.kgXref.refseq'],
                'protAcc': row['hg19.kgXref.protAcc']
            }
    return ensembl_refseq_map


def ensembl_txac_to_refseq(ensembl_ac, loc='rna'):
    """
    Given an ensembl accession number (e.g., ENSG00000182533.6), searches mygene.info for a match and
    returns related NCBI refseq transcripts.
    :param ensembl_ac: the ensembl accession number
    :param loc: what kind of transcript to return: 'genomic' (NC_*), 'rna' (NM_*), or 'protein' (NP_*).
    The special value 'all' returns all three types in a dictionary.
    :return: depending on 'loc', a genomic, rna, or protein refseq accession number for the ensembl transcript, or
    a dictionary of {genomic, rna, protein} elements if loc is 'all'. None if no match is found.
    """
    loc_prefixes = {
        'genomic': 'NC_',
        'rna': 'NM_',
        'protein': 'NP_',
    }

    if loc != 'all' and loc not in loc_prefixes.keys():
        raise ValueError("loc must be 'all' or one of %s, got '%s' instead" % (loc_prefixes.keys(), loc))

    req = requests.get('http://mygene.info/v3/query', {
        'q': 'ensembl.transcript:%s' % ensembl_ac,
        'fields': 'ensembl,refseq'
    })

    resp = req.json()

    if 'total' not in resp or int(resp['total']) <= 0:
        return None

    if loc == 'all':
        return {
            k: [x for x in resp['hits'][0]['refseq'][k] if x.startswith(loc_prefixes[k])][0] for k in loc_prefixes
        }

    return [x for x in resp['hits'][0]['refseq'][loc] if x.startswith(loc_prefixes[loc])][0]
