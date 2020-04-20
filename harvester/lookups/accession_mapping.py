import csv
import itertools
import logging
import os
import sys
import operator as op

import requests
import tqdm
import MySQLdb

from utils_ex.iterables import matched

ensembl_refseq_map = {}


# used for ensembl_txac_to_refseq_ensembldb() to convert ENST IDs to NCBI refseq
def connect_to_ensembldb():
    return MySQLdb.connect(host="ensembldb.ensembl.org", user="anonymous", db="homo_sapiens_core_75_37")


conn = connect_to_ensembldb()


class NoMatchError(RuntimeError):
    pass


def ensembl_txac_to_refseq_ucsc():
    # uses UCSS genome browser export file to map ensembl transcript accessions to refseq transcripts
    with open("/data/ensembl_refseq_slim.tsv") as fp:
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


def ensembl_txac_to_refseq_mygene(ensembl_ac, loc='rna'):
    """
    Given an ensembl accession number (e.g., ENSG00000182533.6), searches mygene.info for a match and
    returns related NCBI refseq transcripts.

    Note that this method unfortunately disagrees with ensembldb, e.g. 'ENST00000288602' produces 'NM_001354609.1',
    whereas ensembldb produces 'NM_004333.4'.

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


def ensembl_txac_to_refseq_ensembldb(ensembl_ac, use_version=True, retries=5, strict_match=False):
    """
    Returns the corresponding NCBI RefSeq accession ID for the given ensemble transcript accession ID, if it exists.

    This method queries the MySQL server at ensembldb.org directly every time.

    >>> ensembl_txac_to_refseq_ensembldb("ENST00000216024")
    'NM_007068.3'
    >>> ensembl_txac_to_refseq_ensembldb("ENST00000206451.6")  # if multiple items exist, the first is returned
    'NM_006263.3'

    :param ensembl_ac: the ensembl transcript accession ID, e.g. "ENST00000216024.2" (the version, i.e., ".2", is optional)
    :param use_version: if true, uses the version given in the ensembl accession ID, if it was specified
    :param retries: the number of times to re-attempt the query before aborting
    :param strict_match: if True, throws an exception if there is not exactly one match
    :return: the matching refseq ID if found, else None
    """

    global conn
    first_pass = True

    while first_pass or retries > 0:
        first_pass = False

        try:
            with conn.cursor() as c:
                # query from https://www.biostars.org/p/106470/#106560 (w/added filter)

                if use_version and '.' in ensembl_ac:
                    ensembl_ac_stable, version = ensembl_ac.split('.')
                else:
                    ensembl_ac_stable, version = None, None

                query = """
                    SELECT transcript.stable_id, xref.display_label
                    FROM transcript, object_xref, xref, external_db
                    WHERE transcript.transcript_id = object_xref.ensembl_id
                     AND object_xref.ensembl_object_type = 'Transcript'
                     AND object_xref.xref_id = xref.xref_id
                     AND xref.external_db_id = external_db.external_db_id
                     AND external_db.db_name = 'RefSeq_mRNA'
                     AND transcript.stable_id=%s
                 """

                if use_version and ensembl_ac_stable is not None and version is not None:
                    numrows = c.execute(query + " AND transcript.version=%s", (ensembl_ac_stable, version))
                else:
                    numrows = c.execute(query, (ensembl_ac,))

                if strict_match and numrows != 1:
                    raise NoMatchError("Expected one row to map to accession %s, got %d" % (ensembl_ac, numrows))

                results = c.fetchone()
                return results[1]  # we return the second field, i.e. display_label

        except MySQLdb.OperationalError as ex:
            # attempt to reconnect if we have any left
            if retries > 0:
                retries -= 1
                conn.close()
                conn = connect_to_ensembldb()
            else:
                raise NoMatchError("Aborted after too many retry attempts resulted in MySQLdb.OperationalError"), None, sys.exc_info()[2]


_ensembl_txac_cache = None
ENSEMBL_TXAC_FILE = "/data/ensembl_txac_to_refseq.tsv"


# noinspection PyTypeChecker
def _cached_ensembl_txac_to_refseq(ensembl_ac, use_version=True):
    """
    Returns the corresponding NCBI RefSeq accession ID for the given ensemble transcript accession ID, if it exists.

    This method uses an in-memory cache, which is created from a file of ensembl-to-refseq mappings the first time it is
    run. If that file doesn't exist, it is retrieved from the MySQL server at ensembldb.org. If *that* fails, then the
    method bails with a RuntimeError.

    >>> _cached_ensembl_txac_to_refseq("ENST00000216024")
    'NM_007068.3'
    >>> _cached_ensembl_txac_to_refseq("ENST00000206451.6")  # if multiple items exist, the first is returned
    'NM_006263.3'

    :param ensembl_ac: the ensembl transcript accession ID, e.g. "ENST00000216024.2" (the version, i.e., ".2", is optional)
    :param use_version: if true, uses the version given in the ensembl accession ID, if it was specified
    :return: the matching refseq ID if found, else None
    """

    # attempt to load resource file from disk
    global _ensembl_txac_cache

    # ensure that cache exists; if it doesn't, create it
    if not _ensembl_txac_cache:
        # check that the ensembl_txac_to_refseq file exists, creating it from ensembldb.org if not
        if not os.path.exists(ENSEMBL_TXAC_FILE):
            global conn

            curpath = os.path.abspath(os.curdir)
            print "Current path is: %s" % curpath

            with conn.cursor() as c:
                retries = 5
                logging.info("Creating ensembl transcript lookup file from source...")

                while retries > 0:
                    try:
                        with open(ENSEMBL_TXAC_FILE, "w") as out_fp:
                            c.execute("""
                                SELECT transcript.stable_id, transcript.version, xref.display_label
                                FROM transcript, object_xref, xref, external_db
                                WHERE transcript.transcript_id = object_xref.ensembl_id
                                 AND object_xref.ensembl_object_type = 'Transcript'
                                 AND object_xref.xref_id = xref.xref_id
                                 AND xref.external_db_id = external_db.external_db_id
                                 AND external_db.db_name = 'RefSeq_mRNA'
                                order by transcript.stable_id, transcript.version;
                            """)

                            results = c.fetchall()
                            out_fp.writelines(["stable_id\tversion\tdisplay_label\n"])
                            out_fp.writelines("\t".join(str(c) for c in x) + "\n" for x in results)
                            break

                    except MySQLdb.OperationalError, ex:
                        retries -= 1
                        conn.close()
                        conn = connect_to_ensembldb()

                if retries > 0:
                    logging.info("done!")
                else:
                    raise RuntimeError("Couldn't contact remote server to get transcript mapping file!")

        # at this point we have the cache file, so now we just need to unpack it into memory
        with open(ENSEMBL_TXAC_FILE) as fp:
            logging.info("Loading cache file from %s..." % ENSEMBL_TXAC_FILE)

            header = next(fp).strip().split('\t')
            reader = csv.DictReader(fp, fieldnames=header, dialect='excel-tab')

            _ensembl_txac_cache = dict(
                (k, [{'refseq': x['display_label'], 'version': x['version']} for x in group])
                for k, group in itertools.groupby(reader, op.itemgetter('stable_id'))
            )

    try:
        if '.' in ensembl_ac:
            ensembl_ac_stable, version = ensembl_ac.split('.')
        else:
            ensembl_ac_stable, version = ensembl_ac, None

        if use_version and version:
            return matched(_ensembl_txac_cache[ensembl_ac_stable], lambda z: z['version'] == version)['refseq']
        else:
            # just get the first one
            return _ensembl_txac_cache[ensembl_ac_stable][0]['refseq']
    except (KeyError, TypeError):
        return None


# we export ensembl_txac_to_refseq for use elsewhere, but can swap the underlying method
ensembl_txac_to_refseq = _cached_ensembl_txac_to_refseq


if __name__ == "__main__":
    import doctest
    doctest.testmod()
