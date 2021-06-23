#!/usr/bin/python
import logging
import re
import subprocess
from collections import defaultdict
from datetime import datetime
from io import StringIO
from pprint import pprint

import psycopg2.extras

from lookups.accession_mapping import NoMatchError, ensembl_txac_to_refseq
from lookups.hgvs_parsing import hgvsparser
from normalizers.gene_enricher import get_gene

# maps between chromosomes and refseq chromosome-level accessions
AC_MAP = {
        '1': 'NC_000001.10',
        '2': 'NC_000002.11',
        '3': 'NC_000003.11',
        '4': 'NC_000004.11',
        '5': 'NC_000005.9',
        '6': 'NC_000006.11',
        '7': 'NC_000007.13',
        '8': 'NC_000008.10',
        '9': 'NC_000009.11',
        '10': 'NC_000010.10',
        '11': 'NC_000011.9',
        '12': 'NC_000012.11',
        '13': 'NC_000013.10',
        '14': 'NC_000014.8',
        '15': 'NC_000015.9',
        '16': 'NC_000016.9',
        '17': 'NC_000017.10',
        '18': 'NC_000018.9',
        '19': 'NC_000019.9',
        '20': 'NC_000020.10',
        '21': 'NC_000021.8',
        '22': 'NC_000022.10',
        'X': 'NC_000023.10',
        '23': 'NC_000023.10',
        'Y': 'NC_000024.9',
    }

class SubmittedVariantTransformer(object):
    """
    Transforms sqlalchemy records into a VCF-compatible format.
    """
    VCF_HEADER_TEMPLATE = u"""##fileformat=VCFv4.0
##fileDate=%(curdate)s
##source=svip_queue
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"""

    @staticmethod
    def as_vcf(records, as_str=False):
        """
        For the given set of records, produces a bare-bones VCF representation.

        as_str: if true, returns a string rather than a file handle
        """

        fp =  StringIO()
        fp.write(SubmittedVariantTransformer.VCF_HEADER_TEMPLATE % {
            'curdate': unicode(datetime.now().strftime("%Y%m%d"))
        })
        fp.writelines(SubmittedVariantTransformer.as_vcf_row(rec) for rec in records)

        if as_str:
            result = fp.getvalue()
            fp.close()
            return result
        else:
            fp.seek(0)
            return fp

    @staticmethod
    def as_vcf_row(rec):
        # CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO
        # QUAL and INFO are both '.', which indicates an empty field(?)
        original_alt = ",".join([x.strip() for x in rec['alt'].strip("[]").replace("None", ".").split(",")])
        return u"\t".join(str(x) for x in [rec['chromosome'], rec['pos'], rec['id'], rec['ref'], original_alt, '.', 'PASS', '.']) + "\n"

def _parse_vep_output(data):
    lines = data.split("\n")
    _metadata = [x for x in lines if re.match(r'^##.*$', x)]
    header   = next((x.strip("#") for x in lines if re.match(r'^#[^#].*$', x)), None)
    records  = [x for x in lines if re.match(r'^[^#].*$', x)]

    if not header:
        raise Exception("VEP output lacks header line, aborting")
    else:
        header = header.split()

    for rec in records:
        cells = dict(zip(header, rec.split()))
        if 'Extra' in cells:
            try:
                cells['Extra'] = dict(x.split("=") for x in cells['Extra'].split(";"))
            except Exception as ex:
                print("Unparseable option line: %s" % cells['Extra'])
                raise ex
        yield cells

class InvalidVEPRecord(Exception):
    def __init__(self, message):
        self.message = message

def harvest(genes):
    """
    "Harvests" pending/errored entries from the SVIP submitted variants queue.

    Each variant is looked up in VEP, with the resulting data being passed to convert()
    to produce a structure that can be inserted into the silo (e.g., the postgres db).

    If 'genes' is specified, filters down the requested variants to just the list of HUGO
    gene symbols given in 'genes'. If 'genes' is None, no filtering is applied.
    """

    # get the postgres silo so we can query the db for submitted variants
    try:
        import harvester
        from silos.postgres_silo import PostgresSilo
        pg_silo = next(silo for silo in harvester.silos if isinstance(silo, PostgresSilo))
    except StopIteration:
        # if it doesn't exist, we can't proceed
        raise Exception("svip_queue harvester can't run without PostgresSilo being in the list of silos")

    with pg_silo._connect() as conn:
        variants_by_gene = defaultdict(list)
        gene_set = set(genes) if genes else None

        with conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as curs:
            curs.execute("select * from svip_submittedvariant where status != 'completed'")
            entries = curs.fetchall()

            logging.info("Processing %d non-complete records: %s" % (len(entries), [x['id'] for x in entries]))

            vcf_data = SubmittedVariantTransformer.as_vcf(entries, as_str=True)

            # TODO: produce set of variants as VCF, feed to VEP, gather output, parse, and pass on
            process = subprocess.Popen(["/app/scripts/run_curator_vep.sh"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            result, stderr_result = process.communicate(input=vcf_data)

            if process.returncode != 0:
                raise Exception("VEP failed w/code %d; stderr: %s" % (process.returncode, stderr_result))

            # build set of submission IDs which should retain only canonical variants
            # we'll use this for filtering below
            canonical_only_ids = set(int(x['id']) for x in entries if x['canonical_only'])

            variants = list(_parse_vep_output(result))

            # filter out non-canonical variants for submissions tha request canonical variants only
            variants = [
                x for x in variants
                if (
                    int(x['Uploaded_variation']) not in canonical_only_ids or
                    ('CANONICAL' in x['Extra'] and x['Extra']['CANONICAL'] == 'YES')
                )
            ]

            # check if any submitted entry wasn't listed
            variant_id_set = set(int(x['Uploaded_variation']) for x in variants)

            for entry in entries:
                if entry['id'] not in variant_id_set:
                    error_msg = 'VEP yielded no records for this entry'
                    logging.info("Skipping %d, reason: %s" % (entry['id'], error_msg))

                    # update the record as erroring out due to a missing symbol
                    curs.execute("""
                    update svip_submittedvariant
                    set processed_on=now(), status='error', error_msg=%s
                    where id=%s
                    """, (error_msg, entry['id']))

            for rec in variants:
                submit_id = int(rec['Uploaded_variation'])

                try:
                    try:
                        this_gene = rec['Extra']['SYMBOL']
                    except KeyError:
                        # skip variants for which no gene symbol was given (e.g., intronic variants)
                        raise InvalidVEPRecord('no gene symbol found in VEP record (consequence: %s)' % rec['Consequence'])

                    if gene_set and this_gene not in gene_set:
                        raise InvalidVEPRecord("skipped because it wasn't in the harvester's gene set")

                    # skip variants who can't be easily named by a protein change
                    if rec['Protein_position'] == '-' or rec['Amino_acids'] == '-':
                        raise  InvalidVEPRecord('protein position and/or amino acids were unspecified (consequence: %s)' % rec['Consequence'])

                except InvalidVEPRecord as ex:
                    logging.info("Skipping %d, reason: %s" % (submit_id, ex.message))
                    # pprint(rec)

                    # update the record as erroring out due to a missing symbol
                    curs.execute("""
                    update svip_submittedvariant
                    set processed_on=now(), status='error', error_msg=%s
                    where id=%s
                    """, (ex.message, submit_id,))
                    continue

                variants_by_gene[this_gene].append(rec)

    # start yielding outside of the postgres connection above, so that the transaction can be committed
    # *before* the postgres silo starts writing to svip_submittedvariant, avoiding a deadlock.
    for gene, variants in variants_by_gene.items():
        yield { 'gene': gene, 'variants': variants }


def _extract_name(variant):
    for part in variant['name'].split():
        if '-' not in part and not part == variant['entrez_name']:
            return part


def accession(hgvs_str):
    if not hgvs_str:
        return None
    return hgvs_str[:(hgvs_str.index(':'))]

def refseq_g_accession(hgvs_g_str):
    chromosome, rest = hgvs_g_str.split(":")
    return "%s:%s" % (AC_MAP[str(chromosome)], rest)

def extract_coordinates(hgvs_g_str):
    var = hgvsparser.parse_g_variant(refseq_g_accession(hgvs_g_str))
    return {
        'chromosome': hgvs_g_str.split(":")[0],
        'start': var.posedit.pos.start.base,
        'end': var.posedit.pos.end.base,
        'ref': var.posedit.edit.ref,
        'alt': var.posedit.edit.alt,
    }

def convert(gene_data):
    """given variant data from VEP, convert to ga4gh."""
    variants = gene_data['variants']

    # we can get some coarse location info, e.g. the chromosome, from the gene symbol itself
    # we'll retrieve that as a failover in case the civic entry is missing that info
    try:
        gene_meta = get_gene(gene_data['gene'])[0]
    except ValueError as ex:
        # this probably means the gene is missing, which means we can't really do anything...
        logging.warn(str(ex))
        return

    for variant in variants:
        coordinates = extract_coordinates(variant['Extra']['HGVSg'])
        gene_name = variant['Extra']['SYMBOL']
        amino_from, amino_to = variant['Amino_acids'].split("/")
        variant_name = "%s%s%s" % (amino_from, variant['Protein_position'], amino_to)

        # notes regarding 'isoform'/'refseq' fields being set to 'None' below:
        # neither field is actually currently used in the frontend, and it's always been a bit unclear
        # what they were supposed to be filled with in the ga4gh g2p schema. imo it's better for them
        # to be missing than incorrect (unless there's a motivating use case to bring them back).

        feature = {
            'geneSymbol': variant['Extra']['SYMBOL'],
            # 'entrez_id': variant['entrez_id'],
            'referenceName': 'GRCh37', # our VEP db is always GRCh37 (at least for now)
            'refseq': None, # accession(variant['Extra']['HGVSc']), # sometimes we get refseq, sometimes we get ensembl...
            'isoform': None, # it's unclear how to get the ensembl transcript identifier from vep, ironically...
            'chromosome': coordinates['chromosome'],
            'start': coordinates['start'],
            'end': coordinates['end'],
            'ref': coordinates['ref'],
            'alt': coordinates['alt'],
            'name': variant_name,
            'description': u'{} {}'.format(gene_name, variant_name),

            'hgvs_g': refseq_g_accession(variant['Extra']['HGVSg']),
            'hgvs_c': variant['Extra']['HGVSc'],
            'hgvs_p': variant['Extra']['HGVSp'],
            # hgvs_ensembl_c is also a possible field, but i don't know what we'd put in it for now

            'biomarker_type': variant['Consequence']
        }

        # if the refseq is still missing but isoform is specified, attempt to convert that into an NCBI refseq
        if not feature['refseq'] and feature['isoform']:
            try:
                feature['refseq'] = ensembl_txac_to_refseq(feature['isoform'])
            except NoMatchError as ex:
                logging.warn(ex)
                feature['refseq'] = None

        feat_assoc = {
            'genes': [gene_data['gene']],
            'features': [feature],
            'feature_names': 'n/a', # evidence_item['name'],
            'association': {
                'variant_name': variant_name,
                'description': '%s %s' % (gene_data['gene'], variant_name), # evidence_item['description'],
                # 'environmentalContexts': [],
                # 'phenotypes': [],
                # 'evidence': [],
                'extras': {
                    'submitted_var_id': int(variant['Uploaded_variation'])
                },
            },
            'source': 'svip_queue',
            'source_url': None
        }

        yield feat_assoc


def harvest_and_convert(genes):
    for gene_data in harvest(genes):
        for feat_assoc in convert(gene_data):
            yield feat_assoc
