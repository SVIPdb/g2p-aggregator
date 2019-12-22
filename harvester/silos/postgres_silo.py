from __future__ import print_function

import itertools
import json
import logging
import operator

import psycopg2
from psycopg2 import sql

from utils_ex.formatting import jsonify_or_none

VERBOSE = False


def populate_args(argparser):
    """add arguments we expect """
    argparser.add_argument('--pg_host',
                           help='''hostname for postgres server''',
                           default='localhost')
    argparser.add_argument('--pg_user',
                           help='''postgres username''',
                           default='postgres')
    argparser.add_argument('--pg_password',
                           help='''password for postgres user''',
                           default='postgres')
    argparser.add_argument('--pg_db',
                           help='''postgres database to populate''',
                           default='db')


class CustomCompare(object):
    def __init__(self, sql_stmt, value=None, values=None):
        self.sql = sql_stmt
        self.value = value
        self.values = values


def compose_params(params):
    """
    Produces a flat list of parameters from a mix of regular parameters and CustomCompare params. For CustomCompare
    entities, if they possess a 'value' only that is passed along, whereas if they have an iterable 'values' set
    that iterable is flatterened.
    :param params:
    :return:
    """
    for elem in params:
        if isinstance(elem, CustomCompare):
            if elem.value:
                yield elem.value
            else:
                for subelem in elem.values:
                    yield subelem
        else:
            yield elem


def array_to_sql(arr):
    """
    Converts an iterable to a postgres array.

    :arr the array to convert
    :return a SQL string suitable for embedding as a value
    """
    return sql.SQL("ARRAY[") + sql.SQL(", ").join(sql.Literal(v) for v in arr) + sql.SQL("]")


def _get_or_insert(curs, target_table, params, key_cols=None, append_cols=None, merge_existing=False):
    """
    Helper method for performing an idempotent insert; it always returns the key of the matched item,
    whether it already existed or had to insert it.
    :param curs: a database cursor
    :param target_table: the table into which to insert the values
    :param params: a set of values that uniquely identify an existing record
    :param key_cols: set of columns to compare for uniqueness before inserting (params is used if not specified)
    :param append_cols: set of jsonb columns to which to append values (set-like append, so duplicates are not recorded)
    :return: the ID of the existing or inserted record
    """

    if VERBOSE:
        logging.info("Attempting insert of: %s" % str(params))

    # create check and insert statements programmatically
    if key_cols is None:
        key_cols = params

    # this checks if there's an existing variant with columns/values specified in key_cols
    # we also check if it already has the sources we'd append via append_cols; if it does, we don't need to update it
    # the complete sql looks something like:
    #  "select id, (sources ? 'civic' and ...) as no_update_needed from api_variant where keyA=valA and..."
    # note: if the key_col is actually a SQL statement (e.g., if we need to do a custom comparison), use it verbatim.
    # otherwise formulate the comparison expression as 'colName=<placeholder>'
    cond = sql.SQL(' and ').join(map(
        lambda x: x[1].sql if isinstance(x[1], CustomCompare) else (sql.Identifier(x[0]) + sql.SQL("=") + sql.Placeholder()),
        key_cols.items()  # x[0] is the column name, x[1] is the value itself
    ))
    check_stmt = sql.SQL('select id, ({no_update_needed}) as no_update_needed from {target_tbl} where {key_compare}').format(
        no_update_needed=sql.SQL(' and ').join(map(
            lambda k: (
                # if it's not a list, checks if the scalar is already in the existing array, or...
                sql.SQL("{} = any({})").format(sql.Literal(append_cols[k]), sql.Identifier(k))
                if not isinstance(append_cols[k], (list, tuple)) else
                # ...checks if the list is a subset of the existing array
                sql.SQL("{}::text[] <@ {}").format(array_to_sql(append_cols[k]), sql.Identifier(k))
            ),
            append_cols.keys()
        )) if append_cols else sql.Literal(True),  # if we have no append_cols, we never need to update
        target_tbl=sql.Identifier(target_table),
        key_compare=cond
    )

    composed_params = list(compose_params(key_cols.values()))

    if VERBOSE:
        logging.info(check_stmt.as_string(curs))
        logging.info("(parameters: %s)" % composed_params)

    # first check if it exists and get its id if so...
    # (note that CustomCompare elements have their value embedded in the object instead of just being the value)
    curs.execute(check_stmt, composed_params)
    result = curs.fetchone()

    if result is not None:
        # TODO: perhaps we merge only if we're adding from a new source, since ostensibly
        #  the same source will have the same info for genes/variants across evidence items.
        #  that would require a bit of coupling to the table spec, specifically checking the sources
        #  before flagging for a merge

        # it does exist, so retrieve its identifying info
        if VERBOSE:
            logging.info("result: %s" % str(result))

        entry_id = result[0]
        update_required = not result[1]

        # and if we have append_cols that need to be updated, do so
        if (append_cols is not None and update_required) or merge_existing:
            # construct the 'append' part of the update if there are set fields that we need to append to
            # we assume our append_cols are arrays values that we're treating as sets (via array_distinct)
            # jsonb_set(a, b, 'null', TRUE) results in a's keys being merged with b's.
            append_part = (
                map(
                    lambda x: (
                        sql.SQL("{}=array_distinct({} || {})").format(
                            sql.Identifier(x),
                            sql.Identifier(x),
                            sql.SQL("ARRAY[") + (sql.Identifier(append_cols[x])) + sql.SQL("]")
                            if not isinstance(append_cols[x], (list, tuple)) else
                            array_to_sql(append_cols[x])
                        )
                    ),
                    append_cols.keys()
                )
            ) if append_cols is not None and update_required else []

            # construct the 'merge' part of the query if merge_existing is true
            # FIXME: arrays should probably be extended instead of having their contents replaced
            if merge_existing:
                merge_items = [p for p in params.items() if p[0] not in key_cols.keys()]
                merge_part = (
                    map(
                        lambda x: (
                            sql.SQL("{}=coalesce({}, {})").format(
                                sql.Identifier(x),
                                sql.Identifier(x),
                                sql.Placeholder()
                            )
                        ),
                        [p[0] for p in merge_items]
                    )
                )
            else:
                merge_items = None
                merge_part = []

            # finally perform the update, which may append to sets or merge fields
            update_stmt = sql.SQL('update {} set {} where id={}').format(
                sql.Identifier(target_table),
                sql.SQL(',').join(append_part + merge_part),
                sql.Literal(entry_id)
            )

            if VERBOSE:
                logging.info(update_stmt.as_string(curs))

            if merge_existing:
                # we fill out the merge columns with new values coalesced with the old ones
                curs.execute(update_stmt, [p[1] for p in merge_items])
            else:
                # we're only appending to columns, in which case the appended items are embedded in the query
                curs.execute(update_stmt)

    else:
        # ...it doesn't exist, so insert it and get its resulting id

        # add in the append columns to the params since we're creating this thing now
        if append_cols:
            for col in append_cols:
                params[col] = '{%s}' % append_cols[col]

        insert_stmt = sql.SQL('insert into {} ({}) values ({}) returning id').format(
            sql.Identifier(target_table),
            sql.SQL(', ').join(map(sql.Identifier, params.keys())),
            sql.SQL(', ').join(sql.Placeholder() * len(params))
        )

        try:
            if VERBOSE:
                logging.info(insert_stmt.as_string(curs))
            curs.execute(insert_stmt, params.values())
            entry_id = curs.fetchone()[0]
        except psycopg2.IntegrityError as ex:
            print("Failed to insert into %s; statement: %s" % (target_table, insert_stmt.as_string(curs)))
            raise ex

    return entry_id


class PostgresSilo:
    """
    Inserts assocation entries into a postgresql database.

    Specifically, breaks down each association entry into genes, variants, and various bits of context and evidence.
    Creates the genes and variants if they don't already exist, then creates an association object keyed to the
    gene/variant. Attaches evidence, environmental context, and phenotype information, if available, to the newly
    created association.
    """

    def __init__(self, args):
        self._pg_host = args.pg_host
        self._pg_db = args.pg_db
        self._pg_user = args.pg_user
        self._pg_password = args.pg_password

        # caches source name to id, populated in save_bulk()
        self.source_to_id = {}

    def __str__(self):
        return "PostgresSilo host:{} db:{}".format(self._pg_host, self._pg_db)

    def _connect(self):
        return psycopg2.connect("host='%s' dbname='%s' user='%s' password='%s'" % (
            self._pg_host, self._pg_db, self._pg_user, self._pg_password
        ))

    def delete_all(self):
        """delete index"""
        try:
            logging.info("PostgresSilo: performing full clear of db")

            # drop everything that we'd be inserting and reset the IDs
            with self._connect() as conn:
                with conn.cursor() as curs:
                    # if you delete all the genes the cascading delete clears everything in the database,
                    # since everything is keyed to a gene at some point
                    curs.execute("delete from api_gene")
                    curs.execute("alter sequence api_environmentalcontext_id_seq restart with 1")
                    curs.execute("alter sequence api_evidence_id_seq restart with 1")
                    curs.execute("alter sequence api_phenotype_id_seq restart with 1")
                    curs.execute("alter sequence api_association_id_seq restart with 1")
                    curs.execute("alter sequence api_variantinsource_id_seq restart with 1")
                    curs.execute("alter sequence api_svipvariant_id_seq restart with 1")
                    curs.execute("alter sequence api_variant_id_seq restart with 1")
                    curs.execute("alter sequence api_gene_id_seq restart with 1")
                conn.commit()

        except Exception as e:
            logging.info("PostgresSilo: delete_all failed: {}".format(e))
            raise e

    def delete_source(self, source):
        """ delete source from index """
        try:
            logging.info("PostgresSilo: clearing records which mention source '%s'..." % source)

            with self._connect() as conn:
                with conn.cursor() as curs:
                    curs.execute(
                        """
                        WITH vis_w_source AS (
                            select vis.id, src.name as name
                            from api_variantinsource vis
                            inner join api_source src on vis.source_id = src.id
                            where src.name=%s
                        )
                        DELETE FROM api_variantinsource B
                        USING vis_w_source C
                        WHERE B.id=C.id
                        """, (source,)
                    )
                    # after removing the source, and thus variant-in-source entries, remove any orphaned variants
                    # (i.e., variants with no remaining supporting entries...)
                    # FIXME: eventually we'll have SVIP variants that may or may not have entries in external databases.
                    #  we'll either need to represent SVIP as another entry in api_source/api_variantinsource, or
                    #  simply allow variants with no public supporting evidence to continue to exist.
                    curs.execute(
                        """delete from api_variant v
                        where not exists (select * from api_variantinsource as avs where avs.variant_id=v.id)"""
                    )
                conn.commit()
            logging.info("...completed clearing for source {}".format(source))

        except Exception as e:
            logging.info("PostgresSilo: delete_source failed: {}".format(e))

    def _save_one(self, curs, feature_association):
        # caveats:
        # 1) this code should be capable of ensuring that the hierarchy can be built incrementally,
        #    e.g., we should insert a previously unseen gene first, then hook up the variant to that new instance
        #    (note 1a: b/c of that, we encounter a variant once per *evidence item*; we might reconsider updating
        #     the variant's information with each evidence item encounter, since that could get ridiculous...)
        #    (note 1b: perhaps we should merge variants upstream and pass this sets of evidence items per variant?)
        #    (note 1c: maybe we only merge a variant's info if we haven't previously merged that source?)
        # 2) we're splitting the dict into multiple tables, so some inserts might succeed even if the whole dict
        #    can't be inserted. do we fail cleanly and revert all changes in that case, or do we just try to insert
        #    what we can?

        # pre-step: ensure the feature_association has everything we need
        if 'gene_identifiers' not in feature_association:
            logging.warn("Feature association lacks essential key 'gene_identifiers', skipping")
            return

        # --------------------------------------------------------------------------------
        # stage 1. insert genes mentioned in the payload if they don't already exist
        # the result of this will be a genes_to_ids lookup table that we can use to replace text references w/db IDs
        # --------------------------------------------------------------------------------

        genes_to_ids = {}
        for gene in feature_association['gene_identifiers']:
            gene_id = _get_or_insert(curs, "api_gene", {
                'entrez_id': int(gene['entrez_id']),
                'ensembl_gene_id': gene['ensembl_gene_id'],
                'location': gene['location'],
                'symbol': gene['symbol'],
            }, append_cols={
                'sources': feature_association['source'],
                'uniprot_ids': gene['uniprot_ids'],
                'aliases': gene['aliases'],
                'prev_symbols': gene['prev_symbols'],
            })
            genes_to_ids[unicode(gene['symbol'])] = gene_id

            if VERBOSE:
                print("Got gene %s w/id: %d" % (gene['symbol'], gene_id))

        # --------------------------------------------------------------------------------
        # stage 2. insert variants mentioned in the features
        # the result of this will be a variants_to_ids lookup table that we can use to replace text references w/db IDs
        # TODO: check if we ever actually get more than one feature per piece of evidence
        # --------------------------------------------------------------------------------

        # FIXME: we need to decide how to deal with older/different genes showing up in references.
        #  do we just link them up to the aliased gene, do we create a new entry, what...?
        def resolve_gene_id(symbol):
            if symbol in genes_to_ids:
                return genes_to_ids[symbol]
            elif len(genes_to_ids) == 1:
                return genes_to_ids.values()[0]
            else:
                raise Exception("geneSymbol %s can't be found in genes_to_ids (%s)" % (symbol, ', '.join(genes_to_ids.keys())) )

        variants_to_ids = {}
        last_variant_id = None
        for feat in feature_association['features']:
            this_gene_id = resolve_gene_id(feat['geneSymbol'])

            # insert a variant keyed to this gene_id
            var_obj = {
                # these two are used as 'key cols' that we use to match variants between sources
                'gene_id': this_gene_id,
                'name': feat['name'],

                'description': feat.get('description'),
                'reference_name': feat.get('referenceName'),
                'refseq': feat.get('refseq'),
                'isoform': feat.get('isoform'),
                'biomarker_type': feat.get('biomarker_type'),

                # position/change info
                'chromosome': feat.get('chromosome'),
                'start_pos': feat.get('start'),
                'end_pos': feat.get('end'),
                'ref': feat.get('ref'),
                'alt': feat.get('alt'),

                # hgvs coordinates
                'hgvs_g': feat.get('hgvs_g'),
                'hgvs_c': feat.get('hgvs_c'),
                'hgvs_p': feat.get('hgvs_p'),

                'dbsnp_ids': feat.get('dbsnp_ids'),
                'myvariant_hg19': feat.get('myvariant_hg19'),

                'mv_info': jsonify_or_none(feat.get('mv_info')),  # optional info from the myvariant_enricher "normalizer",

                'crawl_status': jsonify_or_none(feat.get('crawl_status'))
            }

            if 'sequence_ontology' in feat:
                # in nearly all cases this is included, but there's a civic entry that's missing it for some reason
                seq_ont = feat['sequence_ontology']
                var_obj.update({
                    'so_hierarchy': seq_ont.get('hierarchy'),
                    'soid': seq_ont['soid'],
                    'so_name': seq_ont['name']
                })

            # FIXME: different data sources have different ways of referencing the same variant, sadly
            # we should be able to support the following disambiguations:
            # 1. the same protein-level change (e.g., V600E and Val600Glu or p.Val600Glu)
            # 2. the same SNP (e.g. c.1779T>A)
            # 3. hgvs string matches without reference sequences

            variant_id = _get_or_insert(
                curs, "api_variant", var_obj,
                key_cols={
                    'gene_id': this_gene_id,
                    'name': CustomCompare(
                        sql_stmt=sql.SQL("lower({})=lower({})").format(sql.Identifier('name'), sql.Placeholder()),
                        value=feat['name']
                    ),
                    'hgvs_g': CustomCompare(
                        sql_stmt=sql.SQL("(({id} is null and {val} is null) or {id}={val})").format(
                            id=sql.Identifier('hgvs_g'), val=sql.Placeholder()
                        ),
                        # duplicated b/c SQL.format() ignores keyed substitutions
                        values=[feat.get('hgvs_g'), feat.get('hgvs_g')]
                    ),
                    # 'biomarker_type': var_obj['biomarker_type']
                    # FIXME: the above might identify differences in punctuation or case as different variants...
                    # ^ this is in fact the case, e.g. "missense_variant" (oncokb) vs. "Missense Variant" (civic)
                }, append_cols={'sources': feature_association['source']},
                merge_existing=True
            )

            variants_to_ids[feat['name']] = variant_id
            last_variant_id = variant_id  # keep track of the most-recently-inserted variant (hopefully the only one...)

            if VERBOSE:
                print("Got variant_id: %d" % variant_id)

        # --------------------------------------------------------------------------------
        # stage 3. insert association object that relates variant to evidence
        # results in an ID that we can use to key the evidence to the variant
        # --------------------------------------------------------------------------------

        this_gene_id = resolve_gene_id(feature_association['genes'][0])
        source_id = self.source_to_id[feature_association['source']]

        assoc = feature_association['association']

        # --- stage 3a. look up a variant-to-source mapping, if it exists, or create it if it doesn't
        variant_in_source_id = _get_or_insert(
            curs, "api_variantinsource", {
                'source_id': source_id,
                'variant_id': last_variant_id,
                'variant_url': feature_association.get('source_url'),
                'extras': jsonify_or_none(assoc.get('extras'))
            },
            key_cols={
                'source_id': source_id,
                'variant_id': last_variant_id
            },
            merge_existing=True
        )

        association_obj = {
            'source': feature_association['source'],
            'payload': json.dumps(feature_association),
            'description': assoc.get('description'),
            'drug_labels': assoc.get('drug_labels'),
            'drug_interaction_type': assoc.get('drug_interaction_type'),
            'variant_name': assoc.get('variant_name'),
            'source_link': assoc.get('source_link'),  # this is the actual sample ID
            'variant_in_source_id': variant_in_source_id,

            'evidence_type': assoc.get('evidence_type'),
            'evidence_direction': assoc.get('evidence_direction'),
            'clinical_significance': assoc.get('clinical_significance'),
            'evidence_level': assoc.get('evidence_level'),

            'crawl_status': jsonify_or_none(assoc.get('crawl_status')),
            'extras': jsonify_or_none(assoc.get('extras'))
        }

        curs.execute(
            """
            insert into api_association (%s) values (%s) returning id
            """ % (", ".join(association_obj.keys()), ", ".join(["%s"] * len(association_obj))),
            association_obj.values()
        )

        # get the association object id so we can tie any phenotypes, evidence entries, and env. contexts to it
        assoc_id = curs.fetchone()[0]

        # --------------------------------------------------------------------------------
        # stage 4. insert phenotypes, evidence items, and contexts that are associated with the variant
        # --------------------------------------------------------------------------------

        if 'phenotypes' in assoc:
            for pheno in assoc['phenotypes']:
                curs.execute(
                    """
                    insert into api_phenotype (source, term, pheno_id, family, description, association_id)
                    values (%s, %s, %s, %s, %s, %s)
                    """,
                    (
                        pheno.get('source'),
                        pheno.get('term'),
                        pheno.get('id'),
                        pheno.get('family'),
                        pheno.get('description'),
                        assoc_id
                    )
                )

        if 'evidence' in assoc:
            for evidence in assoc['evidence']:
                evidence_obj = {
                    'publications': evidence['info'].get('publications') if 'info' in evidence else None,
                    'evidenceType_sourceName': evidence['evidenceType'].get('sourceName') if 'evidenceType' in evidence else None,
                    'evidenceType_id': evidence['evidenceType'].get('id') if 'evidenceType' in evidence else None,
                    'association_id': assoc_id
                }

                curs.execute(
                    """
                    insert into api_evidence (%s) values (%s)
                    """ % (", ".join('"%s"' % k for k in evidence_obj.keys()), ", ".join(["%s"] * len(evidence_obj))),
                    evidence_obj.values()
                )

        if 'environmentalContexts' in assoc:
            for env_context in assoc['environmentalContexts']:
                curs.execute(
                    """
                    insert into api_environmentalcontext (source, term, envcontext_id, usan_stem, description, association_id)
                    values (%s, %s, %s, %s, %s, %s)
                    """,
                    (
                        env_context.get('source'),
                        env_context.get('term'),
                        env_context.get('id'),
                        env_context.get('usan_stem'),
                        env_context.get('description'),
                        assoc_id
                    )
                )

    def save_bulk(self, feature_association_generator):
        """ write many feature_associations to the db """

        conn = self._connect()

        with conn.cursor() as curs:
            curs.execute("select * from api_source")
            self.source_to_id = dict((x[1], x[0]) for x in curs.fetchall())

            if len(self.source_to_id) == 0:
                # FIXME: ideally we should ensure this elsewhere, but for now let's make sure we have some sources
                curs.execute("INSERT INTO public.api_source (id, name, display_name) VALUES (1, 'civic', 'CIViC')")
                curs.execute("INSERT INTO public.api_source (id, name, display_name) VALUES (2, 'oncokb', 'OncoKB')")
                curs.execute("INSERT INTO public.api_source (id, name, display_name) VALUES (3, 'clinvar', 'ClinVar')")
                curs.execute("INSERT INTO public.api_source (id, name, display_name) VALUES (4, 'cosmic', 'COSMIC')")

                print("Inserted civic, oncokb, clinvar, cosmic into sources, since there weren't any...")


        # split transactions by source (for now)
        for source, source_feats in itertools.groupby(feature_association_generator, key=operator.itemgetter('source')):
            print("=> Getting source %s..." % source)

            for gene, gene_feats in itertools.groupby(source_feats, key=operator.itemgetter('genes')):
                print("=> Processing gene %s..." % gene)

                total_inserted = 0
                skipped = 0

                with conn:
                    with conn.cursor() as curs:
                        # FIXME: consider chunking into multiple transactions so as not to overload the trans. buffer
                        # we may also run harvesters in parallel, in which case we may prefer less contention between
                        # long-running transactions
                        for feature_association in gene_feats:
                            total_inserted += 1
                            try:
                                self._save_one(curs, feature_association)
                            except Exception as ex:
                                skipped += 1
                                continue

                    logging.info("Inserted %d entries for gene %s, skipped %d" % (total_inserted, gene, skipped))

                    conn.commit()

    def save(self, feature_association):
        """ write dict to a series of tables"""

        with self._connect() as conn:
            with conn.cursor() as curs:
                self._save_one(curs, feature_association)
            conn.commit()
