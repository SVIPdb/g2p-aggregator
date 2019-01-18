from __future__ import print_function

import itertools
import json
import logging
import operator

import psycopg2
from psycopg2 import sql


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


def _get_or_insert(curs, target_table, params, key_cols=None, append_cols=None):
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

    # initially we don't know the gene_id, or even if it exists yet
    gene_id = None

    # the call to every() below determines if there's a single column in the appends that doesn't have the source
    # already set. if there is, we have to peform a full update to all append columns
    cond = sql.SQL(' and ').join(map(lambda x: sql.Identifier(x) + sql.SQL("=") + sql.Placeholder(), key_cols.keys()))
    check_stmt = sql.SQL('select id, ({}) as no_update_needed from {} where {}').format(
        sql.SQL(' and ').join(map(
            lambda k: sql.SQL("{} ? {}").format(sql.Identifier(k), sql.Literal(append_cols[k])),
            append_cols.keys()
        )) if append_cols else "true", # if we have no append_cols, we never need to update
        sql.Identifier(target_table),
        cond
    )

    if VERBOSE:
        logging.info(check_stmt.as_string(curs))

    # first check if it exists and get its id if so...
    curs.execute(check_stmt, key_cols.values())
    result = curs.fetchone()

    if result is not None:
        if VERBOSE:
            logging.info("result: %s" % str(result))

        gene_id = result[0]
        update_required = not result[1]

        if append_cols is not None and update_required:
            # it did exist, so we need to update its append_cols columns with whatever we're supposed to append to them
            # update mytable set somecol=somecol||{"newval":null} effectively appends newval to the set somecol
            update_stmt = sql.SQL('update {} set {} where id={}').format(
                sql.Identifier(target_table),
                sql.SQL(',').join(map(
                    lambda x: (
                        sql.SQL("{}=jsonb_set({}, '{}', 'null', TRUE)").format(
                            sql.Identifier(x),
                            sql.Identifier(x),
                            sql.SQL("{") + sql.SQL(append_cols[x]) + sql.SQL("}")
                        )
                    ),
                    append_cols.keys()
                )),
                sql.Literal(gene_id)
            )
            if VERBOSE:
                logging.info(update_stmt.as_string(curs))
            curs.execute(update_stmt)

    else:
        # ...it doesn't exist, so insert it and get its resulting id
        # add in the append columns to the params since we're creating this thing now
        for col in append_cols:
            params[col] = '{"%s":null}' % append_cols[col]

        insert_stmt = sql.SQL('insert into {} ({}) values ({}) returning id').format(
            sql.Identifier(target_table),
            sql.SQL(', ').join(map(sql.Identifier, params.keys())),
            sql.SQL(', ').join(sql.Placeholder() * len(params))
        )
        if VERBOSE:
            logging.info(insert_stmt.as_string(curs))
        curs.execute(insert_stmt, params.values())
        gene_id = curs.fetchone()[0]

    return gene_id


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
                    curs.execute("delete from api_environmentalcontext")
                    curs.execute("delete from api_evidence")
                    curs.execute("delete from api_phenotype")
                    curs.execute("delete from api_association")
                    curs.execute("delete from api_variant")
                    curs.execute("delete from api_gene")
                    curs.execute("alter sequence api_environmentalcontext_id_seq restart with 1")
                    curs.execute("alter sequence api_evidence_id_seq restart with 1")
                    curs.execute("alter sequence api_phenotype_id_seq restart with 1")
                    curs.execute("alter sequence api_association_id_seq restart with 1")
                    curs.execute("alter sequence api_variant_id_seq restart with 1")
                    curs.execute("alter sequence api_gene_id_seq restart with 1")
                conn.commit()

        except Exception as e:
            logging.info("PostgresSilo: delete_all failed: {}".format(e))

    def delete_source(self, source):
        """ delete source from index """
        try:
            logging.info("PostgresSilo: clearing records which mention source '%s'..." % source)

            with self._connect() as conn:
                with conn.cursor() as curs:
                    curs.execute("delete from api_environmentalcontext")
                    curs.execute("delete from api_evidence")
                    curs.execute("delete from api_phenotype")
                    curs.execute("delete from api_variant where sources ? %s", (source,))
                    # curs.execute("delete from api_gene where sources ? %s", (source,))
                conn.commit()
            logging.info("...completed clearing for source {}".format(source))

        except Exception as e:
            logging.info("PostgresSilo: delete_source failed: {}".format(e))

    @staticmethod
    def _save_one(curs, feature_assocation):
        # caveats:
        # 1) this code should be capable of ensuring that the hierarchy can be built incrementally,
        #    e.g., we should insert a previously unseen gene first, then hook up the variant to that new instance
        # 2) we're splitting the dict into multiple tables, so some inserts might succeed even if the whole dict
        #    can't be inserted. do we fail cleanly and revert all changes in that case, or do we just try to insert
        #    what we can?

        # --------------------------------------------------------------------------------
        # stage 1. insert genes mentioned in the payload if they don't already exist
        # the result of this will be a genes_to_ids lookup table that we can use to replace text references w/db IDs
        # --------------------------------------------------------------------------------

        genes_to_ids = {}
        for gene in feature_assocation['gene_identifiers']:
            gene_id = _get_or_insert(curs, "api_gene", {
                'entrez_id': int(gene['entrez_id']),
                'ensembl_gene_id': gene['ensembl_gene_id'],
                'uniprot_ids': gene['uniprot_ids'],
                'location': gene['location'],
                'symbol': gene['symbol']
            }, append_cols={'sources': feature_assocation['source']})
            genes_to_ids[gene['symbol']] = gene_id

            if VERBOSE:
                print("Got gene %s w/id: %d" % (gene['symbol'], gene_id))

        # --------------------------------------------------------------------------------
        # stage 2. insert variants mentioned in the features
        # the result of this will be a variants_to_ids lookup table that we can use to replace text references w/db IDs
        # TODO: check if we ever actually get more than one feature per piece of evidence
        # --------------------------------------------------------------------------------

        variants_to_ids = {}
        last_variant_id = None
        for feat in feature_assocation['features']:
            # insert a variant keyed to this gene_id
            var_obj = {
                'gene_id': genes_to_ids[feat['geneSymbol']],
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
                'hgvs_g': feat.get('hgvs_g')
            }

            if 'sequence_ontology' in feat:
                # in nearly all cases this is included, but there'sa civic entry that's missing it for some reason
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
                    'gene_id': genes_to_ids[feat['geneSymbol']],
                    'name': feat['name'],
                    # 'biomarker_type': var_obj['biomarker_type']
                    # FIXME: the above might identify differences in punctuation or case as different variants...
                    # ^ this is in fact the case, e.g. "missense_variant" (oncokb) vs. "Missense Variant" (civic)
                }, append_cols={'sources': feature_assocation['source']}
            )

            variants_to_ids[feat['name']] = variant_id
            last_variant_id = variant_id  # keep track of the most-recently-inserted variant (hopefully the only one...)

            if VERBOSE:
                print("Got variant_id: %d" % variant_id)

        # --------------------------------------------------------------------------------
        # stage 3. insert association object that relates variant to evidence
        # results in an ID that we can use to key the evidence to the variant
        # --------------------------------------------------------------------------------

        assoc = feature_assocation['association']

        curs.execute(
            """
            insert into api_association
            (source, source_url, payload, description, drug_labels, variant_name, source_link, evidence_label, response_type, evidence_level, gene_id, variant_id)
            values (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s) returning id
            """,
            (
                feature_assocation['source'],
                feature_assocation['source_url'],
                json.dumps(feature_assocation),
                assoc.get('description'),
                assoc.get('drug_labels'),
                assoc.get('variant_name'),
                assoc.get('source_link'),
                assoc.get('evidence_label'),
                assoc.get('response_type'),
                assoc.get('evidence_level'),
                genes_to_ids[feature_assocation['genes'][0]],
                last_variant_id
            )
        )

        # get the association object id so we can tie any phenotypes, evidence entries, and env. contexts to it
        assoc_id = curs.fetchone()[0]

        # we can't continue unless we created an association object
        # FIXME: replace with a proper exception and handler later
        assert assoc_id is not None

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
                curs.execute(
                    """
                    insert into api_evidence
                    (type, description, publications, "evidenceType_sourceName", "evidenceType_id", association_id)
                    values (%s, %s, %s, %s, %s, %s)
                    """,
                    (
                        evidence.get('type', 'unknown'),  # one of predictive, diagnostic, prognostic, or predisposing; we need to determine this from the payload somehow
                        evidence.get('description'),
                        evidence['info'].get('publications') if 'info' in evidence else None,
                        evidence['evidenceType'].get('sourceName') if 'evidenceType' in evidence else None,
                        evidence['evidenceType'].get('id') if 'evidenceType' in evidence else None,
                        assoc_id
                    )
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

        # split transactions by source (for now)
        for source, source_feats in itertools.groupby(feature_association_generator, key=operator.itemgetter('source')):
            print("=> Getting source %s..." % source)

            with conn:
                with conn.cursor() as curs:
                    # FIXME: consider chunking into multiple transactions so as not to overload the trans. buffer
                    # we may also run harvesters in parallel, in which case we may prefer less contention between
                    # long-running transactions
                    for feature_association in source_feats:
                        try:
                            self._save_one(curs, feature_association)
                        except Exception as ex:
                            print(str(ex))
                            import ipdb
                            ipdb.set_trace()
                            
                conn.commit()

    def save(self, feature_association):
        """ write dict to a series of tables"""

        with self._connect() as conn:
            with conn.cursor() as curs:
                self._save_one(curs, feature_association)
            conn.commit()
