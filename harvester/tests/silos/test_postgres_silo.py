import argparse

from psycopg2 import sql

from silos import postgres_silo
from silos.postgres_silo import PostgresSilo


def test_postgres_array_conversion():
    silo = _make_silo()
    from silos.postgres_silo import array_to_sql

    with silo._connect() as conn:
        # test simple int arrays
        result = array_to_sql([1, 2, 3]).as_string(conn)
        assert result == 'ARRAY[1, 2, 3]'

        result = array_to_sql(['civic', 'oncokb', 'clinvar']).as_string(conn)
        assert result == """ARRAY['civic', 'oncokb', 'clinvar']"""


def test_postgres_array_operations():
    silo = _make_silo()
    from silos.postgres_silo import array_to_sql

    with silo._connect() as conn:
        with conn.cursor() as curs:
            items = ['civic', 'oncokb', 'clinvar']

            # attempt to perform a subset operation, which should be true in this case
            curs.execute(
                sql.SQL("select {}::text[] <@ {}").format(
                    array_to_sql(items),
                    array_to_sql(items + ['other', 'entries'])
                )
            )
            result = curs.fetchone()[0]
            assert result is True

            # another subset operation, this time false
            curs.execute(
                sql.SQL("select {}::text[] <@ {}").format(
                    array_to_sql(items),
                    array_to_sql(['civic', 'blah'])
                )
            )
            result = curs.fetchone()[0]
            assert result is False


def test_postgres_array_merging():
    silo = _make_silo()
    from silos.postgres_silo import array_to_sql

    with silo._connect() as conn:
        with conn.cursor() as curs:
            # attempt union, then distincts
            curs.execute(
                sql.SQL("select array_distinct({} || {})").format(
                    array_to_sql('[1, 2, 2]'),
                    array_to_sql('[2, 3, 4]')
                )
            )
            result = curs.fetchone()[0]
            assert result == [1, 2, 3, 4]


def _make_silo():
    argparser = argparse.ArgumentParser()
    postgres_silo.populate_args(argparser)
    args = argparser.parse_args(['--pg_host', 'db', '--pg_db', 'svip_api'])
    silo = PostgresSilo(args)
    return silo
