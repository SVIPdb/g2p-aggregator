#!/usr/bin/env python

# steps:
# 1. gather list of variants from SVIP db in a format that VEP can accept (e.g., HGVS.g)
# 2. run VEP, produce outputs
# 3. iterate through variants in the output, updating the db where neccessary

# ideally we should use sqlalchemy, because i'm tired of writing fragile psycopg2 code directly

# the command will eventually be similar to the following
# ./exec_vep.sh ./vep -id 'NC_000007.13:g.140453136A>T' --port 3337 --no_stats --database --json --output_file STDOUT
import os

from sqlalchemy import create_engine, MetaData, Table, select
from sqlalchemy.ext.declarative import declarative_base

engine = create_engine('postgresql://%(user)s:%(password)s@%(host)s/%(db)s' % {
    'user': os.environ.get('POSTGRES_USER'),
    'password': os.environ.get('POSTGRES_PASSWORD'),
    'host': os.environ.get('POSTGRES_HOST'),
    'db': os.environ.get('POSTGRES_DB'),
}, echo=False)

Base = declarative_base()
meta = MetaData(bind=engine)

# table reflection
variants = Table('api_variant', meta, autoload=True)

if __name__ == '__main__':
    # 1. produce set of variants with which to query VEP
    with engine.begin() as conn:
        results = conn.execute(select([variants.c.hgvs_g]))
        with open("input_variants.txt", "w") as fp:
            for row in results:
                fp.write("%s\n" % row['hgvs_g'])

    # 2. query VEP

    # 3. iterate over VEP response, update variant entries w/VEP data
