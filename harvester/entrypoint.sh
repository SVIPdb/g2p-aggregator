#!/usr/bin/env bash

# echo "POSTGRES_PASSWORD: $POSTGRES_PASSWORD"
# echo "POSTGRES_DB: $POSTGRES_DB"
# echo "POSTGRES_HOST: $POSTGRES_HOST"
# echo "POSTGRES_USER: $POSTGRES_USER"

# wait for postgres to wake up
until PGPASSWORD="$POSTGRES_PASSWORD" psql -d  "$POSTGRES_DB" -h "$POSTGRES_HOST" -U "$POSTGRES_USER" -c '\q' 2> /dev/null; do
    >&2 echo "Postgres is unavailable at $POSTGRES_HOST - sleeping"
    sleep 1
done

# and do a test harvest, dumping the results to a file
python harvester.py -ds --silos file --file_output_dir=/app/output_new \
  --harvesters oncokb civic jax --genes BRAF --phases all
