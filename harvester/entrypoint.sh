#!/usr/bin/env bash

declare -A SILOS
SILOS[postgres]="--silos postgres --pg_host ${POSTGRES_HOST} --pg_db ${POSTGRES_DB} --pg_user ${POSTGRES_USER} --pg_pass ${POSTGRES_PASSWORD}"
SILOS[elastic]="--silos elastic --elastic_search elasticsearch"
SILOS[file]="--silos file --file_output_dir=/output/raw"

CHOSEN_SILO="postgres"

if [[ ${CHOSEN_SILO} == "postgres" ]]; then
    # wait for postgres to wake up
    until PGPASSWORD="$POSTGRES_PASSWORD" psql -d  "$POSTGRES_DB" -h "$POSTGRES_HOST" -U "$POSTGRES_USER" -c '\q' 2> /dev/null; do
        >&2 echo "Postgres is unavailable at $POSTGRES_HOST - sleeping"
        sleep 1
    done
fi

# and do a test harvest, dumping the results to a file
python harvester.py -d \
  ${SILOS[$CHOSEN_SILO]} \
  --harvesters oncokb civic --genes BRAF EGFR
  # --harvesters oncokb civic jax --genes BRAF --phases all
