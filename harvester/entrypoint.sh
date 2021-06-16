#!/usr/bin/env bash

HARVESTERS="oncokb civic clinvar cosmic"

# declare profiles here, so we don't have to have 1000 entrypoint scripts
case "$1" in
	all)
	    HARVESTERS="cosmic clinvar civic"
            GENES=$( tail -n +2 "../data/genesets/cosmic/cancer_gene_census.csv" | cut -d ',' -f 1 | xargs )
            EXTRA_ARGS="--gene_chunk_size 15"
            ;;

        cosmic)
            HARVESTERS="cosmic"
            GENES=$( tail -n +2 "../data/genesets/cosmic/cancer_gene_census.csv" | cut -d ',' -f 1 | xargs )
            EXTRA_ARGS="--gene_chunk_size 15"
            ;;

        clinvar)
            HARVESTERS="clinvar"
            GENES=$( tail -n +2 "../data/genesets/cosmic/cancer_gene_census.csv" | cut -d ',' -f 1 | xargs )
            EXTRA_ARGS="--gene_chunk_size 15"
            ;;

        civic)
            HARVESTERS="civic"
            GENES=$( tail -n +2 "../data/genesets/cosmic/cancer_gene_census.csv" | cut -d ',' -f 1 | xargs )
            EXTRA_ARGS="--gene_chunk_size 15"
            ;;

        onkokb)
            HARVESTERS="oncokb"
            GENES=$( tail -n +2 "../data/genesets/cosmic/cancer_gene_census.csv" | cut -d ',' -f 1 | xargs )
            EXTRA_ARGS="--gene_chunk_size 15"
            ;;
        test)
            GENES="BCL9L BCLAF1 BCOR BCORL1 BCR BIRC3 BIRC6 BLM BMP5 BMPR1A BRAF BRCA1 BRCA2 BRD3 BRD4'"
            EXTRA_ARGS="--log=INFO"
            ;;
        simple)
            GENES="RAC1 ABL2"
            EXTRA_ARGS="--log=INFO"
            ;;
        merged)
            echo "Updating combined geneset from OncoKB and CiVIC"
            if ../data/genesets/merged/acquire_merged_list.sh; then
                HARVESTERS="oncokb civic clinvar cosmic"
                GENES=$( cat "../data/genesets/merged/merged_gene_list_latest.txt" | xargs )
            else
                echo "Failed to acquire merged geneset"
                exit 1
            fi
            ;;

        debug)
            GENES="RAC1"
            POSTGRES_DB="svip_api_debug"
            EXTRA_ARGS="--log=INFO"
            ;;

        *)
            echo $"Usage: $0 {cosmic|simple|debug}"
            exit 1
esac

# shift off the command that we parsed previously
shift 1

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

python harvester.py -ds \
  ${SILOS[$CHOSEN_SILO]} \
  --harvesters ${HARVESTERS} --genes ${GENES} \
  ${EXTRA_ARGS} $@
