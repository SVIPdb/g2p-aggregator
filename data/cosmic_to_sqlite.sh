#!/usr/bin/env bash

SOURCE_TSV=${1:-"CosmicMutantExport_sorted.tsv"}
TARGET_DBFILE=${2:-"CosmicMutantExport.sqlite"}

echo "* Reading from ${SOURCE_TSV} and writing to db file ${TARGET_DBFILE}..."

SCHEMA='CREATE TABLE variants (
    "Gene name" VARCHAR,
    "Accession Number" VARCHAR,
    "Gene CDS length" FLOAT,
    "HGNC ID" FLOAT,
    "Sample name" VARCHAR,
    "ID_sample" FLOAT,
    "ID_tumour" FLOAT,
    "Primary site" VARCHAR,
    "Site subtype 1" VARCHAR,
    "Site subtype 2" VARCHAR,
    "Site subtype 3" VARCHAR,
    "Primary histology" VARCHAR,
    "Histology subtype 1" VARCHAR,
    "Histology subtype 2" VARCHAR,
    "Histology subtype 3" VARCHAR,
    "Genome-wide screen" BOOLEAN,
    "GENOMIC_MUTATION_ID" VARCHAR,
    "LEGACY_MUTATION_ID" VARCHAR,
    "MUTATION_ID" FLOAT,
    "Mutation CDS" VARCHAR,
    "Mutation AA" VARCHAR,
    "Mutation Description" VARCHAR,
    "Mutation zygosity" BOOLEAN,
    "LOH" BOOLEAN,
    "GRCh" FLOAT,
    "Mutation genome position" VARCHAR,
    "Mutation strand" VARCHAR,
    "SNP" BOOLEAN,
    "Resistance Mutation" VARCHAR,
    "FATHMM prediction" VARCHAR,
    "FATHMM score" FLOAT,
    "Mutation somatic status" VARCHAR,
    "Pubmed_PMID" FLOAT,
    "ID_STUDY" FLOAT,
    "Sample Type" VARCHAR,
    "Tumour origin" VARCHAR,
    "Age" FLOAT,
    "HGVSP" VARCHAR,
    "HGVSC" VARCHAR,
    "HGVSG" VARCHAR
);'

# create the initial schema
sqlite3 "${TARGET_DBFILE}" "${SCHEMA}"

# create a local tempdir, since we need a lot of disk space and we're probably in a mount with a lot...
mkdir -p tmp
TMP_DIR=$( realpath ./tmp )

# load up our table and create indices on it (this can take time...)
sqlite3 -csv -separator $'\t' -batch "${TARGET_DBFILE}" <<EOF
.mode tabs
PRAGMA temp_store_directory = '${TMP_DIR}';
.import "${SOURCE_TSV}" variants
create index gene_idx on variants ("Gene name");
EOF
