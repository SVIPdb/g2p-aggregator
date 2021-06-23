#!/usr/bin/env bash

TMPFILE="/tmp/intermediate.txt"

VEP_CACHE="/opt/data"
VEP_VERSION="104"
ASSEMBLY="GRCh37"

# created via the following command:
# gzip --decompress --to-stdout Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz | \
#  bgzip --index --index-name Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.bgz.tbi -@4 -t \
#  > Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.bgz
FASTA_BGZ="Homo_sapiens.${ASSEMBLY}.75.dna.primary_assembly.fa.bgz"

vep --format vcf \
  --offline --merged --use_given_ref --canonical --force_overwrite --hgvs --hgvsg --symbol --af_gnomad --uniprot \
  --max_af --sift b --variant_class \
  --fasta ${VEP_CACHE}/homo_sapiens_merged/${VEP_VERSION}_${ASSEMBLY}/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
  --dir_cache ${VEP_CACHE} -o ${TMPFILE} && \
  filter_vep -i ${TMPFILE} --filter "Consequence is not upstream_gene_variant and Consequence is not downstream_gene_variant"
