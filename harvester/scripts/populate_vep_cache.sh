#!/usr/bin/env bash

# ensure that the VEP cache exists at /opt/data; if not, download it
# and do a little maintenance (e.g., build the index, unzip the FASTA file since for the life of me i can't get it to work otherwise...)

# this script will do the following:
# 1) installs VEP cache, plugins, and fasta files
# 2) unzips fa.gz so we don't have to worry about VEP's weird gz/bgz issue
# 3) indexes it using samtools

VEP_CACHE="/opt/data"
VEP_VERSION="104"
ASSEMBLY="GRCh37"
FASTA_FILE="${VEP_CACHE}/homo_sapiens_merged/${VEP_VERSION}_${ASSEMBLY}/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"

# if the /opt/data folder already has the cache of the correct version in it, and ${FASTA_FILE}.fai already exists,
# just exit with success
if [[ -d "${VEP_CACHE}/homo_sapiens_merged/${VEP_VERSION}_${ASSEMBLY}" && -f "${FASTA_FILE}" && -f "${FASTA_FILE}.fai" ]]; then
  >&2 echo "Cache already exists, aborting"
  exit 1
fi

/opt/vep/src/ensembl-vep/INSTALL.pl --AUTO cfp --NO_UPDATE \
  --CACHE_VERSION ${VEP_VERSION} --CACHEDIR ${VEP_CACHE} \
  --PLUGINS all \
  --SPECIES homo_sapiens_merged --ASSEMBLY ${ASSEMBLY} && \
gzip -d -c ${FASTA_FILE}.gz > ${FASTA_FILE}  && \
samtools faidx ${FASTA_FILE}
