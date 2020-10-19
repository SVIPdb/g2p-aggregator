#!/usr/bin/env bash

# this script will do the following:
# 1. build a custom VEP image that correctly sets up the local user
# 2. ensure the VEP cache for homo_sapiens_merged_vep_99_GRCh37.tar.gz exists
# 3. accepts VEP-readable identifiers on stdin, emitting a JSON document with annotations as a result
# 4. inserts JSON representations into the database

docker build -t local_vep ./image

# point this where you want VEP to download its data
VEP_DATA_DIR="/data/resources/ensembldb"

docker run -it \
  --user=`id -u vep`:`id -g vep` \
  -e LOCAL_USER_ID="$( id -u $USER )" \
  -v "${VEP_DATA_DIR}:/opt/vep/.vep" \
  local_vep \
  "$@"
