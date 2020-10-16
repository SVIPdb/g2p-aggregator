#!/usr/bin/env bash

SOURCE_URL=$1
DEST_FILENAME=$2
COSMIC_CREDS=${3:-'.cosmic_creds'}

if [ -z ${SOURCE_URL} -o -z ${DEST_FILENAME} -o -z ${COSMIC_CREDS} ]; then
    echo "Usage: $0 <source URL> <dest filename> [credentials file]"
    exit 1
fi

if [ ! -f ${COSMIC_CREDS} ]; then
    echo "COSMIC cached credentials file ${COSMIC_CREDS} doesn't exist; please run ./get_cosmic_creds.sh first"
    exit 1
fi

# presumable .cosmic_creds exists now, so we can proceed
echo "Resolving temporary download URL using ${COSMIC_CREDS}..."

RESOLVED_URL=$( curl -s -S -H "Authorization: Basic $( cat ${COSMIC_CREDS} )" ${SOURCE_URL} | jq -r .url )

if [ -z "${RESOLVED_URL}" ]; then
    echo "Couldn't resolve download URL, please check your credentials file and that you have jq installed"
    exit 2
fi

echo "Downloading from ${RESOLVED_URL}..."
curl -L -o "${DEST_FILENAME}" "${RESOLVED_URL}"
