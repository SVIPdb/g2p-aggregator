#!/usr/bin/env bash

function die() {
    echo "$1"
    exit ${2:-1}
}

VARIATION_XML_PATH="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/clinvar_variation/"
VARIATION_XML_FILE="ClinVarVariationRelease_00-latest.xml.gz"
VARIATION_MD5_FILE="${VARIATION_XML_FILE}.md5"

# grab the latest MD5 from ClinVarVariation
latest_md5=$( curl -s -S -f ${VARIATION_XML_PATH}${VARIATION_XML_FILE}.md5 | cut -d ' ' -f 1 ) || die "Couldn't retrieve MD5 from server"

# check if we have an existing MD5 file and if its contents match ours
if [ -f "${VARIATION_MD5_FILE}" ] && ( diff "${VARIATION_MD5_FILE}" <( echo $latest_md5 ) 2>&1 >/dev/null  ); then
    die "local MD5 file exists with matching checksum to server's version, aborting"
fi

# download the full release, write the checksum out
curl -f -L -O "${VARIATION_XML_PATH}${VARIATION_XML_FILE}" || die "Couldn't download ${VARIATION_XML_FILE}"
echo $latest_md5 > ${VARIATION_MD5_FILE}

# also delete the clinvar harvester's variant count cache file since it's invalid now
rm .var_counts_cache.json 2>/dev/null
