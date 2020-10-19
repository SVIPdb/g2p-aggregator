#!/usr/bin/env bash

set -o pipefail

# makes a folder named YYYY-MM-DD, then downloads gene lists from
# oncokb and civic into that directory.

# in addition to the original data from each source, two lists named
# (source)_gene_list.txt are produced which contain just the HUGO gene symbols
# in lexicographical order.

if [ ! -f ../../.oncokb_apikey ]; then
    ( cd ../.. ; ./get_oncokb_creds.sh )
fi

ONCOKB_KEY=$( cat ../../.oncokb_apikey )
CUR_DIR="cached_$( date +%F )"

mkdir -p ${CUR_DIR}
pushd ${CUR_DIR}

# download and process CIViC genelist
curl --show-error -O https://civicdb.org/downloads/nightly/nightly-GeneSummaries.tsv
cut -d$'\t' -f 3 nightly-GeneSummaries.tsv | tail +2 | sort > civic_gene_list.txt

# download and process OncoKB genelist (only includes curated genes)
curl --show-error -H "Authorization: Bearer ${ONCOKB_KEY}" https://www.oncokb.org/api/v1/utils/cancerGeneList | jq > oncokb_genes.json
jq -r '.[].hugoSymbol' oncokb_genes.json | sort > oncokb_gene_list.txt || (
    echo "Error occurred when parsing OncoKB data:"
    echo $( cat oncokb_genes.json | jq -r '{.status, .title, .detail}')
)

# merge the two lists, producing a de-duplicated result
( cat civic_gene_list.txt ; cat oncokb_gene_list.txt ) | sort -u > merged_gene_list.txt

# output some stats
echo "CIViC Genes: $( wc -l civic_gene_list.txt )"
echo "OncoKB Genes: $( wc -l oncokb_gene_list.txt )"
echo "Merged Genes: $( wc -l merged_gene_list.txt )"

popd

# symlink merged_gene_list_latest.txt to the new list so that we can refer to it in other places, e.g. the harvester's entrypoint
rm merged_gene_list_latest.txt 2>/dev/null
ln -s $( realpath "${CUR_DIR}/merged_gene_list.txt" ) merged_gene_list_latest.txt
