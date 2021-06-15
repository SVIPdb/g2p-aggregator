# we rely on some bash features (e.g., <( cmd )), so we use bash as the makefile's shell
SHELL=bash

# COSMIC Mutation Data; can download from https://grch37-cancer.sanger.ac.uk/cosmic/download
# Note: we want GrCh37 build for now since most databases use this build.
# https://grch37-cancer.sanger.ac.uk/cosmic/files?data=/files/grch37/cosmic/v81/CosmicMutantExport.tsv.gz
COSMIC_LOOKUPGEN_SCRIPT=../harvester/lookups/cosmic_lookup_table.py

BENCHMARK_DATASET=oncokb

# previous version: 20170208
CGI_VERSION=20180117

all: oncokb_all_actionable_variants.tsv cgi_biomarkers_per_variant.tsv \
		cosmic_lookup_table.tsv CosmicMutantExport.sqlite \
		clinvar/ClinVarVariationRelease_00-latest.xml \
		cgi_mut_benchmarking.tsv oncokb_mut_benchmarking.tsv \
		benchmark_results.txt \
		genesets/cosmic/cancer_gene_census.csv \
		pmkb_interpretations.json

## Do bulk download of OncoKB and set up TSV for ingestion.
oncokb_all_actionable_variants.tsv:
	curl -o oncokb_all_actionable_variants_$(shell date +"%Y%m%d").tsv http://oncokb.org/api/v1/utils/allActionableVariants.txt
	ln -s oncokb_all_actionable_variants_$(shell date +"%Y%m%d").tsv oncokb_all_actionable_variants.tsv

## Do bulk download of CGI and set up TSV for ingestion.
cgi_biomarkers_per_variant.tsv:
	mkdir -p cgi_biomarkers_$(CGI_VERSION)
	cd cgi_biomarkers_$(CGI_VERSION) && \
	curl -Ok https://www.cancergenomeinterpreter.org/data/cgi_biomarkers_$(CGI_VERSION).zip && \
	unzip cgi_biomarkers_$(CGI_VERSION).zip
	cp cgi_biomarkers_$(CGI_VERSION)/cgi_biomarkers_per_variant.tsv .

catalog_of_validated_oncogenic_mutations.tsv:
	curl -L -o catalog_of_validated_oncogenic_mutations_latest.zip 'https://www.cancergenomeinterpreter.org/data/catalog_of_validated_oncogenic_mutations_latest.zip?ts=20180216'
	unzip catalog_of_validated_oncogenic_mutations_latest.zip catalog_of_validated_oncogenic_mutations.tsv

## COSMIC requires credentials to download its datafiles; this script will cache those credentials
## this target runs a script that notifies the user that they must first register for COSMIC, then
## instructs them on how to create .cosmic_creds, which allows for automated COSMIC downloads
.cosmic_creds:
	./get_cosmic_creds.sh

CosmicMutantExport.tsv.gz: .cosmic_creds
	./get_cosmic_datafile.sh 'https://cancer.sanger.ac.uk/cosmic/file_download/GRCh37/cosmic/v92/CosmicMutantExport.tsv.gz' CosmicMutantExport.tsv.gz

## Create COSMIC lookup table suitable for annotating variants with chrom/start/end/ref/alt
cosmic_lookup_table.tsv: CosmicMutantExport.tsv.gz
	gunzip -c CosmicMutantExport.tsv.gz | tail -n +2 | python2 ${COSMIC_LOOKUPGEN_SCRIPT} create_table > temp.tsv
	(head -n 1 temp.tsv ; tail -n +2 temp.tsv | sort | uniq) > cosmic_lookup_table.tsv
	rm -f temp.tsv

## produce a version of COSMIC indexed by gene, making it much faster to process entries later on
CosmicMutantExport.sqlite: CosmicMutantExport.tsv.gz
	./cosmic_to_sqlite.sh <( gunzip -c CosmicMutantExport.tsv.gz | tail -n +2 ) CosmicMutantExport.sqlite

## ClinVar XML download
clinvar/ClinVarVariationRelease_00-latest.xml:
	cd ./clinvar && ./get_clinvar_latest.sh && gunzip ClinVarVariationRelease_00-latest.xml.gz

## Create lookup tables for proteins
proteins:
	mkdir -p ./proteins/
	wget -r -l1 -np -nd "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/" -P ./proteins -A "human.*.protein.faa.gz"

## Create cgi_mut_benchmarking.tsv file.
cgi_mut_benchmarking.tsv:
	curl -s https://www.cancergenomeinterpreter.org/data/cgi_biomarkers_latest.zip -o cgi_biomarkers_latest.zip && \
	unzip cgi_biomarkers_latest.zip && \
	grep MUT cgi_biomarkers.tsv | \
	egrep -v "oncogenic mutation|insertion|deletion|wildtype|promoters|_variant" | \
	cut -f4 | sed 's/\:/	/g' > cgi_mut_benchmarking.tsv

## Create oncokb_mut_benchmarking.tsv file.
oncokb_mut_benchmarking.tsv: oncokb_all_actionable_variants.tsv
	cut -f1-2 oncokb_all_actionable_variants.tsv | \
	egrep -v "del|ins|Fusion|Oncogenic|mutations|dup|Wildtype|Amplification|Truncating|splice" \
	> oncokb_mut_benchmarking.tsv

## Benchmark file.
benchmark_results.txt: cosmic_lookup_table.tsv $(BENCHMARK_DATASET)_mut_benchmarking.tsv
	 python2 ${COSMIC_LOOKUPGEN_SCRIPT} benchmark --lookup-file cosmic_lookup_table.tsv \
	 	--benchmark-file $(BENCHMARK_DATASET)_mut_benchmarking.tsv > benchmark_results.txt

## OncoKB requires credentials to access the API; this prompts the user to register,
## then prompts the user for their API key, which is stored in .oncokb_apikey
.oncokb_apikey:
	./get_oncokb_creds.sh

## COSMIC census geneset, used as our canonical cancer gene list during harvesting
genesets/cosmic/cancer_gene_census.csv: .cosmic_creds
	./get_cosmic_datafile.sh 'https://cancer.sanger.ac.uk/cosmic/file_download/GRCh37/cosmic/v92/cancer_gene_census.csv' genesets/cosmic/cancer_gene_census.csv

## Merged geneset, used when we want to retrieve all the genes in OncoKB + CIViC
genesets/merged/merged_gene_list_latest.txt: .oncokb_apikey
	cd genesets/merged ; ./acquire_merged_list.sh

pmkb_interpretations.json:
	curl -L -O 'https://s3-us-west-2.amazonaws.com/g2p-0.7/unprocessed-files/pmkb_interpretations.json'

clean:
	rm -rf \
	  .oncokb_apikey .cosmic_creds \
	  cosmic_lookup_table.tsv CosmicMutantExport.tsv.gz CosmicMutantExport_sorted.tsv \
	  cgi_mut_benchmarking.tsv oncokb_mut_benchmarking.tsv benchmark.txt \
	  cgi_biomarkers_$(CGI_VERSION) cgi_biomarkers_per_variant.tsv \
	  oncokb_all_actionable_variants*.tsv \
	  pmkb_interpretations.json \
	  genesets/cosmic/cancer_gene_census.csv \
	  genesets/merged/merged_gene_list_latest.txt
