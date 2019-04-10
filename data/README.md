# data

## overview

This directory contains the 'raw' files from various sources necessary for g2p harvesting.

## contents

* Makefile: produces summaries, indices, etc. for various files, and fetches a few.
  - *FIXME: add better documentation for the Makefile's targets*

* [pmkb_interpretations.json](https://s3-us-west-2.amazonaws.com/g2p-0.7/unprocessed-files/pmkb_interpretations.json)
  * manual download
  * raw interpretations from cornell linda at standardmolecular dot com

* [molecularmatch_trials.json](https://s3-us-west-2.amazonaws.com/g2p-0.7/unprocessed-files/molecularmatch_trials.json)
  * pre-harmonization download of molecularmatch trials
  * recreate  via `python harvester.py --silos file --harvesters molecularmatch_trials  --phases harvest`

* cgi
  * [cgi_biomarkers_per_variant.tsv](https://www.cancergenomeinterpreter.org/biomarkers)
  * [catalog_of_validated_oncogenic_mutations.tsv](https://www.cancergenomeinterpreter.org/mutations)
  * manual download
  * raw evidence

* [allAnnotatedVariants.txt](http://oncokb.org/api/v1/utils/allAnnotatedVariants.txt), [allActionableVariants.txt](http://oncokb.org/api/v1/utils/allActionableVariants.txt)
  * download by clicking link, or going to [oncokb.org/#/dataAccess](http://oncokb.org/#/dataAccess)
  * used to run harvester on oncokb
  * 'oncokb\_' prefix added to files for clarity

* [data_mutations_extended_1.0.1.txt,data_clinical_1.0.1.txt](https://www.synapse.org/#!Synapse:syn7851250 )
  * manual registration and download
  * cohort for GENIE analysis notebook

* [harvester/cosmic_lookup_table.tsv](https://grch37-cancer.sanger.ac.uk/cosmic/files?data=/files/grch37/cosmic/v81/CosmicMutantExport.tsv.gz)
  * manual download
  * pre processing done by harvester/Makefile including [harvester/oncokb_all_actionable_variants.tsv](http://oncokb.org/api/v1/utils/allActionableVariants.txt)

* [msk_impact_2017/*](http://www.cbioportal.org/study?id=msk_impact_2017#summary)
  * manual download  from 'Download Data' button at top of page
  * cohort for GENIE analysis notebook

* non_alt_loci_set.json
  * wget ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/json/non_alt_loci_set.json
  * gene symbol validation

* COSMIC
  * CosmicMutantExport.tsv.gz: a set of samples from COSMIC. must be manually downloaded, used to produce other reference files. the Makefile has this to say about it:
      > COSMIC Mutation Data; can download from https://grch37-cancer.sanger.ac.uk/cosmic/download
      > Note: we want GrCh37 build for now since most databases use this build.
      > https://grch37-cancer.sanger.ac.uk/cosmic/files?data=/files/grch37/cosmic/v81/CosmicMutantExport.tsv.gz
  * CosmicMutantExport.tsv: extracted version of the above, produced during make all.
  * CosmicMutantExport_sorted.tsv: COSMIC samples sorted by gene, then by alteration.

* disease_alias.tsv: used by harvester/normalizers/disease_normalizer.py to map aliases of diseases to a single disease name.
* ensembl_refseq_slim.tsv: provides (incorrect?) mappings from ensembl transcript accession IDs to NCBI refseq accession IDs
