#!/bin/bash

# go to /opt/sw and:
# git clone https://github.com/alireza202/TOBI.git

rf=/opt/ref
sw=/opt/sw

### batch_multi_array variable
bamdir=/scratch/Bam_files
main_outputdir=/Results/GBM
list_file=/scratch/Bam_files/running.txt
script=${sw}/TOBI/gbm_pipeline/run_pipeline.sh

### paths needed for annotation and filtering
Annotation_Filtering=${sw}/TOBI/gbm_pipeline/do_annotation_filtering.sh
BcfTools=/usr/local/bin/bcftools
java7=/usr/bin/java
SNPEFF_HOME=${sw}/snpEff-v3.6
SNPEFF="-jar $SNPEFF_HOME/snpEff.jar -c $SNPEFF_HOME/snpEff.config"
SNPSIFT="-jar $SNPEFF_HOME/SnpSift.jar"
vcfEffOnePerLine=${sw}/TOBI/gbm_pipeline/scripts/vcfEffOnePerLine.pl
PythonParsing=${sw}/TOBI/gbm_pipeline/scripts/parse_tsv.py
vcf2report=${sw}/TOBI/gbm_pipeline/scripts/vcf2report.py
filter_indel_techn_biol=${sw}/TOBI/gbm_pipeline/scripts/filter_indel_techn_biol.R

### references
ref=${rf}/homo_sapiens_assembly19/Homo_sapiens_assembly19.fasta
dbSnp138=${rf}/dbSnp138.vcf
clinvar=${rf}/clinvar_20140303.vcf
cosmic_variants=${rf}/CosmicVariants_v66_20130725.vcf 
super_normal=${rf}/219normals.cosmic.hitless100.noExactMut.mutless5000.all_samples.vcf
cbio=${rf}/cbio.fix.sort.vcf
dbNSFP=${rf}/dbNSFP2.4.txt.gz
dbNSFP_header=${rf}/header_fields.txt
