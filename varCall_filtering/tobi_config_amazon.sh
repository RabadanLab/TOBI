#!/bin/bash

# go to a volume (e.g. /Results) and:
# git clone https://github.com/alireza202/TOBI.git

rf=/opt/ref
sw=/opt/sw
tobi=/Results/TOBI

### batch_multi_array variable
bamdir=/scratch/Bam_files
main_outputdir=/Results/GBM
list_file=/scratch/Bam_files/running.txt
script=${tobi}/varCall_filtering/run_pipeline.sh

### paths needed for annotation and filtering
Annotation_Filtering=${tobi}/varCall_filtering/do_annotation_filtering.sh
BcfTools=/usr/local/bin/bcftools
java7=/usr/bin/java
SNPEFF_HOME=${sw}/snpEff-v3.6
SNPEFF="-jar $SNPEFF_HOME/snpEff.jar -c $SNPEFF_HOME/snpEff.config"
SNPSIFT="-jar $SNPEFF_HOME/SnpSift.jar"
vcfEffOnePerLine=${tobi}/varCall_filtering/scripts/vcfEffOnePerLine.pl
PythonParsing=${tobi}/varCall_filtering/scripts/parse_tsv.py
vcf2report=${tobi}/varCall_filtering/scripts/vcf2report.py
filter_indel_techn_biol=${tobi}/varCall_filtering/scripts/filter_indel_techn_biol.R

### references
ref=${rf}/homo_sapiens_assembly19/Homo_sapiens_assembly19.fasta
dbSnp138=${rf}/dbSnp138.vcf
clinvar=${rf}/clinvar_20140303.vcf
cosmic_variants=${rf}/CosmicVariants_v66_20130725.vcf 
super_normal=${rf}/219normals.cosmic.hitless100.noExactMut.mutless5000.all_samples.vcf
cbio=${rf}/cbio.fix.sort.vcf
dbNSFP=${rf}/dbNSFP2.4.txt.gz
dbNSFP_header=${rf}/header_fields.txt
