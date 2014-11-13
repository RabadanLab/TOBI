#!/bin/bash

### batch_multi_array variable
bamdir=/ifs/scratch/c2b2/rr_lab/ar3177/TCGA/BAM_files
main_outputdir=/ifs/scratch/c2b2/rr_lab/ar3177/Results/GBM
list_file=/ifs/scratch/c2b2/rr_lab/ar3177/TCGA/BAM_files/tmp_running.txt
script=/ifs/home/c2b2/rr_lab/ar3177/bin/aws/varCall_filtering/run_pipeline.sh

### paths needed for annotation and filtering
Annotation_Filtering=/ifs/home/c2b2/rr_lab/ar3177/bin/aws/varCall_filtering/do_annotation_filtering.sh
BcfTools=/ifs/home/c2b2/rr_lab/shares/bin/samtools-0.1.19/bcftools/bcftools
java7=/nfs/apps/java/1.7.0_25/bin/java
SNPEFF_HOME=/ifs/scratch/c2b2/rr_lab/shares/snpEff-v3.6
SNPEFF="-jar $SNPEFF_HOME/snpEff.jar -c $SNPEFF_HOME/snpEff.config"
SNPSIFT="-jar $SNPEFF_HOME/SnpSift.jar"
vcfEffOnePerLine=/ifs/home/c2b2/rr_lab/ar3177/bin/TOBI/varCall_filtering/scripts/vcfEffOnePerLine.pl
PythonParsing=/ifs/home/c2b2/rr_lab/ar3177/bin/TOBI/varCall_filtering/scripts/parse_tsv.py
vcf2report=/ifs/home/c2b2/rr_lab/ar3177/bin/TOBI/varCall_filtering/scripts/vcf2report.py
filter_indel_techn_biol=/ifs/home/c2b2/rr_lab/ar3177/bin/TOBI/varCall_filtering/scripts/filter_indel_techn_biol.R

### references
ref=/ifs/home/c2b2/rr_lab/ar3177/shares_scratch/ref/hg19/downloads/homo_sapiens_assembly19/Homo_sapiens_assembly19.fasta
dbSnp138=/ifs/scratch/c2b2/rr_lab/shares/snpEff-v2.1/dbSnp138.vcf
clinvar=/ifs/scratch/c2b2/rr_lab/shares/snpEff-v3.5/clinvar_20140303.vcf
cosmic_variants=/ifs/scratch/c2b2/rr_lab/shares/ref/COSMIC/CosmicVariants_v66_20130725.vcf 
super_normal=/ifs/scratch/c2b2/rr_lab/jw2983/mutation_only_tumor/meganormal/219normals.cosmic.hitless100.noExactMut.mutless5000.all_samples.vcf
cbio=/ifs/scratch/c2b2/rr_lab/shares/ref/vcfs/cbio.fix.sort.vcf
dbNSFP=/ifs/scratch/c2b2/rr_lab/shares/snpEff-v3.6/data/dbNSFP2.4.txt.gz
dbNSFP_header=/ifs/home/c2b2/rr_lab/ar3177/bin/TOBI/varCall_filtering/scripts/header_fields.txt
