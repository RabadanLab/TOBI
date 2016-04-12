#!/bin/bash

# creates a "true somatic call" file from a somatic MAF from TCGA

input_file=${1}
output_file=${2}

# selecting Chromosome  Start_position  Reference_Allele  Tumor_Seq_Allele1  Tumor_Seq_Allele2  Tumor_Sample_Barcode
# 1st awk removes insertions and deletions 
# 2nd awk selects TCGA barcode: TCGA-Tissue source site-Participant-tissue type
sed '/^#/d' ${input_file} | awk '(($10!="DEL") && ($10!="INS"))' | cut -f1,5,6,11,13,16 | sed '1d' > tmp_maf.txt #| 
	#awk -F $'\t' 'BEGIN {OFS = FS}{split($6,a,"-"); $6=a[1]"-"a[2]"-"a[3]"-"a[4];}1' > tmp_maf.txt 
head -n1 tmp_maf.txt

echo "gene_symbol"$'\t'"chr"$'\t'"start_position"$'\t'"reference_allele"$'\t'"variant_allele"$'\t'"case_id" > ${output_file}
cat tmp_maf.txt >> ${output_file}; head -n2 ${output_file}

#also generate tumor barcode list file
cut -f6 ${output_file} | sort -u | sed '1d' > patient.${output_file}
rm tmp_maf.txt