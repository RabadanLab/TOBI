#!/bin/bash

bamdir=/ifs/scratch/c2b2/rr_lab/ar3177/TCGA/BAM_files
outputdir=/ifs/scratch/c2b2/rr_lab/ar3177/Results/running

list_files=running.txt
#list_files=all_normal.txt

# Number of cases
num_bam=$(cat ${bamdir}/${list_files} | wc | awk '{print $1}')

for j in $(seq 1 1 ${num_bam}) 
#for j in 1
do
	# Get the j th bam file
	case_name=$(cat ${bamdir}/${list_files} | sed -n ${j}p)

	for i in {1..73}
	do
#		mv ${outputdir}/${case_name}/output_folder/vcffiles_${i}/final_raw_${i}.vcf.all.annotations.tsv ${outputdir}/${case_name}/output_folder/vcffiles_${i}/old_final_raw_${i}.vcf.all.annotations.tsv
		filename=${outputdir}/${case_name}/output_folder/vcffiles_${i}/raw_${i}.vcf.all.annotations_filt_indel_techn_biol.tsv
		if [ ! -f ${filename} ];
		then
			echo ${filename}
		fi
	done
done