#!/bin/bash

# This is the batch script to run the samples given in "list_file"
# Make sure bam files are in the correct folder

bamdir=/ifs/scratch/c2b2/rr_lab/ar3177/TCGA/BAM_files
outputdir=/ifs/scratch/c2b2/rr_lab/ar3177/Results/running

# flag=BAF
flag=B

mem=8G
tim=1

list_file=running.txt

script=/ifs/home/c2b2/rr_lab/ar3177/bin/TOBI/gbm_pipeline/run_pipeline.sh

# Number of cases
num_bam=$(cat ${bamdir}/${list_file} | wc -l)

for j in $(seq 1 1 ${num_bam}) 
#for j in 1
do
	# Get the j th bam file
	case_name=$(cat ${bamdir}/${list_file} | sed -n ${j}p)
	short_case_name=$(echo ${case_name} | awk -F'-' '{print $2 $3}')
	BAM_file=${bamdir}/${case_name}/*/*.bam
	
	mkdir -p ${outputdir}/${case_name}
	mkdir -p ${outputdir}/${case_name}/logs
	mkdir -p ${outputdir}/${case_name}/output_folder

	for i in {1..73}
#	for i in 1
	do
		qsub -V -e ${outputdir}/${case_name}/logs -o ${outputdir}/${case_name}/logs -N A${i}-${short_case_name} -cwd -l mem=${mem},time=${tim}:: \
			${script} --bam $BAM_file --index $i --patient ${case_name} --steps $flag --outputdir ${outputdir}/${case_name}/output_folder --memory 6
	done
done