#!/bin/bash

dir=/ifs/home/c2b2/rr_lab/ar3177/scratch/Results/GBM
outputdir=/ifs/home/c2b2/rr_lab/ar3177/scratch/Results

not_present=/ifs/home/c2b2/rr_lab/ar3177/scratch/Results/Finding_lost_mutations/missed_TCGA_06_5410.txt

# not_present.txt is a list of TCGA mutations not present in the final mutations table from our pipeline. It consists of "gene" "case_name" "chr" "pos".


while read gene case_name chr pos
do
	# Find the corresponsing region
	for i in {1..73}
	do
		set=$(cat ${dir}/partitions.txt | sed -n ${i}p)
		CHR=$(echo ${set} | cut -f1 -d' ')
		POS1=$(echo ${set} | cut -f2 -d' ')
		POS2=$(echo ${set} | cut -f3 -d' ')
			
		if [ "$chr" = "$CHR" ] && [ $pos -ge $POS1 ] && [ $pos -le $POS2 ]
		then
			num=${i}
		fi
	done
	
	if [ -z "$num" ]
	then
		echo "No match was found for" $case_name $chr $pos
	fi
	
	# Finding the corresponding QUAL and DP, etc
	line=$(cat ${outputdir}/${case_name}/output_folder/vcffiles_${num}/raw_${num}.vcf.all.annotations.rep.vcf | grep ${pos})
	
#	ID=$(echo ${line} | cut -f3  -d' ')
#	echo -e ${gene}'\t'${case_name}'\t'${num}'\t'${chr}'\t'${pos}'\t'${ID}
	
	if [ -z "${line}" ]
	then
	    echo ${chr} ${pos} "cannot be found in" ${outputdir}/${case_name}/output_folder/vcffiles_${num}/raw_${num}.vcf.all.annotations.rep.vcf
	else
		QUAL=$(echo ${line} | cut -f6  -d' ')
		DP=$(echo ${line} | awk -F'DP=' '{print $2}' | awk -F';' '{print $1}')
		DP4=$(echo ${line} | awk -F'DP4=' '{print $2}' | awk -F';' '{print $1}')
		MQ=$(echo ${line} | awk -F'MQ=' '{print $2}' | awk -F';' '{print $1}')
		PV4=$(echo ${line} | awk -F'PV4=' '{print $2}' | awk -F';' '{print $1}')
		EFF=$(echo ${line} | awk -F'EFF=' '{print $2}' | awk -F'(' '{print $1}')
		COMMON=$(echo ${line} | awk -F'COMMON=' '{print $2}' | awk -F';' '{print $1}')
		G5=$(echo ${line} | awk -F'G5=' '{print $2}' | awk -F';' '{print $1}')
		G5A=$(echo ${line} | awk -F'G5A=' '{print $2}' | awk -F';' '{print $1}')
		MEGA=$(echo ${line} | awk -F'MEGANORMAL_ID=' '{print $2}' | awk -F';' '{print $1}')
		
		type=SNP
		echo -e ${gene}'\t'${case_name}'\t'${num}'\t'${chr}'\t'${pos}'\t'${type}'\t'${QUAL}'\t'${DP}'\t'${DP4}'\t'${MQ}'\t'${PV4}'\t'${EFF}'\t'${COMMON}'\t'${G5}'\t'${G5A}'\t'${MEGA}
	fi
done <  <( cat ${not_present} )
