#!/bin/bash

bamdir=/ifs/scratch/c2b2/rr_lab/ar3177/TCGA/BAM_files
outputdir=/ifs/scratch/c2b2/rr_lab/ar3177/Results/GBM/Done
finaldir=/ifs/scratch/c2b2/rr_lab/ar3177/Results/GBM/final_tables

stepstr=12
list_file=running.txt
#list_file=with_src.txt

prefix=raw
suffix=vcf.all.annotations_filt_indel_techn_biol.tsv

# concatenating tsv files of each case
if [[ $stepstr == *1* ]]; then
	# number of new cases
	num_bam=$(cat ${bamdir}/${list_file} | wc -l)
	
	
	for j in $(seq 1 1 ${num_bam})
#	for j in 1
	do
		# get the j th bam file
		case_name=$(cat ${bamdir}/${list_file} | sed -n ${j}p)
		echo ${case_name}
		
		# delete the final file if it exists (to avoid appending it)
		rm -f ${outputdir}/${case_name}/output_folder/all_mutations_${case_name}.txt;
		
		# putting the header
		cat ${outputdir}/${case_name}/output_folder/vcffiles_2/${prefix}_2.${suffix} | head -1 > \
		    ${outputdir}/${case_name}/output_folder/all_mutations_${case_name}.txt
		
		list_files=""
		for i in {1..73}
		do
			
			# check if the tsv file is empty
			num_lines=$(cat ${outputdir}/${case_name}/output_folder/vcffiles_${i}/${prefix}_${i}.${suffix} | wc -l)

			if [[ $num_lines > 1 ]]
			then
				# getting all the contents without the headers
				echo -ne "${i} "
				cat ${outputdir}/${case_name}/output_folder/vcffiles_${i}/${prefix}_${i}.${suffix} \
					| sed '1d' \
					>> ${outputdir}/${case_name}/output_folder/all_mutations_${case_name}.txt
			fi
		done
		
		echo -e "\n"
	done
fi

# concatenating all cases
if [[ $stepstr == *2* ]]
then
	today=$(date | awk '{print $2 "-" $3 "-" $6}')
	
	output_file=${finaldir}/all_GBM_mutations_${today}_filt_indel_techn_biol.txt
	
	# number of all cases
	num_bam=$(cat ${bamdir}/${list_file} | wc -l)
	
	# delete the final file if it exists (to avoid appending it)
	rm -f ${output_file};
	
	for i in $(seq 1 1 ${num_bam})
	do
		# get the i th bam file
		case_name=$(cat ${bamdir}/${list_file} | sed -n ${i}p)
		
		# putting the header
		if [ ${i} -eq 1 ]
		then
			head -1 ${outputdir}/${case_name}/output_folder/all_mutations_${case_name}.txt > ${output_file}
		fi
		
		# getting all the contents without the headers
		cat ${outputdir}/${case_name}/output_folder/all_mutations_${case_name}.txt | sed '1d' \
			>> ${output_file}
	done
	
fi


