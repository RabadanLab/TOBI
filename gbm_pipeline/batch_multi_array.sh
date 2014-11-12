#!/bin/bash

# This is the batch script to run the case names given in "list_file"
# Bam files are in folders with case names
# Make sure bam files are in the correct folder
# main_outputdir will contain a folder for each case, in each folder two
# folders of "logs" and "output_folder". In each "output_folder",
# there will be a vccfiles_# folder for each segment of chromosome.

# To read all the variables needed
# Make sure to modify the config file accordingly

config=/Results/TOBI/gbm_pipeline/tobi_config_amazon.sh

if [[ -f ${config} ]]
then
        # on Amazon
        amazon=1
else
        # on HPC
        config=/ifs/home/c2b2/rr_lab/ar3177/bin/aws/gbm_pipeline/tobi_config_hpc.sh
fi

source ${config}

# flag=BAF
flag=B

if [[ ${flag} == *A* ]]
then
	mem=10G
else
	mem=4G
fi

java_mem=6
tim=1

# Number of cases
num_bam=$(cat ${list_file} | wc -l)

#for j in $(seq 1 1 ${num_bam}) 
for j in 1
do
	# Get the j th bam file
	# Adapt this part accordingly for your bam file and case names
	case_name=$(cat ${list_file} | sed -n ${j}p)
	short_case_name=$(echo ${case_name} | awk -F'-' '{print $2 $3}')
	BAM_file=${bamdir}/${case_name}/*/*.bam
	
	mkdir -p ${main_outputdir}/${case_name}
	mkdir -p ${main_outputdir}/${case_name}/logs
	mkdir -p ${main_outputdir}/${case_name}/output_folder
	
	# This qsubs all chromosome parts. Change the numbers after -t, if needed
	if [[ ${amazon} == 0 ]] 
	then
		qsub -t 30-30 -V -e ${main_outputdir}/${case_name}/logs -o ${main_outputdir}/${case_name}/logs \
			-N A-${short_case_name} -cwd -l mem=${mem},time=${tim}:: \
			${script} --bam ${BAM_file} --annot_filt ${Annotation_Filtering} --patient ${case_name} --ref ${ref} \
			--steps ${flag} --outputdir ${main_outputdir}/${case_name}/output_folder --memory ${java_mem} \
			--config_file ${config}
	else
		qsub -t 30-30 -V -e ${main_outputdir}/${case_name}/logs -o ${main_outputdir}/${case_name}/logs \
			-N A-${short_case_name} -cwd \
			${script} --bam ${BAM_file} --annot_filt ${Annotation_Filtering} --patient ${case_name} --ref ${ref} \
			--steps ${flag} --outputdir ${main_outputdir}/${case_name}/output_folder --memory ${java_mem} \
			--config_file ${config}
	fi
done