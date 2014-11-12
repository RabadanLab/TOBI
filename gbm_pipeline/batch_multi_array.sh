#!/bin/bash

# This is the batch script to run the case names given in "list_file"
# Bam files are in folders with case names
# Make sure bam files are in the correct folder
# main_outputdir will contain a folder for each case, in each folder two
# folders of "logs" and "output_folder". In each "output_folder",
# there will be a vccfiles_# folder for each segment of chromosome.

# To read all the variables needed
# Make sure to modify the config file accordingly

helpmessage=$( cat <<EOF
Usage example:

For hpc:
	$0 --config /ifs/home/c2b2/rr_lab/ar3177/bin/TOBI/gbm_pipeline/tobi_config_hpc.sh --steps B --bam 1 -s 1 -e 73 --cluster hpc

For amazon:
	$0 --config /Results/TOBI/gbm_pipeline/tobi_config_hpc.sh --steps AF -s 30 -e 30 --cluster amazon

Required Arguments:

  --config	config file
  --flag	BAF or any other combination
  --bam		1 for the first sample, 'all' for all samples
  -s		start index
  -i		end index
  --cluster amazon or hpc

Notes:

This scripts assumes that you have only one tumor bam file which is mapped to human hg19.
If this is not so, you must change the genome reference as well as the genome partition.
This scripts also assumes that Samtools, bgzip, tabix are in your PATH.
Also, this script uses hard-wired path to SnpEff.
You should change these to reflect their paths on your system.

EOF
)

# If no arguments, echo help message and quit
if [ $# == 0 ]
then
	echo "$helpmessage"
	exit;
fi

### getopts 
while [ $# -gt 0 ]
do
	if [  "$1" == "-h" -o "$1" == "-help" -o "$1" == "--help" ]; then
		shift; 
		echo "$helpmessage"
		exit;
	elif [  "$1" == "--config" ]; then
		shift; 
		config=$1; 
		shift
	elif [  "$1" == "-steps" -o "$1" == "--steps" ]; then
		shift; 
		flag=$1; 
		shift
	elif [  "$1" == "-bam" -o "$1" == "--bam" ]; then
		shift; 
		bam=$1; 
		shift
	elif [  "$1" == "-s" -o "$1" == "--start" ]; then
		shift; 
		s=$1; 
		shift
	elif [  "$1" == "-e" -o "$1" == "--end" ]; then
		shift; 
		e=$1; 
		shift
	elif [  "$1" == "-cluster" -o "$1" == "--cluster" ]; then
		shift; 
		cluster=$1; 
		shift
	else	# if unknown argument, just shift
		shift
	fi
done

source ${config}

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

if [[ bam == all ]]
then
	bam="seq 1 1 ${num_bam}"
fi

for j in ${bam} 
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
	if [[ ${cluster} == hpc ]] 
	then
		qsub -t $s-$e -V -e ${main_outputdir}/${case_name}/logs -o ${main_outputdir}/${case_name}/logs \
			-N A-${short_case_name} -cwd -l mem=${mem},time=${tim}:: \
			${script} --bam ${BAM_file} --annot_filt ${Annotation_Filtering} --patient ${case_name} --ref ${ref} \
			--steps ${flag} --outputdir ${main_outputdir}/${case_name}/output_folder --memory ${java_mem} \
			--config_file ${config}
	elif [[ ${cluster} == amazon ]]
		qsub -t $s-$e -V -e ${main_outputdir}/${case_name}/logs -o ${main_outputdir}/${case_name}/logs \
			-N A-${short_case_name} -cwd \
			${script} --bam ${BAM_file} --annot_filt ${Annotation_Filtering} --patient ${case_name} --ref ${ref} \
			--steps ${flag} --outputdir ${main_outputdir}/${case_name}/output_folder --memory ${java_mem} \
			--config_file ${config}
	else
		echo "Error: Input the correct cluster name (amazon or hpc)."
	fi
done