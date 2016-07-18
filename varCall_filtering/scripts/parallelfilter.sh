#!/bin/bash

while [ $# -gt 0 ]
do
	if [ "$1" == "-outputdir" ]; then
		shift; 
		output=$1; 
		shift
	elif [ "$1" == "-inputdir" ]; then
		shift; 
		inputdir=$1; 
		shift
	elif [ "$1" == "-whichscript" ]; then
		shift; 
		whichscript=$1; 
		shift
	elif [ "$1" == "-source_dir" ]; then
		shift; 
		source_dir=$1; 
		shift
	elif [ "$1" == "-case_name" ]; then
		shift; 
		case_name=$1; 
		shift
	else	# if unknown argument, just shift
		shift
	fi
done

#only works for human atm. change to work for all later
case $SGE_TASK_ID in
	1) c="1";;
	2) c="2";;
	3) c="3";;
	4) c="4";;
	5) c="5";;
	6) c="6";;
	7) c="7";;
	8) c="8";;
	9) c="9";;
	10) c="10";;
	11) c="11";;
	12) c="12";;
	13) c="13";;
	14) c="14";;
	15) c="15";;
	16) c="16";;
	17) c="17";;	
	18) c="18";;
	19) c="19";;
	20) c="20";;
	21) c="21";;
	22) c="22";;
	23) c="X";;
	24) c="Y";;
	25) c="MT";;
esac
echo ${source_dir}/vcf2report.py --inputdir ${inputdir} --output ${output} --case_name ${case_name}.${c}.recode

python ${source_dir}/vcf2report.py --inputdir ${inputdir} --output ${output} --case_name ${case_name}.${c}.recode

Rscript ${source_dir}${whichscript} ${output}/filter/${case_name}.${c}.recode_not_filt.tsv \
${output}/filter/${case_name}.${c}_filt_indel_techn_biol.tsv
    
    
    
    
    