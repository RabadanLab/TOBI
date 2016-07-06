#!/bin/bash

while [ $# -gt 0 ]
do
	if [ "$1" == "-outputdir" ]; then
		shift; 
		outputdir=$1; 
		shift
	elif [ "$1" == "-snpeff" ]; then
		shift; 
		snpeff=$1; 
		shift
	elif [ "$1" == "-case_name" ]; then
		shift; 
		case_name=$1; 
		shift
	elif [ "$1" == "-input" ]; then
		shift; 
		input=$1; 
		shift
	elif [ "$1" == "-dbnsfp" ]; then
		shift; 
		dbnsfp=$1; 
		shift
	elif [ "$1" == "-header" ]; then
		shift; 
		header=$1; 
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


echo "[region] "$c

java -Xmx6G -jar ${snpeff}/SnpSift.jar dbnsfp ${dbnsfp} \
-v -f ${header} ${outputdir}/annotate/${case_name}.${c}.eff.vcf \
	> ${outputdir}/annotate/${case_name}.${c}.eff.vcf.tmp
	
rm ${outputdir}/annotate/${case_name}.${c}.eff.vcf

mv ${outputdir}/annotate/${case_name}.${c}.eff.vcf.tmp ${outputdir}/annotate/${case_name}.${c}.eff.vcf



