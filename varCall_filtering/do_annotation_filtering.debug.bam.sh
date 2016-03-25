#!/bin/bash

### default variables
java_memory=6
filter=on
debug=0                                         # bool to turn on debugging

###
helpmessage=$( cat <<EOF
Usage:

do_annotation_filtering.sh -i inputfile -s AF -f on

Required Arguments:

  - inputfile including path

Optional Arguments:
  - steps (AF)
  - filter (on or off)

EOF
)

# If no arguments, echo help message and quit
if [ $# == 0 ]
then
	echo "$helpmessage"
	exit;
fi

while [ $# -gt 0 ]
do
	if [  "$1" == "-h" -o "$1" == "-help" -o "$1" == "--help" ]; then
		shift;
		echo "$helpmessage"
		exit;
	elif [  "$1" == "-i" -o "$1" == "--input-file" ]; then
		shift; 
		inputfile=$(basename $1)
		outputdir=$(dirname $1)
		shift
	elif [  "$1" == "-s" -o "$1" == "--steps" ]; then
		shift; 
		stepstr=$1; 
		shift
	elif [  "$1" == "--config_file" ]; then
		shift; 
		config=$1; 
		shift
	elif [  "$1" == "-f" -o "$1" == "--filter" ]; then
		shift; 
		filter=$1; 
		shift
	elif [  "$1" == "--memory" ]; then
		shift; 
		java_memory=$1; 
		shift
	elif [  "$1" == "-debug" -o "$1" == "--debug" ]; then 
		shift; 
		debug=$1; 
		shift
	else	
		# if unknown argument, just shift
		shift
	fi
done

echo "[start]"
echo "[pwd] "`pwd`
echo "[date] "`date`
echo "[input_file] "$inputfile
echo "[output_dir] "$outputdir
echo "[java_memory] "$java_memory
echo "[steps] "$stepstr
echo "[filter] "$filter
echo "[config file] "$config

# Reading in the config file
source ${config}

if [[ $stepstr == *A* ]]
then
	echo "[STEP2] annotation"
	### SnpEff
	dbNSFP_header=$(cat ${dbNSFP_header})

	${java7} -Xmx${java_memory}G $SNPEFF GRCh37.71 -noStats -v -lof \
		-canon -no-downstream -no-intergenic -no-intron -no-upstream -no-utr ${outputdir}/${inputfile} \
		> ${outputdir}/${inputfile}.eff.vcf
	${java7} -Xmx${java_memory}G $SNPSIFT annotate $dbSnp138 -v ${outputdir}/${inputfile}.eff.vcf \
		> ${outputdir}/${inputfile}.eff.dbSNP.vcf
	${java7} -Xmx${java_memory}G $SNPSIFT annotate $clinvar -v ${outputdir}/${inputfile}.eff.dbSNP.vcf \
		> ${outputdir}/${inputfile}.eff.dbSNP.clinvar.vcf
	${java7} -Xmx${java_memory}G $SNPSIFT annotate $cosmic_variants -v ${outputdir}/${inputfile}.eff.dbSNP.clinvar.vcf \
		> ${outputdir}/${inputfile}.eff.dbSNP.clinvar.cosmic.vcf
	${java7} -Xmx${java_memory}G  $SNPSIFT annotate $super_normal -v ${outputdir}/${inputfile}.eff.dbSNP.clinvar.cosmic.vcf \
		> ${outputdir}/${inputfile}.eff.dbSNP.clinvar.cosmic.meganormal1.vcf
	${java7} -Xmx${java_memory}G  $SNPSIFT annotate $cbio -v ${outputdir}/${inputfile}.eff.dbSNP.clinvar.cosmic.meganormal1.vcf \
		> ${outputdir}/${inputfile}.eff.dbSNP.clinvar.cosmic.meganormal1.cbio.vcf
	${java7} -Xmx${java_memory}G  $SNPSIFT dbnsfp $dbNSFP -v -f $dbNSFP_header ${outputdir}/${inputfile}.eff.dbSNP.clinvar.cosmic.meganormal1.cbio.vcf \
		> ${outputdir}/${inputfile}.eff.dbSNP.clinvar.cosmic.meganormal1.cbio.dbNSFP.vcf
	
	# To separate lines with multiple EFF into different lines
	echo "separate lines with multiple EFF into different lines"
	cat ${outputdir}/${inputfile}.eff.dbSNP.clinvar.cosmic.meganormal1.cbio.dbNSFP.vcf | ${vcfEffOnePerLine} > ${outputdir}/${inputfile}.all.annotations.vcf
	
	echo "[date 3] "`date`;
fi

if [[ $stepstr == *F* ]]
then
	echo "[STEP3] filtering"
	
	# Replacing ++ characters
	cat ${outputdir}/${inputfile}.all.annotations.vcf | sed 's/GERP++/GERP/g' > ${outputdir}/${inputfile}.all.annotations.rep.vcf
	
	echo "python vcf2report and parse_tsv"
	
	case_name=$(echo ${outputdir} | awk -F'/' '{print $(NF-2)}')  #cjm, change for appropriate naming depth
	echo ${case_name}
	
	# Replacing # and ' characters
	cat ${outputdir}/${inputfile}.all.annotations.rep.vcf | \
		${vcf2report} 0 | \
		${PythonParsing} ${case_name} | \
		sed -e "s/#//g;s/'//g" > \
		${outputdir}/${inputfile}.all.annotations_not_filt.tsv
	
	if [[ ${filter} == on ]]
	then
		echo "Applying filters"
		# Applying all filters together
		Rscript ${filter_indel_techn_biol} \
			${outputdir}/${inputfile}.all.annotations_not_filt.tsv \
			${outputdir}/${inputfile}.all.annotations_filt_indel_techn_biol.tsv
		gzip ${outputdir}/${inputfile}.all.annotations_not_filt.tsv #cjm, making more space
		gzip ${outputdir}/${inputfile} #cjm 
	fi
	
	echo "[date 4] "`date`;
fi
#debug option currently not implemented... #cjm
echo "debug option is " $debug
if [ $debug -eq 0 ]; then #cjm
	echo "removing extra files"
        rm -f ${outputdir}/${inputfile}*eff.vcf #cjm
        rm -f ${outputdir}/${inputfile}*dbSNP.vcf #cjm
        rm -f ${outputdir}/${inputfile}*clinvar.vcf #cjm
        rm -f ${outputdir}/${inputfile}*cosmic.vcf #cjm
        rm -f ${outputdir}/${inputfile}*meganormal1.vcf #cjm
        rm -f ${outputdir}/${inputfile}*cbio.vcf #cjm
        rm -f ${outputdir}/${inputfile}*indel.tsv #cjm
        rm -f ${outputdir}/${inputfile}*techn.tsv #cjm
        rm -f ${outputdir}/${inputfile}*dbNSFP.vcf #cjm
        rm -f ${outputdir}/${inputfile}*all.annotations.rep.vcf #cjm, this is first file generated in step F
 #cjm
fi #cjm

echo [deltat]
echo [End] `date`
