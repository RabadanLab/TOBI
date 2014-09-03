#!/bin/bash

bamdir=/ifs/scratch/c2b2/rr_lab/ar3177/TCGA/BAM_files

SNPEFF_HOME=/ifs/scratch/c2b2/rr_lab/shares/snpEff-v3.6
SNPEFF="-jar $SNPEFF_HOME/snpEff.jar -c $SNPEFF_HOME/snpEff.config"
SNPSIFT="-jar $SNPEFF_HOME/SnpSift.jar"
cosmic_variants=/ifs/scratch/c2b2/rr_lab/shares/ref/COSMIC/CosmicVariants_v66_20130725.vcf 
super_normal=/ifs/scratch/c2b2/rr_lab/jw2983/mutation_only_tumor/meganormal/219normals.cosmic.hitless100.noExactMut.mutless5000.all_samples.vcf
cbio=/ifs/scratch/c2b2/rr_lab/shares/ref/vcfs/cbio.fix.sort.vcf
dbNSFP=/ifs/scratch/c2b2/rr_lab/shares/snpEff-v3.6/data/dbNSFP2.4.txt.gz
dbNSFP_header=$(cat /ifs/home/c2b2/rr_lab/ar3177/bin/TOBI/gbm_pipeline/scripts/header_fields.txt)

java7=/ifs/scratch/c2b2/rr_lab/shares/jdk1.7.0_55/bin/java
VcfQuery=/ifs/scratch/c2b2/rr_lab/shares/vcftools/bin/vcf-query

vcfEffOnePerLine=/ifs/home/c2b2/rr_lab/ar3177/bin/TOBI/gbm_pipeline/scripts/vcfEffOnePerLine.pl
PythonParsing=/ifs/home/c2b2/rr_lab/ar3177/bin/TOBI/gbm_pipeline/scripts/parse_tsv.py
vcf2report=/ifs/home/c2b2/rr_lab/ar3177/bin/TOBI/gbm_pipeline/scripts/vcf2report.py

filter_indel=/ifs/home/c2b2/rr_lab/ar3177/bin/TOBI/gbm_pipeline/scripts/filter_indel.R
filter_techn=/ifs/home/c2b2/rr_lab/ar3177/bin/TOBI/gbm_pipeline/scripts/filter_techn.R
filter_biol=/ifs/home/c2b2/rr_lab/ar3177/bin/TOBI/gbm_pipeline/scripts/filter_biol.R


java_memory=6
#stepstr=BAF
filter=on

helpmessage=$( cat <<EOF
Usage:

~/myScripts/do_annotation_filtering.sh -i inputfile -s AF -f on

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
	elif [  "$1" == "-f" -o "$1" == "--filter" ]; then
		shift; 
		filter=$1; 
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
echo "[steps] "$stepstr
echo "[filter] "$filter


if [[ $stepstr == *A* ]]
then
	echo "[STEP2] annotation"
	### SnpEff

	${java7} -Xmx${java_memory}G $SNPEFF GRCh37.71 -noStats -v -lof \
		-canon -no-downstream -no-intergenic -no-intron -no-upstream -no-utr ${outputdir}/${inputfile} \
		> ${outputdir}/${inputfile}.eff.vcf
	${java7} -Xmx${java_memory}G $SNPSIFT annotate $SNPEFF_HOME/dbSnp138.vcf -v ${outputdir}/${inputfile}.eff.vcf \
		> ${outputdir}/${inputfile}.eff.dbSNP.vcf
	${java7} -Xmx${java_memory}G $SNPSIFT annotate $SNPEFF_HOME/clinvar_20140303.vcf -v ${outputdir}/${inputfile}.eff.dbSNP.vcf \
		> ${outputdir}/${inputfile}.eff.dbSNP.clinvar.vcf
	${java7} -Xmx${java_memory}G $SNPSIFT annotate $cosmic_variants -v ${outputdir}/${inputfile}.eff.dbSNP.clinvar.vcf \
		> ${outputdir}/${inputfile}.eff.dbSNP.clinvar.cosmic.vcf
	${java7} -Xmx${java_memory}G  $SNPSIFT annotate $super_normal -v ${outputdir}/${inputfile}.eff.dbSNP.clinvar.cosmic.vcf \
		> ${outputdir}/${inputfile}.eff.dbSNP.clinvar.cosmic.meganormal1.vcf
#	${java7} -Xmx${java_memory}G  $SNPSIFT annotate $super_normal2 -v ${outputdir}/${inputfile}.eff.dbSNP.clinvar.cosmic.meganormal1.vcf \
#		> ${outputdir}/${inputfile}.eff.dbSNP.clinvar.cosmic.meganormal1.2.vcf
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
	
	case_name=$(echo ${outputdir} | awk -F'/' '{print $(NF-2)}')
	echo ${case_name}
	
	# Replacing # and ' characters
	cat ${outputdir}/${inputfile}.all.annotations.rep.vcf | \
		${vcf2report} 0 | \
		${PythonParsing} ${case_name} | \
		sed -e "s/#//g;s/'//g" > \
		${outputdir}/${inputfile}.all.annotations_not_filt.tsv
	
	echo "Applying filters"
	R --slave -f ${filter_indel} --args \
		${outputdir}/${inputfile}.all.annotations_not_filt.tsv \
		${outputdir}/${inputfile}.all.annotations_filt_indel.tsv
	R --slave -f ${filter_techn} --args \
		${outputdir}/${inputfile}.all.annotations_filt_indel.tsv \
		${outputdir}/${inputfile}.all.annotations_filt_indel_techn.tsv
	R --slave -f ${filter_biol} --args \
		${outputdir}/${inputfile}.all.annotations_filt_indel_techn.tsv \
		${outputdir}/${inputfile}.all.annotations_filt_indel_techn_biol.tsv
	
	echo "[date 4] "`date`;
fi

echo [deltat]
echo [End] `date`