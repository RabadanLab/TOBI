#!/bin/bash

#TOBI
# checks that final TOBI file has value, if not, can check if individual pre-cursor files exist

count=0
data=$1
output=$2
echo ${data}
row=""
for i in `ls ${data}`; do
	echo ${i}
	if [ -d $data/${i} ]; then
		count=`expr $count + 1`
		for c in `seq 1 1 74`; do
			location=$data/${i}/savi_${c}/out/finalfilter.vcf
			ann=$data/${i}/output_folder/vcffiles_${c}/raw_${c}.vcf.all.annotations_filt_indel_techn_biol.tsv
			if [ ! -s $ann ]; then
				echo -e "${i}"$'\t'"${c}.ann"$'\t'"no_file" >> ${output} #DX_status.txt
			elif [ $( cat "$ann" | grep -v chrom | wc -l ) == 0 ]; then
				echo -e "${i}"$'\t'"${c}.ann"$'\t'"0" >> ${output} #DX_status.txt
			fi
		done	
	fi
done

echo $count
exit 0

