#!/bin/bash

#TOBI
# checks that final TOBI file has value, if not, can check if individual pre-cursor files exist

count=0
data=$1
list_file=$2
echo ${data}"; "${list_file}
row=""
for i in `ls ${data} | grep -f ${list_file}`; do
#	echo ${i}
	if [ -d $data/${i} ]; then
		count=`expr $count + 1`
			#ann=$data/${i}/output_folder/vcffiles_${c}/raw_${c}.vcf.all.annotations_filt_indel_techn_biol.tsv
			ann=$data/${i}/output_folder/${i}.oxoG.snp.capture.tcga.vcf.all.annotations_filt_indel_techn_biol.tsv
		wc -l "$ann"
		awk '{print NF}' "$ann" | sort | uniq -c | sort -nr
		cut -f1 "$ann" | sort | uniq -c | sort -nr 
			if [ ! -s $ann ]; then
				echo -e "${i}"$'\t'"filt_indel_techn_biol.tsv"$'\t'"no_file" >> DX_status.txt
			elif [ $( cat "$ann" | grep -v chrom | wc -l ) == 0 ]; then
				echo -e "${i}"$'\t'"filt_indel_techn_biol.tsv"$'\t'"0" >> DX_status.txt
			fi
	fi
done

echo $count
exit 0

