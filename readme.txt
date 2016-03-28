######cjm edits for TOBI BAM files Mar 28, 2016 ####

########################################################################################
TTTTTTTTTTTTTTTTTTTTTTT          OOOOOOOOO          BBBBBBBBBBBBBBBBB        IIIIIIIIII
T:::::::::::::::::::::T        OO:::::::::OO        B::::::::::::::::B       I::::::::I
T:::::::::::::::::::::T      OO:::::::::::::OO      B::::::BBBBBB:::::B      I::::::::I
T:::::TT:::::::TT:::::T     O:::::::OOO:::::::O     BB:::::B     B:::::B     II::::::II
TTTTTT  T:::::T  TTTTTT     O::::::O   O::::::O       B::::B     B:::::B       I::::I  
        T:::::T             O:::::O     O:::::O       B::::B     B:::::B       I::::I  
        T:::::T             O:::::O     O:::::O       B::::BBBBBB:::::B        I::::I  
        T:::::T             O:::::O     O:::::O       B:::::::::::::BB         I::::I  
        T:::::T             O:::::O     O:::::O       B::::BBBBBB:::::B        I::::I  
        T:::::T             O:::::O     O:::::O       B::::B     B:::::B       I::::I  
        T:::::T             O:::::O     O:::::O       B::::B     B:::::B       I::::I  
        T:::::T             O::::::O   O::::::O       B::::B     B:::::B       I::::I  
      TT:::::::TT           O:::::::OOO:::::::O     BB:::::BBBBBB::::::B     II::::::II
      T:::::::::T            OO:::::::::::::OO      B:::::::::::::::::B      I::::::::I
      T:::::::::T              OO:::::::::OO        B::::::::::::::::B       I::::::::I
      TTTTTTTTTTT                OOOOOOOOO          BBBBBBBBBBBBBBBBB        IIIIIIIIII

TOBI: Tumor Only Boosting Identification of Driver Mutations

Input: BAM files as described in "list file" from step0; maf file to generate "list file" #WXS bam files in `list_file`.

Ver. 1.2.BAM: Mar 28, 2016
cjmadubata modified from Alireza Roshan Ghias's code (Ver. 1.1: Nov 07, 2014 https://github.com/alireza202/TOBI.git TOBI)

########################################################################################
###varCall_filtering###
Step0. Run: varCall_filtering/maf_to_tobi_som_var_list.sh {path_to_maf} {output_true_somatic} 
	to generate list of  true somatic varints {output_true_somatic} 
	and list of case names patient.{output_true_somatic} (aka list_file)

Step1. Update tobi_config file for your samples, particularly location of list_file

VCF files are assumed to be in same folder with different sample names at beginning of file
followed by a common suffix
	e.g. path/to/vcf/files/TCGA-XX-XXXX-*suffix
	each TCGA-XX-XXXX represents one patient and one line in "list_file"

Step 2. Check the batch run file
	varCall_filtering/batch_multi_array.bam_input.sh

Step 5. Run all samples using -BAF flags:
	varCall_filtering/batch_multi_array.bam_input.sh --config /path/to/config_file 
	--steps AF --bam all -s 0 -e 0 --cluster hpc --filter on --vcfsuf oxoG.snp.capture.tcga.vcf
	# --bam -s -e flags are not used for VCF analysis

Step 6. Check number of lines in final file sizes, using the script below:
	varCall_filtering/tobi_out_empty.sh {main_outputdir_from_config} {list_file}	
	
If some cases did not run, find the problem, and run them again.

Step 7. Merge all final files for each case, and then all together using the script 
below:

	varCall_filtering/merge_all_tsvs.sh {list_file} {main_outputdir_from_config} {/path/to/output/table} {label_for_this_TOBI_run} 	

Step 7a. If using COSMIC after v66, need to replace "cnt" column with "cosmic_nsamp" for next TOBI steps
	sed -i '1s/cnt/cosmic_nsamp/' {/path/to/output/table}

### machine_learning ###
Step 8. Pre-processing using R. Needs customization each time.

	machine_learning/pre_processing_glioma.151230.R
