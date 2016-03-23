######cjm edits for TOBI####

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

Input: vcf files as described in "list file" from step0; maf file to generate "list file" #WXS bam files in `list_file`.

On Amazon, go to a volume (e.g. /Results) and:

	git clone https://github.com/alireza202/TOBI.git TOBI

Ver. 1.2: Mar 22, 2016
cjmadubata modified from Alireza Roshan Ghias's code (Ver. 1.1: Nov 07, 2014)

########################################################################################

Step0. Run: varCall_filtering/maf_to_tobi_som_var_list.sh {path_to_maf} {output_true_somatic} 
	to generate list of  true somatic varints {output_true_somatic} 
	and list of case names patient.{output_true_somatic} (aka list_file)

Step1. Update tobi_config file for your samples 

VCF files are assumed to be in same folder with different sample names at beginning of file
followed by a common suffix
	e.g. path/to/vcf/files/TCGA-XX-XXXX-*suffix
	each TCGA-XX-XXXX represents one patient and one line in "list_file"

Step 2. Check the batch run file
	varCall_filtering/batch_multi_array.LGG_TCGA.sh_suffix 

Step 5. Run all samples using -AF flags:
	varCall_filtering/batch_multi_array.LGG_TCGA.sh_suffix --config /path/to/config_file 
	--steps AF --bam all -s 0 -e 0 --cluster hpc --filter on --
	# --bam -s -e flags are not used for VCF analysis

Step 6. Check the final file sizes, before and after filtering, using the script below:

	check_final_run.sh
	
If some cases where not done, find the problem, and run them again.

Step 7. Merge all final files for each case, and then all together using the script 
below:

	merge_all_tsvs.sh

Step 8. Pre-processing using R. Needs customization each time.

	machine_learning/pre_processing.R
