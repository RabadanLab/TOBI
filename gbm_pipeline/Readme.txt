###########################################################################################                                                                    
TTTTTTTTTTTTTTTTTTTTTTT     OOOOOOOOO     BBBBBBBBBBBBBBBBB   IIIIIIIIII
T:::::::::::::::::::::T   OO:::::::::OO   B::::::::::::::::B  I::::::::I
T:::::::::::::::::::::T OO:::::::::::::OO B::::::BBBBBB:::::B I::::::::I
T:::::TT:::::::TT:::::TO:::::::OOO:::::::OBB:::::B     B:::::BII::::::II
TTTTTT  T:::::T  TTTTTTO::::::O   O::::::O  B::::B     B:::::B  I::::I  
        T:::::T        O:::::O     O:::::O  B::::B     B:::::B  I::::I  
        T:::::T        O:::::O     O:::::O  B::::BBBBBB:::::B   I::::I  
        T:::::T        O:::::O     O:::::O  B:::::::::::::BB    I::::I  
        T:::::T        O:::::O     O:::::O  B::::BBBBBB:::::B   I::::I  
        T:::::T        O:::::O     O:::::O  B::::B     B:::::B  I::::I  
        T:::::T        O:::::O     O:::::O  B::::B     B:::::B  I::::I  
        T:::::T        O::::::O   O::::::O  B::::B     B:::::B  I::::I  
      TT:::::::TT      O:::::::OOO:::::::OBB:::::BBBBBB::::::BII::::::II
      T:::::::::T       OO:::::::::::::OO B:::::::::::::::::B I::::::::I
      T:::::::::T         OO:::::::::OO   B::::::::::::::::B  I::::::::I
      TTTTTTTTTTT           OOOOOOOOO     BBBBBBBBBBBBBBBBB   IIIIIIIIII

TOBI: Tumor Only Boosting Identification of Driver Mutations

Input: WXS bam files in `list_file`.

Alireza Roshan Ghias
Ver. 1.1: Nov 07, 2014

###########################################################################################

Step1. Update tobi_config file for your problem.

Bam files are assumed to be in separate folders with the main folder named after the case
and embedded in another folder, like:

	TCGA-00-0000/*/*.bam

Step 2. Update the list_file.txt file to contain the case names.

Step 3. Check the batch run file to contain the right information.

Step 4. Run one sample to check things are in the right place. Start with B flag. 
Run the batch file. Check the vcf file, and if needed, the log file. Continue with A, 
and F separately.

Step 5. Run all samples.

Step 6. Check the final file sizes, before and after filtering, using the script below:

	check_final_run.sh
	
If some cases where not done, find the problem, and run them again.

Step 7. Merge all final files for each case, and then all together using the script below:

	merge_all_tsvs.sh

Step 8. Pre-processing using R. Needs customization each time.

	/Users/ar3177/Dropbox/Columbia/RabadanLab/Projects/GBM/GBM/matching-somatic-in-tumors-2.R

Step 9. Then check if the somatic and the total number of mutations are within the range for each sample. Discard the Abnormal samples. Put their name in the Abnormal.txt file.

Step 10. Export the filtered file.