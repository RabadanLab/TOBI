#########################################################################################
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

ADD DESC HERE ~

Ver. 1.2: April 12, 2016
cjmadubata & tchu modified from Alireza Roshan Ghias's code 
(Ver. 1.1: Nov 07, 2014 https://github.com/alireza202/TOBI.git TOBI)

dependencies:
	- Python 2.7.11
	- Perl v5.10.1
	- R v3.1.2
	- Java 1.7.0_25
	- samtools 0.1.19
	- bcftools 0.1.19
	- VCFtools v0.1.10.1
	- snpEff v3.6 & dbNSFP (https://sites.google.com/site/jpopgen/dbNSFP)
	- snpSift v3.6

#########################################################################################
###varCall_filtering###

inputs at each step:
	V (variant calling): indexed .bam files in a folder. Files must have .bam extension 
		and filename cannot start with a number.
	A (annotation): .vcf files in a folder. Files must have .vcf extension and filename
		cannot start with a number. If starting from this step, please format vcf to
		match bcftools output.
	F (filter): .vcf files in a folder. Files must have .vcf extension and filename 
		cannot start with a number. 
	
usage: TOBIvaf.py [-h] [--inputdir INPUTDIR] [--output OUTPUT]
                   [--config CONFIG] [--steps STEPS] [--cluster {hpc,amazon}]
                   [--debug] [--cleanup] [--ref REF] [--start START]
                   [--end END] [--snpeff SNPEFF] [--annovcf ANNOVCF]
                   [--dbnsfp DBNSFP] [--vcftype {default,TCGA}]
                   
TOBIv1.2: Tumor Only Boosting Identification of Driver Mutations. All arguments
can be specified in a config file. (See included varCall.config file as an
example). 

optional arguments:
  -h, --help            	show this help message and exit
  --inputdir INPUTDIR   	[REQUIRED] directory for bam/vcf files.
  --output OUTPUT       	[REQUIRED] output directory.
  --config CONFIG       	Config file specifying command line arguments.
            				Arguments specified in the command line overwrite config 
            				file arguments.
  --steps STEPS         	[REQUIRED] Specify which steps of pipeline to run. V:
                        	variant calling A: annotate F: filter eg. --steps AF
  --cluster {hpc,amazon}	[REQUIRED] Specify which cluster to run on. hpc: run
                        	on an SGE hpc cluster amazon: CURRENTLY UNIMPLEMENTED
  --debug               	Debug/verbose flag. Default: False
  --cleanup             	Delete temporary debug files. Default True
  --ref REF             	[REQUIRED - VCF] Reference genome file.
  --start START         	Start index used for testing. Default 1
  --end END             	End index used for testing. Default 74
  --snpeff SNPEFF       	[REQUIRED - ANNOTATE] Directory where snpEff is
  --annovcf ANNOVCF     	[REQUIRED - ANNOTATE] A comma separated list of .vcf
                        	files to annotate with.
  --dbnsfp DBNSFP       	[REQUIRED - ANNOTATE] Path to dbNSFP file
  --vcftype {default,TCGA}	Specifies vcf type specically for TCGA filtering

#########################################################################################
### machine_learning ###
Step 8. Pre-processing using R. Needs customization each time.

usage: TOBIml.py [-h] [--input INPUT] [--output OUTPUT] [--somatic SOMATIC]
                 [--log LOG] [--check_missed CHECK_MISSED] [--suffix SUFFIX]
                 [--vcftype {default,TCGA}] [--train_size TRAIN_SIZE]
                 [--verbose]
                 {preprocess,machinelearning}

TOBIv1.2: Tumor Only Boosting Identification of Driver Mutations. Machine
learning step.

positional arguments:
  {preprocess,machinelearning}
						preprocess: preprocessing step; 
						machinelearning: machine learning step
						
optional arguments:
  -h, --help            show this help message and exit
  --input INPUT         [REQUIRED] input file
  --output OUTPUT       [REQUIRED] output file
  --somatic SOMATIC     [REQUIRED] formatted file containing somatic variants
  --log LOG             Optional argument to specify a log to pipe stdout and
                        stderr to
  --check_missed CHECK_MISSED
                        [PP ARG] checking which mutations in important genes
                        are missed by filtering
  --suffix SUFFIX       [ML ARG] a label specific to this particular run (e.g.
                        <date>_<disease>)
  --vcftype {default,TCGA}
                        Specifies vcf type specically for TCGA filtering
  --train_size TRAIN_SIZE
                        [ML ARG] number of patients you want in the training
                        set.
  --verbose             verbose flag
