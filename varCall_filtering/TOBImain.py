import argparse
import os
import subprocess
import scripts.helpers as helpers
import re

def get_arg():
    """Get Arguments"""

    prog_description = """TOBIv1.2 ADD DESC HERE"""
    parser = argparse.ArgumentParser(description=prog_description)
    
    #main arguments
    parser.add_argument(
        '--inputdir',
        help = "directory for bam/vcf files"
        )
    parser.add_argument(
        '--output',
        help = "output directory"
        )
    parser.add_argument(
        '--config',
        help = """optional config file specifying command line arguments.
            Arguments specified in the config file will overwrite command
            line arguments."""                
        )
    parser.add_argument(
        '--steps',
        type = str,
        help = """Specify which steps of pipeline to run. (REQUIRED)
            B: variant calling
            A: annotate 
            F: filter
            eg. --steps AF"""
        )
    parser.add_argument(
        '--cluster',
        choices = ['hpc','amazon'],
        help = """Specify which cluster to run on. (REQUIRED)
            hpc: run on an SGE hpc cluster
            amazon: run on amazon's clusters (NOT IMPLEMENTED)"""
        )
    parser.add_argument(
        '--debug',
        default = False,
        action = 'store_true',
        help = "Debug flag. NOT IMPLEMENTED"                
        )
    
    #arguments for varcall
    parser.add_argument(
        '--ref',
        help = "Reference genome file. Required for VCF calling"
        )
    parser.add_argument(
        '--start',
        type = int,
        default = 1,
        help = "Start index. Default 1"
        )
    parser.add_argument(
        '--end',
        type = int,
        default = 74,
        help = "End index. Default 74"
        )
    
    #arguments for annotate
    parser.add_argument(
        '--snpeff',
        help = "directory where snpEff is"
        )
    parser.add_argument(
        '--annovcf',
        help = "comma separated list of .vcf files for annotation."
        )
    parser.add_argument(
        '--dbnsfp',
        help = "path to dbNSFP file"
        )
    
    software = os.path.dirname(os.path.realpath(__file__))
    parser.add_argument(
        "--scripts",
        default=software,
        help="""location of scripts dir (directory where this script resides 
        - use this option only if qsub-ing with the Oracle Grid Engine)"""
        )
    parser.add_argument(
        "--verbose", 
        action="store_true", 
        help="verbose mode: echo commands, etc (default: off)"
        )

    args = parser.parse_args()
    # add key-value pairs to the args dict
    vars(args)['cwd'] = os.getcwd()

    return args

def main():
    args = get_arg()
    
    if args.config:
        args = helpers.parse_config(args)
    if args.debug:
        print(args)
        
    #vcf calling
    if "B" in args.steps or "b" in args.steps:
        input_filenames = helpers.get_filenames(args.inputdir,"bam") 

        if not os.path.exists(args.output):
            os.makedirs(args.output + "/logs") 
                
        vcf_call(input_filenames,args)
        #set inputdir as vcf_call's output
        args.inputdir = args.output

    #annotate
    if "A" in args.steps or "a" in args.steps:
        input_filenames = helpers.get_filenames(args.inputdir,"vcf") 
        #for each case name/bam file, make output directories
        if not os.path.exists(args.output):
            os.makedirs(args.output + "/logs")
                
        annotate(input_filenames, args)
    
    #filter
    if bool(re.search(".*F.*",args.steps)):
        pass

def vcf_call(input_filenames, args):
    source_dir = os.path.dirname(os.path.realpath(__file__))
    for case_name in input_filenames:
        #run mpileup
        proc = subprocess.Popen(
            helpers.mpileup_cmdgen(args,case_name,source_dir), 
            shell=True,
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE
            ) 
        #wait for job completion 
        proc.wait()
        
        #concat raw_n.vcf files into new file
        proc = subprocess.Popen(
            helpers.vcf_concat_cmdgen(args,case_name), 
            shell=True,
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE
            ) 
        proc.wait()
        
        proc = subprocess.Popen(
            "vcf-sort " + case_name + ".vcf.gz > " 
                + args.output+"/"+case_name+".sorted.vcf.gz", 
            shell=True,
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE
            ) 
        proc.wait()
        os.remove(args.output+"/"+case_name+".vcf.gz")
        os.rename(args.output+"/"+case_name+".sorted.vcf.gz",
                  args.output+"/"+case_name+".vcf.gz")
        
        #purge raw_n.vf files
        #ADD FLAG HERE
        if not args.debug:
            helpers.purge(args.output, "raw_\d*\.vcf")
        
    return 

def annotate(input_filenames, args):
    #source_dir = os.path.dirname(os.path.realpath(__file__))
    annovcf = args.annovcf.replace("\n","").split(',')
    dbnsfp_header="SIFT_score,SIFT_converted_rankscore,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred,LRT_score,LRT_converted_rankscore,LRT_pred,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MutationAssessor_score,MutationAssessor_rankscore,MutationAssessor_pred,FATHMM_score,FATHMM_rankscore,FATHMM_pred,RadialSVM_score,RadialSVM_rankscore,RadialSVM_pred,LR_score,LR_rankscore,LR_pred,Reliability_index,CADD_raw,CADD_raw_rankscore,CADD_phred,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,phyloP46way_primate,phyloP46way_primate_rankscore,phyloP46way_placental,phyloP46way_placental_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phastCons46way_primate,phastCons46way_primate_rankscore,phastCons46way_placental,phastCons46way_placental_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,SiPhy_29way_pi,SiPhy_29way_logOdds,SiPhy_29way_logOdds_rankscore,LRT_Omega,UniSNP_ids,1000Gp1_AC,1000Gp1_AF,1000Gp1_AFR_AC,1000Gp1_AFR_AF,1000Gp1_EUR_AC,1000Gp1_EUR_AF,1000Gp1_AMR_AC,1000Gp1_AMR_AF,1000Gp1_ASN_AC,1000Gp1_ASN_AF,ESP6500_AA_AF,ESP6500_EA_AF"
    for case_name in input_filenames:      
        if args.cluster == "hpc":
            #snpEff annotate
            proc = subprocess.Popen(
                helpers.snpeff_cmdgen(args,case_name),
                shell=True,
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE
                ) 
            proc.wait()

            #snpSIFT annotate for each vcf provided
            for vcf in annovcf:
                proc = subprocess.Popen(
                    helpers.snpsift_cmdgen(args,case_name,vcf),
                    shell=True,
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.PIPE
                    )
                proc.wait()
                #hacky way to create temp file fix this later
                
                os.remove(args.output+"/"+case_name+".eff.vcf")
                os.rename(args.output+"/"+case_name+".eff.vcf.tmp",
                          args.output+"/"+case_name+".eff.vcf"
                          )
            #snpSIFT dbnsfp
            proc = subprocess.Popen(
                helpers.snpdbnsfp_cmdgen(args,case_name,args.dbnsfp,dbnsfp_header),
                shell=True,
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE
                )
            proc.wait()

def filter(input_filenames,args):
    pass
            
if __name__ == "__main__":
    main()
