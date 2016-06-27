#!/usr/bin/env python

import argparse
import os
import subprocess
import varCall_filtering.scripts.helpers as helpers
import varCall_filtering.scripts.vcf2report as vcf2report
import varCall_filtering.scripts.parse_tsv as parse_tsv
import sys


def get_arg():
    #Get Arguments

    prog_description = """TOBIv1.2: 
        Tumor Only Boosting Identification of Driver Mutations
        All arguments can be specified in a config file. (See
        included varCall.config file as an example)."""
    parser = argparse.ArgumentParser(description=prog_description)
    
    #main arguments
    parser.add_argument(
        '--inputdir',
        help = "[REQUIRED] directory for bam/vcf files. "
        )
    parser.add_argument(
        '--output',
        help = "[REQUIRED] output directory."
        )
    parser.add_argument(
        '--config',
        help = """config file specifying command line arguments.
            Arguments specified in the command line overwrite config 
            file arguments."""                
        )
    parser.add_argument(
        '--steps',
        type = str,
        help = """[REQUIRED] Specify which steps of pipeline to run.
            V: variant calling
            A: annotate 
            F: filter
            M: merge
            eg. --steps AF"""
        )
    parser.add_argument(
        '--cluster',
        choices = ['hpc','amazon'],
        help = """[REQUIRED] Specify which cluster to run on.
            hpc: run on an SGE hpc cluster
            amazon: CURRENTLY UNIMPLEMENTED"""
        )
    parser.add_argument(
        '--debug',
        default = False,
        action = 'store_true',
        help = "Debug/verbose flag. Default: False"                
        )
    parser.add_argument(
        '--cleanup',
        default = True,
        action = 'store_false',
        help = "Delete temporary debug files. Default True"                
        )
    
    #arguments for varcall
    parser.add_argument(
        '--ref',
        help = "[REQUIRED - VCF] Reference genome file."
        )
    parser.add_argument(
        '--start',
        type = int,
        default = 1,
        help = "Start index used for testing. Will not work in config. Default 1"
        )
    parser.add_argument(
        '--end',
        type = int,
        default = 74,
        help = "End index used for testing. Will not work in config. Default 74"
        )
    
    #arguments for annotate
    parser.add_argument(
        '--snpeff',
        help = "[REQUIRED - ANNOTATE] Directory where snpEff is"
        )
    parser.add_argument(
        '--annovcf',
        help = """[REQUIRED - ANNOTATE] A comma separated list of .vcf files 
            to annotate with."""
        )
    parser.add_argument(
        '--dbnsfp',
        help = "[REQUIRED - ANNOTATE] Path to dbNSFP file"
        )
    
    #arguments for filter
    parser.add_argument(
        "--vcftype", 
        choices = ['default','TCGA'], 
        help="Specifies vcf type specically for TCGA filtering"
        )
    
    #arguments for merge
    parser.add_argument(
        "--mergename",  
        help="[REQUIRED - MERGE] Name for final merged file"
        )
    
    #print help if no arguments are given
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
        
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
        
    helpers.check_main_args(args)
        
    #vcf calling
    if "V" in args.steps or "v" in args.steps:
        helpers.check_varcall_args(args)
        input_filenames = helpers.get_filenames(args.inputdir,"bam") 

        if not os.path.exists(args.output+"/vcfcall"):
            os.makedirs(args.output + "/vcfcall/logs") 
        
        helpers.multithread(vcf_call,args,input_filenames)
        
        args.inputdir = args.output +"/vcfcall"

    #annotate
    if "A" in args.steps or "a" in args.steps:
        helpers.check_anno_args(args)
        input_filenames = helpers.get_filenames(args.inputdir,"vcf") 
        #for each case name/bam file, make output directories
        if not os.path.exists(args.output+"/annotate"):
            os.makedirs(args.output + "/annotate/logs")
        
        helpers.multithread(annotate,args,input_filenames)
        #set inputdir as annotate's output
        args.inputdir = args.output +"/annotate"
    
    #filter
    if "F" in args.steps or "f" in args.steps:
        helpers.check_filt_args(args)
        input_filenames = helpers.get_filenames(args.inputdir,"vcf") 
        #for each case name/bam file, make output directories
        if not os.path.exists(args.output+"/filter"):
            os.makedirs(args.output + "/filter/logs")
            
        helpers.multithread(filter_vcf,args,input_filenames)
        #set inputdir as filter's output
        args.inputdir = args.output +"/filter"
    
    #merge
    if "M" in args.steps or "m" in args.steps:
        helpers.check_merge_args(args)
        input_filenames = helpers.get_filenames(args.inputdir,"tsv")
        if len(input_filenames) <= 1:
            sys.exit("[ERROR] One or less files to merge.")
        if not os.path.exists(args.output+"/final"):
            os.makedirs(args.output + "/final/logs")
        merge_tsv(input_filenames,args)
        

def vcf_call(case_name, args):
    source_dir = os.path.dirname(os.path.realpath(__file__))
    #run mpileup
    helpers.runShellCmd(helpers.mpileup_cmdgen(args,case_name,source_dir))
    
    #concat raw_n.vcf files into new file
    helpers.runShellCmd(helpers.vcf_concat_cmdgen(args,case_name))
    
    #sort vcf file
    helpers.runShellCmd(
            "vcf-sort -c " + args.output+"/vcfcall/"+case_name + ".vcf  > " 
            + args.output+"/vcfcall/"+case_name+".sorted.vcf")
    
    #hacky way to create temp file fix this later
    os.remove(args.output+"/vcfcall/"+case_name+".vcf")
    os.rename(args.output+"/vcfcall/"+case_name+".sorted.vcf",
              args.output+"/vcfcall/"+case_name+".vcf")
    #purge raw_n.vf files
    if args.cleanup:
        helpers.purge(args.output+"/vcfcall/"+case_name, "raw_\d*\.vcf")

def annotate(case_name, args):
    source_dir = os.path.dirname(os.path.realpath(__file__))
    annovcf = args.annovcf.replace("\n","").split(',')
    dbnsfp_header="SIFT_score,SIFT_converted_rankscore,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred,LRT_score,LRT_converted_rankscore,LRT_pred,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MutationAssessor_score,MutationAssessor_rankscore,MutationAssessor_pred,FATHMM_score,FATHMM_rankscore,FATHMM_pred,RadialSVM_score,RadialSVM_rankscore,RadialSVM_pred,LR_score,LR_rankscore,LR_pred,Reliability_index,CADD_raw,CADD_raw_rankscore,CADD_phred,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,phyloP46way_primate,phyloP46way_primate_rankscore,phyloP46way_placental,phyloP46way_placental_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phastCons46way_primate,phastCons46way_primate_rankscore,phastCons46way_placental,phastCons46way_placental_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,SiPhy_29way_pi,SiPhy_29way_logOdds,SiPhy_29way_logOdds_rankscore,LRT_Omega,UniSNP_ids,1000Gp1_AC,1000Gp1_AF,1000Gp1_AFR_AC,1000Gp1_AFR_AF,1000Gp1_EUR_AC,1000Gp1_EUR_AF,1000Gp1_AMR_AC,1000Gp1_AMR_AF,1000Gp1_ASN_AC,1000Gp1_ASN_AF,ESP6500_AA_AF,ESP6500_EA_AF" 
    if args.cluster == "hpc":
        #split vcf into multiple vcfs
        for chrom in range(1,23) + ['X', 'Y', 'MT']:
            helpers.runShellCmd("vcftools --recode --recode-INFO-all --vcf "
                + args.inputdir + "/" + case_name +".vcf" 
                +" --out " + args.output +"/annotate/"+ case_name +"."+ str(chrom) +" --chr " + str(chrom)) 
        #snpEff annotate arrya job
        helpers.runShellCmd(helpers.snpeffarray_cmdgen(args,case_name,source_dir)) 
        #snpSIFT annotate for each vcf provided
        for vcf in annovcf:
            helpers.runShellCmd(helpers.snpsift_cmdgen(args,case_name,vcf))
            #hacky way to create & replace temp file. fix this later
            os.remove(args.output+"/annotate/"+case_name+".eff.vcf")
            os.rename(args.output+"/annotate/"+case_name+".eff.vcf.tmp",
                      args.output+"/annotate/"+case_name+".eff.vcf"
                      )
        #snpSIFT dbnsfp
        helpers.runShellCmd(
            helpers.snpdbnsfp_cmdgen(args,case_name,args.dbnsfp,dbnsfp_header)
            )
        if args.cleanup:
            os.remove(args.output+"/annotate/"+case_name+".eff.vcf")
            
        #split effects into one effect per line
        helpers.runShellCmd(helpers.oneEff_cmdgen(args,case_name,source_dir))
        
        if args.cleanup:
            os.remove(args.output+"/annotate/"+case_name+".eff.all.vcf")

def filter_vcf(case_name,args):
    if args.debug:
        print("[Preprocessing file]")
    source_dir = os.path.dirname(os.path.realpath(__file__))
    #read in file
    inputfile = open(args.inputdir+"/"+case_name+".vcf","r")
    case_file = inputfile.read()
    inputfile.close()
    #find GERP++ and replace with GERP
    case_file = case_file.replace('GERP++','GERP')
    #vcf2report
    case_file = vcf2report.convert(case_file)
    #parse_tsv case_name
    case_file = parse_tsv.convert(case_file,case_name)
    #get rid of # and ' characters
    case_file = case_file.replace("#","")
    case_file = case_file.replace("'","")
    #write to new file
    outputfile = open(args.output+"/filter/"+case_name+"_not_filt.tsv","w")
    outputfile.write(case_file)
    outputfile.close()
    #apply R script filters
    if args.vcftype == "default":
        script = "/varCall_filtering/scripts/filter_indel_techn_biol.pediatric.R "
    elif args.vcftype == "TCGA":
        script ="/varCall_filtering/scripts/filter_indel_techn_biol.R "
    else:
        sys.exit("[ERROR]: Invalid '--vcftype' argument: " + args.vcftype)
    
    cmd = "qsub -V -b y -sync y -N " + case_name \
        + " -l mem=10G,time=2:: -pe smp 2 " \
        + "-e " + args.output + "/filter/logs/"+ case_name +".e " \
        + "-o " + args.output+"/filter/logs/"+case_name+".o " \
        +"Rscript " + source_dir + script \
        + args.output+"/filter/"+case_name+"_not_filt.tsv " \
        + args.output+"/filter/"+case_name+"_filt_indel_techn_biol.tsv" 
    if args.debug:
        print(cmd)        
    proc = subprocess.Popen(
        cmd,
        shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
    if args.debug:
        while True:
            nextline = proc.stdout.readline()
            if nextline == '' and proc.poll() != None:
                break
            sys.stdout.write(nextline)
            sys.stdout.flush()
        while True:
            nextline = proc.stderr.readline()
            if nextline == '' and proc.poll() != None:
                break
            sys.stderr.write(nextline)
            sys.stderr.flush()
    proc.wait()
    if args.cleanup:
        os.remove(args.output+"/filter/"+case_name+"_not_filt.tsv")

def merge_tsv(input_filenames,args):
    first = True
    for case_name in input_filenames:
        #save first file into output
        inputfile = open(args.inputdir+"/"+case_name+".tsv","r")
        case_file = inputfile.read()
        inputfile.close()
        if first:
            first = False
            outputfile = open(args.output+"/final/"+args.mergename+".merged.tsv","a")
            outputfile.write(case_file[:case_file.rfind('\n')])
        #remove header from next files and merge into output
        else:
            outputfile.write(case_file[case_file.find('\n'):case_file.rfind('\n')])
    outputfile.close()
if __name__ == "__main__":
    main()
