#!/usr/bin/env python

import argparse
import os
import varCall_filtering.scripts.helpers as helpers
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
    helpers.runShellCmd("qsub -l mem=4G,time=2:: -sync y -b y "
                        + "-N vcf_call_" + case_name  
                        + " -o " + args.output +"/vcfcall/logs/" 
                        + " -e " + args.output +"/vcfcall/logs/ "
                        +source_dir+"/varCall_filtering/scripts/vaf_vcfcall.py "
                        + "--ref " + args.ref 
                        + " --start "+ str(args.start)
                        + " --end " + str(args.end)
                        + " --case_name " + case_name
                        + " --debug"
                        + " --inputdir " + args.inputdir
                        + " --output " + args.output)
    

def annotate(case_name, args):
    source_dir = os.path.dirname(os.path.realpath(__file__))
    args.annovcf = args.annovcf.replace("\n","")
    helpers.runShellCmd("qsub -l mem=8G,time=2:: -sync y -b y "
                        + "-N annotate_" + case_name 
                        + " -o " + args.output +"/annotate/logs/"
                        + " -e " + args.output +"/annotate/logs/ "
                        +source_dir+"/varCall_filtering/scripts/vaf_annotate.py "
                        + "--snpeff " + args.snpeff 
                        + " --dbnsfp " + args.dbnsfp
                        + " --case_name " + case_name
                        + " --debug"
                        + " --inputdir " + args.inputdir
                        + " --output " + args.output
                        + " --cluster " + args.cluster
                        + " --annovcf "+ args.annovcf)
    

def filter_vcf(case_name,args):
    source_dir = os.path.dirname(os.path.realpath(__file__))
    helpers.runShellCmd("qsub -l mem=12G,time=2:: -sync y -b y "
                        + "-N filter" + case_name  
                        + " -o " + args.output +"/filter/logs"
                        + " -e " + args.output +"/filter/logs "
                        +source_dir+"/varCall_filtering/scripts/vaf_filter.py "
                        + "--vcftype " + args.vcftype 
                        + " --case_name " + case_name
                        + " --debug"
                        + " --inputdir " + args.inputdir
                        + " --output " + args.output)
    

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
