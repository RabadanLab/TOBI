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
        
        #purge raw_n.vf files
        #ADD FLAG HERE
        if True:
            helpers.purge(args.output, "raw_\d*\.vcf")
        
    return 

def annotate(input_filenames, args):
    source_dir = os.path.dirname(os.path.realpath(__file__))
    annovcf = args.annovcf.replace("\n","").split(',')
    for case_name in input_filenames:      
        if args.cluster == "hpc":
            proc = subprocess.Popen(
                helpers.snpeff_cmdgen(args,case_name),
                shell=True,
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE
                ) 
            #wait for job completion 
            proc.wait()
            #for each annotation ref given, run snpSift
            for vcf in annovcf:
                proc = subprocess.Popen(
                    helpers.snpsift_cmdgen(args,case_name,vcf),
                    shell=True,
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.PIPE
                    )
                proc.wait()
            
if __name__ == "__main__":
    main()
