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
        '--plist',
        help = "list of patient IDs"
        )
    parser.add_argument(
        '--vcfsuffix',
        help = "suffix of vcf files that follows TCGA name"
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
        
    input_filenames = []
    #get list of .bam filenames in input directory and remove '.bam' suffix
    for (dirpath, dirnames, filenames) in os.walk(args.inputdir):
        for filename in filenames:  
            if bool(re.search("^.*\.bam$",filename)):
                input_filenames.append(filename[:-4])
        break  
    
    for case_name in input_filenames:
        #for each case name/bam file, make output directories
        if not os.path.exists(args.output + "/" + case_name):
            os.makedirs(args.output + "/" + case_name)
            os.makedirs(args.output + "/" + case_name + "/logs")
            os.makedirs(args.output + "/" + case_name + "/output_folder") 
  
    #vcf calling
    vcf_call(input_filenames,args)

    #annotate
    #annotate(input_filenames, args)

def vcf_call(input_filenames, args):
    for case_name in input_filenames:
        source_dir = os.path.dirname(os.path.realpath(__file__))
        pileup_cmd = helpers.mpileup_cmdgen(
            args.start,
            args.end,
            args.ref,
            args.inputdir,
            case_name,
            args.output,
            source_dir
            )
        
        if args.debug:
            print('[start]' + str(args.start))
            print('[end]' + str(args.end))
            print('[ref]' + args.ref)
            print('[input]' + args.inputdir)
            print('[output]' + args.output)
            print(pileup_cmd)
        
        proc = subprocess.Popen(
            pileup_cmd, 
            shell=True,
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE
            ) 
        #wait for job completion 
        proc.wait()
        
        proc = subprocess.Popen(
            helpers.vcf_concat_cmdgen(args,case_name), 
            shell=True,
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE
            ) 
        proc.wait()
        
        helpers.purge(args.output + "/" + case_name \
            + "/output_folder", "raw_\d*\.vcf")
        
    return 

def annotate(input_filenames, args):
    for case_name in input_filenames:      
        if args.cluster == "hpc":
            cmd = helpers.annotate_cmdgen(case_name,args)
            proc = subprocess.Popen(
                cmd,
                shell=True,
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE
                ) 
                #wait for job completion 
            proc.wait()
            
if __name__ == "__main__":
    main()
