import argparse
import os
import subprocess
import cmd_gen

def get_arg():
    """Get Arguments"""

    prog_description = """TOBIv1.2 ADD DESC HERE"""
    parser = argparse.ArgumentParser(description=prog_description)


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
    
    parser.add_argument(
        '--inputdir',
        help = "directory for vcf files"
        )
    parser.add_argument(
        '--output',
        help = "output directory"
        )
    parser.add_argument(
        '--plist',
        help = "list of patient IDs"
        )
    parser.add_argument(
        '--config',
        help = "optional config file specifying command line arguments."                
        )
    parser.add_argument(
        '--steps ',
        type = str,
        help = """Specify which steps of pipeline to run. (REQUIRED)
            B: NOT IMPLEMENTED
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
        '--vcfsuffix',
        help = "suffix of vcf files that follows TCGA name"
        )
    parser.add_argument(
        '--debug',
        default = False,
        action = 'store_true',
        help = "Debug flag. NOT IMPLEMENTED"                
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

    # print args
    if args.debug == True:
        print(args)
        print
        
    return args

def main():
    args = get_arg()
    input_filenames = []
    #get list of filenames in input directory
    for (dirpath, dirnames, filenames) in os.walk(args.inputdir):
        input_filenames.extend(filenames)
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
    annotate(input_filenames, args)

def vcf_call(input_filenames, args):
    for case_name in input_filenames:
        pileup_cmd = cmd_gen.mpileup(
            args.start,
            args.end,
            args.ref,
            args.inputdir + "/" + case_name,
            args.output
            )
        proc = subprocess.Popen(
        pileup_cmd, 
        shell=True,
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE
        ) 
        #wait for job completion 
        proc.wait()
    return 

def annotate(input_filenames, args):
    for case_name in input_filenames:      
        if args.cluster == "hpc":
            cmd = "qsub -V -N filt-" + case_name \
                + "-l mem=10G,time=2:: -pe smp 2 " \
                + "-o " + args.output + "/" + case_name + "/logs " \
                + "-e " + args.output + "/" + case_name + "/logs " \
                + "-m as " \
                + "do_annotation_filtering.sh_TCGA_protected.sh " \
                + "--input-file " + args.inputdir + "/" + case_name \
                + " -s " + args.steps \
                + " --memory 6 --filter on " \
                + "--config_file " + args.config \
                + " --debug " + args.debug
            print(cmd)
            proc = subprocess.Popen(
                cmd,
                shell=True,
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE
                ) 
                #wait for job completion 
            proc.wait()