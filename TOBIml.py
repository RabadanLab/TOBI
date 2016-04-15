#!/usr/bin/env python

import argparse
import os
import sys
import subprocess

def get_arg():
    prog_description = """TOBIv1.2: 
        Tumor Only Boosting Identification of Driver Mutations.
        Machine learning step. <ADD DESC HERE>"""
    parser = argparse.ArgumentParser(description = prog_description)
    
    parser.add_argument(
        "step",
        choices = ['preprocess','machinelearning'],
        help = "pp: preprocessing step; ml: machine learning step"
        )
    parser.add_argument(
        '--input',
        help = "[REQUIRED] input file"
        )
    parser.add_argument(
        '--output',
        help = "[REQUIRED] output file"
        )
    parser.add_argument(
        '--somatic',
        help = "[REQUIRED] formatted file containing somatic variants"
        )
    parser.add_argument(
        '--log',
        help = 'Optional argument to specify a log to pipe stdout and stderr to'
        )
    parser.add_argument(
        '--check_missed',
        help = '[PP ARG] checking which mutations in important genes are missed by filtering'
        )
    parser.add_argument(
        '--suffix',
        help = '[ML ARG] a label specific to this particular run (e.g. <date>_<disease>)'
        )
    parser.add_argument(
        "--vcftype", 
        choices = ['default','TCGA'], 
        default = "default",
        help="Specifies vcf type specically for TCGA filtering"
        )
    parser.add_argument(
        '--train_size',
        help = '[ML ARG] number of patients you want in the training set.',
        )
    parser.add_argument(
        '--verbose',
        action = "store_true",
        default = False,
        help = "verbose flag"
        )
    
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    
    args = parser.parse_args()
    return args
    
def check_main_args(args):
    if args.input == None:
        sys.exit("[ERROR]: Missing required '--inputdir' argument.")
    if args.output == None:
        sys.exit("[ERROR]: Missing required '--output' argument.")
    if args.somatic == None:
        sys.exit("[ERROR]: Missing required '--somatic' argument.")
        
def check_ml_args(args):
    if args.suffix == None:
        sys.exit("[ERROR]: Missing required '--suffix' argument.")
    if args.train_size == None:
        sys.exit("[ERROR]: Missing required '--train_size' argument.")

def pp_cmdgen(args,source_dir):
    if args.log == None:
        pipe = ""
    else:
        pipe = ">> " + args.log + " 2>&1"
        
    cmd = source_dir +"/machine_learning/pre_processing.R " \
        + args.input + " " \
        + args.output + " " \
        + args.somatic + " " \
        + source_dir + " " \
        + args.vcftype + " " \
        + pipe
    if args.verbose:
        print(cmd)
    return cmd
    
def ml_cmdgen(args,source_dir):
    if args.log == None:
        pipe = ""
    else:
        pipe = ">> " + args.log + " 2>&1"
    cmd = source_dir + "/machine_learning/machine_learning.R " \
        + args.input + " " \
        + args.output + " " \
        + args.somatic + " " \
        + args.suffix + " " \
        + args.train_size + " " \
        + source_dir + " " \
        + args.vcftype + " " \
        + pipe
    if args.verbose:
        print(cmd)
    return cmd

def runShellCmd(cmd):
    proc = subprocess.Popen(
        cmd, 
        shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        ) 
    if proc.wait() != 0:
        sys.exit("[ERROR] process '" + cmd + "' terminated with non-zero exit code")
        

def main():
    args = get_arg()
    source_dir = os.path.dirname(os.path.realpath(__file__))
    check_main_args(args)
    
    #preprocessing
    if args.step == "preprocess":
        runShellCmd(pp_cmdgen(args,source_dir))
    
    #machine learning
    if args.step == "machinelearning":
        check_ml_args(args)
        runShellCmd(ml_cmdgen(args,source_dir))

if __name__ == "__main__":
    main()
    