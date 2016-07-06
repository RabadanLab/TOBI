#!/usr/bin/env python

import os
import argparse
import helpers as helpers
import sys

def get_arg():
#arguments for varcall
    parser = argparse.ArgumentParser()
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
    parser.add_argument(
        '--case_name',
        help = "case name"
        )
    parser.add_argument(
        '--debug',
        default = False,
        action = 'store_true',
        help = "Debug/verbose flag. Default: False"                
        )
    parser.add_argument(
        '--inputdir',
        help = "input"
        )
    parser.add_argument(
        '--output',
        help = "output"
        )
    parser.add_argument(
        '--cleanup',
        default = True,
        action = 'store_false',
        help = "Delete temporary debug files. Default True"                
        )
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
        
    args = parser.parse_args()
    # add key-value pairs to the args dict
    vars(args)['cwd'] = os.getcwd()

    return args

def main():
    args = get_arg()
    source_dir = os.path.dirname(os.path.realpath(__file__))
    case_name = args.case_name
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
        
if __name__ == "__main__":
    main()