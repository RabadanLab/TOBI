#!/usr/bin/env python
import argparse
import subprocess
import vcf2report as vcf2report
import parse_tsv as parse_tsv
import sys
import os

def get_arg():
#arguments for varcall
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--vcftype", 
        choices = ['default','TCGA'], 
        help="Specifies vcf type specically for TCGA filtering"
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
    case_name = args.case_name
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
        script = "/filter_indel_techn_biol.pediatric.R "
    elif args.vcftype == "TCGA":
        script ="/filter_indel_techn_biol.R "
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
        
if __name__ == "__main__":
    main()