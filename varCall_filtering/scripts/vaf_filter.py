#!/usr/bin/env python
import argparse
import sys
import os
import helpers as helpers

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
    
    helpers.split_vcf(args,case_name,"filter")
    #apply R script filters
    if args.vcftype == "default":
        script = "/filter_indel_techn_biol.pediatric.R "
    elif args.vcftype == "TCGA":
        script ="/filter_indel_techn_biol.R "
    else:
        sys.exit("[ERROR]: Invalid '--vcftype' argument: " + args.vcftype)
    
    helpers.runShellCmd(helpers.filterarray_cmdgen(args, case_name, source_dir, script))
    
    first = True
    for chrom in range(1,23) + ['MISC']:
        #save first file into output
        inputfile = open(args.output+"/filter/"+case_name+"."+str(chrom)+"_filt_indel_techn_biol.tsv","r")
        case_file = inputfile.read()
        inputfile.close()
        if first:
            first = False
            outputfile = open(args.output+"/filter/"+case_name+"_filt_indel_techn_biol.tsv","a")
            outputfile.write(case_file[:case_file.rfind('\n')])
        #remove header from next files and merge into output
        else:
            outputfile.write(case_file[case_file.find('\n'):case_file.rfind('\n')])
        if args.cleanup:
            os.remove(args.output+"/filter/"+case_name+"."+str(chrom)+"_filt_indel_techn_biol.tsv")
            os.remove(args.output+"/filter/"+case_name+"."+str(chrom)+".recode_not_filt.tsv")
            os.remove(args.output+"/filter/"+case_name+"."+str(chrom)+".recode.vcf")
            os.remove(args.output+"/filter/"+case_name+"."+str(chrom)+".log")
    outputfile.close()
        
if __name__ == "__main__":
    main()