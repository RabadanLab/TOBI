#!/usr/bin/env python
import argparse
import helpers as helpers
import sys
import os

def get_arg():
#arguments for varcall
    parser = argparse.ArgumentParser()
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
    parser.add_argument(
        '--cluster',
        choices = ['hpc','amazon'],
        help = """[REQUIRED] Specify which cluster to run on.
            hpc: run on an SGE hpc cluster
            amazon: CURRENTLY UNIMPLEMENTED"""
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
    annovcf = args.annovcf.split(',')
    case_name = args.case_name 
    dbnsfp_header="SIFT_score,SIFT_converted_rankscore,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred,LRT_score,LRT_converted_rankscore,LRT_pred,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MutationAssessor_score,MutationAssessor_rankscore,MutationAssessor_pred,FATHMM_score,FATHMM_rankscore,FATHMM_pred,RadialSVM_score,RadialSVM_rankscore,RadialSVM_pred,LR_score,LR_rankscore,LR_pred,Reliability_index,CADD_raw,CADD_raw_rankscore,CADD_phred,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,phyloP46way_primate,phyloP46way_primate_rankscore,phyloP46way_placental,phyloP46way_placental_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phastCons46way_primate,phastCons46way_primate_rankscore,phastCons46way_placental,phastCons46way_placental_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,SiPhy_29way_pi,SiPhy_29way_logOdds,SiPhy_29way_logOdds_rankscore,LRT_Omega,UniSNP_ids,1000Gp1_AC,1000Gp1_AF,1000Gp1_AFR_AC,1000Gp1_AFR_AF,1000Gp1_EUR_AC,1000Gp1_EUR_AF,1000Gp1_AMR_AC,1000Gp1_AMR_AF,1000Gp1_ASN_AC,1000Gp1_ASN_AF,ESP6500_AA_AF,ESP6500_EA_AF" 
    if args.cluster == "hpc":
        #split vcf into multiple vcfs
        if args.debug:
                print("[Splitting vcf file by chromosome]")
        for chrom in range(1,23) + ['X', 'Y', 'MT']:
            if args.debug:
                print("[Chromosome " + str(chrom) + "]")
            helpers.runShellCmd("vcftools --recode --recode-INFO-all --vcf "
                + args.inputdir + "/" + case_name +".vcf" 
                +" --out " + args.output +"/annotate/"+ case_name +"."+ str(chrom) +" --chr " + str(chrom)) 
        #snpEff annotate arrya job
        helpers.runShellCmd(helpers.snpeffarray_cmdgen(args,case_name,source_dir)) 
        #snpSIFT annotate for each vcf provided
        for vcf in annovcf:
            helpers.runShellCmd(helpers.snpsiftarray_cmdgen(args,case_name,vcf,source_dir))
            #hacky way to create & replace temp file. fix this later
            #os.remove(args.output+"/annotate/"+case_name+".eff.vcf")
            #os.rename(args.output+"/annotate/"+case_name+".eff.vcf.tmp",
            #          args.output+"/annotate/"+case_name+".eff.vcf"
            #          )
        #snpSIFT dbnsfp
        helpers.runShellCmd(
            helpers.snpdbnsfparray_cmdgen(args,case_name,args.dbnsfp,source_dir,dbnsfp_header)
            )
        #if args.cleanup:
            #os.remove(args.output+"/annotate/"+case_name+".eff.vcf")
        
        helpers.runShellCmd(helpers.vcf_snp_concat_cmdgen(args,case_name))
        
        #split effects into one effect per line
        helpers.runShellCmd(helpers.oneEff_cmdgen(args,case_name,source_dir))
        
        if args.cleanup:
            os.remove(args.output+"/annotate/"+case_name+".eff.all.vcf")
            helpers.purge(args.output+"/annotate/",case_name+".*.eff.vcf")
            
        
if __name__ == "__main__":
    main()