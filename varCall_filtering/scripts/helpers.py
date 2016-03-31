import ConfigParser
import os
import re

def mpileup_cmdgen(args,case_name,source_dir):
    cmd = "qsub -sync y -t " + str(args.start) + "-" + str(args.end) \
        + " -V " \
        + "-N " + case_name \
        + " -e " + args.output + "/logs/"+ case_name +".vcfcall.e " \
        + "-o " + args.output + "/logs/"+ case_name +".vcfcall.o " \
        + "-cwd -l mem=10G,time=1:: " \
        + source_dir + "/parallel_pileup.sh" \
        + " --bam " + args.inputdir + "/" + case_name + ".bam"\
        + " --ref " + args.ref \
        + " --outputdir " + args.output 
    if(args.debug):
        print('[Performing mpileup]')
        print(cmd)
    return cmd

def snpeff_cmdgen(args,case_name):
    cmd = "qsub -V -b y -sync y -N " + case_name \
        + " -l mem=10G,time=2:: -pe smp 2 " \
        + "-e " + args.output + "/logs/"+ case_name +".snpeff.e " \
        + "-o " + args.output+"/"+case_name+".eff.vcf " \
        + "java -Xmx6G " \
        + "-jar "+ args.snpeff+"/snpEff.jar -c "+ args.snpeff+"/snpEff.config" \
        + " GRCh37.71 " \
        + "-noStats -v -lof -canon " \
        + "-no-downstream -no-intergenic -no-intron -no-upstream -no-utr " \
        + args.inputdir+"/"+case_name+".vcf.gz"\
        + " > " + args.output + "/logs/"+ case_name +".snpeff.o"
    if(args.debug):
        print('[Annotating with snpEff]')
        print(cmd)
    return cmd

def snpsift_cmdgen(args,case_name,vcf):
    cmd = "qsub -V -b y -sync y -N " + case_name \
        + " -l mem=10G,time=2:: -pe smp 2 " \
        + "-e " + args.output + "/logs/"+ case_name +".snpeff.e " \
        + "-o " + args.output+"/"+case_name+".eff.vcf.tmp " \
        + "java -Xmx6G " \
        + "-jar "+ args.snpeff+"/SnpSift.jar annotate -v " \
        + vcf + " " + args.output+"/"+case_name+".eff.vcf " \
        + " > " + args.output + "/logs/"+ case_name +".snpeff.o"
    if(args.debug):
        print(cmd)
    return cmd

def snpdbnsfp_cmdgen(args,case_name,dbnsfp,header):
    cmd = "qsub -V -b y -sync y -N " + case_name \
        + " -l mem=10G,time=2:: -pe smp 2 " \
        + "-e " + args.output + "/logs/"+ case_name +".snpeff.e " \
        + "-o " + args.output+"/"+case_name+".eff.vcf " \
        + "java -Xmx6G " \
        + "-jar "+ args.snpeff+"/SnpSift.jar dbnsfp " \
        + dbnsfp + " -v -f " + header + " " \
        + args.output+"/"+case_name+".eff.vcf " \
        + " > " + args.output + "/logs/"+ case_name +".snpeff.o"
    if(args.debug):
        print(cmd)
    return cmd

def vcf_concat_cmdgen(args,case_name):
    vcflist = []
    for i in range(args.start,args.end):
        vcfname = args.output + "/raw_" + str(i) + ".vcf"
        vcflist.append(vcfname)
    vcfstr = " ".join(vcflist)
    cmd = "vcf-concat "+ vcfstr + " |gzip -c > " \
        + args.output + "/" + case_name + ".vcf.gz"
    if args.debug:
        print('[Concatenating vcf files and sorting]')
        print(cmd)
    return cmd

def get_filenames(inputdir,extension):
    input_filenames = []
    #get list of .bam/.vcf filenames in input directory 
    #and remove '.bam'/'.vcf' suffix
    if extension == "bam":
        pattern = "^.*\.bam$"
        snip_val = -4
    if extension == "vcf":
        pattern = "^.*\.vcf.gz$"
        snip_val = -7
    for (dirpath, dirnames, filenames) in os.walk(inputdir):
        for filename in filenames:  
            if bool(re.search(pattern,filename)):
                input_filenames.append(filename[:snip_val])
        break
    return input_filenames

def purge(directory, pattern):
    for f in os.listdir(directory):
        if re.search(pattern, f):
            os.remove(os.path.join(directory, f))

def parse_config(args):
    Config = ConfigParser.ConfigParser()
    Config.read(args.config)
    for i in [(args.inputdir, 'inputdir', 'main'), 
              (args.output, 'output', 'main'), 
              (args.steps, 'steps', 'main'), 
              (args.cluster, 'cluster', 'main'),
              (args.ref, 'ref', 'varcall'),
              (args.snpeff, 'snpeff', 'annotate'),
              (args.annovcf, 'annovcf', 'annotate'),
              (args.dbnsfp, 'dbnsfp', 'annotate')
              ]:
        if not i[0] and i[1] in ConfigSectionMap(Config, i[2]):
            vars(args)[i[1]] = ConfigSectionMap(Config, i[2])[i[1]]
    return args

def ConfigSectionMap(Config, section):
    """Process config file"""
    # https://wiki.python.org/moin/ConfigParserExamples
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1