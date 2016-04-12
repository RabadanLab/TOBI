import ConfigParser
import os
import re
import sys

def mpileup_cmdgen(args,case_name,source_dir):
    #generate command for mpileup
    cmd = "qsub -sync y -t " + str(args.start) + "-" + str(args.end) \
        + " -V " \
        + "-N " + case_name \
        + " -e " + args.output + "/vcfcall/logs/"+ case_name +".vcfcall.e " \
        + "-o " + args.output + "/vcfcall/logs/"+ case_name +".vcfcall.o " \
        + "-cwd -l mem=10G,time=1:: " \
        + source_dir + "/parallel_pileup.sh" \
        + " --bam " + args.inputdir + "/" + case_name + ".bam"\
        + " --ref " + args.ref \
        + " --outputdir " + args.output + "/vcfcall"
    if(args.debug):
        print('[Performing mpileup]')
        print(cmd)
    return cmd

def vcf_concat_cmdgen(args,case_name):
    #generate command for vcf-concat
    vcflist = []
    for i in range(args.start,args.end+1):
        vcfname = args.output + "/vcfcall/raw_" + str(i) + ".vcf"
        vcflist.append(vcfname)
    vcfstr = " ".join(vcflist)
    cmd = "vcf-concat "+ vcfstr + " > " \
        + args.output + "/vcfcall/" + case_name + ".vcf"
    if args.debug:
        print('[Concatenating vcf files and sorting]')
        print(cmd)
    return cmd

def snpeff_cmdgen(args,case_name):
    #generate command for snpeff 
    cmd = "qsub -V -b y -sync y -N " + case_name \
        + " -l mem=10G,time=2:: -pe smp 2 " \
        + "-e " + args.output + "/annotate/logs/"+ case_name +".snpeff.e " \
        + "-o " + args.output+"/annotate/"+case_name+".eff.vcf " \
        + "java -Xmx6G " \
        + "-jar "+ args.snpeff+"/snpEff.jar -c "+ args.snpeff+"/snpEff.config" \
        + " GRCh37.71 " \
        + "-noStats -v -lof -canon " \
        + "-no-downstream -no-intergenic -no-intron -no-upstream -no-utr " \
        + args.inputdir+"/"+case_name+".vcf"\
        + " > " + args.output + "/annotate/logs/"+ case_name +".snpeff.o"
    if(args.debug):
        print('[Annotating with snpEff]')
        print(cmd)
    return cmd

def snpsift_cmdgen(args,case_name,vcf):
    #generate command for snpsift annotate
    cmd = "qsub -V -b y -sync y -N " + case_name \
        + " -l mem=10G,time=2:: -pe smp 2 " \
        + "-e " + args.output + "/annotate/logs/"+ case_name +".snpeff.e " \
        + "-o " + args.output+"/annotate/"+case_name+".eff.vcf.tmp " \
        + "java -Xmx6G " \
        + "-jar "+ args.snpeff+"/SnpSift.jar annotate -v " \
        + vcf + " " + args.output+"/annotate/"+case_name+".eff.vcf " \
        + " > " + args.output + "/annotate/logs/"+ case_name +".snpeff.o"
    if(args.debug):
        print(cmd)
    return cmd

def snpdbnsfp_cmdgen(args,case_name,dbnsfp,header):
    #generate command for snpsift dbnsfp
    cmd = "qsub -V -b y -sync y -N " + case_name \
        + " -l mem=10G,time=2:: -pe smp 2 " \
        + "-e " + args.output + "/annotate/logs/"+ case_name +".snpeff.e " \
        + "-o " + args.output+"/annotate/"+case_name+".eff.all.vcf " \
        + "java -Xmx6G " \
        + "-jar "+ args.snpeff+"/SnpSift.jar dbnsfp " \
        + dbnsfp + " -v -f " + header + " " \
        + args.output+"/annotate/"+case_name+".eff.vcf " \
        + " > " + args.output + "/annotate/logs/"+ case_name +".snpeff.o"
    if(args.debug):
        print(cmd)
    return cmd

def oneEff_cmdgen(args,case_name,source_dir):
    cmd = "cat "+ args.output + "/annotate/" + case_name + ".eff.all.vcf | " \
        + source_dir+"/scripts/vcfEffOnePerLine.pl > " \
        + args.output +"/annotate/"+case_name+ ".all.annotations.vcf"
    if(args.debug):
        print(cmd)
    return cmd

def get_filenames(inputdir,extension):
    input_filenames = []
    #get list of .bam/.vcf filenames in input directory 
    #and remove '.bam'/'.vcf' suffix
    if extension == "bam":
        pattern = "^.*\.bam$"
    if extension == "vcf":
        pattern = "^.*\.vcf$"
    for (dirpath, dirnames, filenames) in os.walk(inputdir):
        for filename in filenames:  
            if bool(re.search(pattern,filename)):
                input_filenames.append(filename[:-4])
        break
    return input_filenames

def purge(directory, pattern):
    #function to purge files matching a pattern
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
              (args.dbnsfp, 'dbnsfp', 'annotate'),
              (args.vcftype, 'vcftype', 'filter')
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

def check_main_args(args):
    if args.inputdir == None:
        sys.exit("[ERROR]: Missing required '--inputdir' argument.")
    if args.output == None:
        sys.exit("[ERROR]: Missing required '--output' argument.")
    if args.steps == None:
        sys.exit("[ERROR]: Missing required '--steps' argument.")
    if args.cluster == None:
        sys.exit("[ERROR]: Missing rquired '--cluster' argument.")
        
def check_varcall_args(args):
    if args.ref == None:
        sys.exit("[ERROR]: Missing required '--ref' argument.")
    if args.end <= args.start:
        sys.exit("[ERROR]: '--end' argument >= '--start' argument")

def check_anno_args(args):
    if args.snpeff == None:
        sys.exit("[ERROR]: Missing required '--snpeff' argument.")
    if args.dbnsfp == None:
        sys.exit("[ERROR]: Missing required '--dbnsfp' argument.")

def check_filt_args(args):
    if args.vcftype == None:
        sys.exit("[ERROR]: Missing required '--vcftype' argument.")