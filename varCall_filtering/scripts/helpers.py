import ConfigParser

def mpileup_cmdgen(start,end,ref,inputbam,output):
    cmd = "qsub -sync y -t " + start + "-" + end \
        + " -V " \
        + " -e " + output + '/' + inputbam + "/logs " \
        + " -o " + output + '/' + inputbam + "/logs " \
        + "-cwd -l mem=10G,time=1:: " \
        + "parallel_pileup.sh" \
        + " --bam " + inputbam \
        + " --ref " + ref \
        + " --outputdir " + output + '/' + inputbam + "/output_folder "
    return cmd

def annotate_cmdgen(case_name,args):
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
    if(args.debug):
        print(cmd)
    return cmd

def parse_config(args):
    Config = ConfigParser.ConfigParser()
    Config.read(args.config)
    for i in [(args.inputdir, 'inputdir', 'main'), 
              (args.output, 'output', 'main'), 
              (args.steps, 'steps', 'main'), 
              (args.cluster, 'cluster', 'main'),
              (args.ref, 'ref', 'varcall'),
              (args.start, 'start', 'varcall'),
              (args.end, 'end', 'varcall'),
              (args.plist, 'plist', 'annotate'),
              (args.vcfsuffix, 'vcfsuffix', 'annotate')
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