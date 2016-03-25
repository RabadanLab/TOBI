def mpileup(start,end,ref,inputbam,output):
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

