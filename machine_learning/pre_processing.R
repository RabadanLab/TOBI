#!/usr/bin/env Rscript
# runs pre-processing, by labelling variants as truly somatic or no
rm(list = ls())

library(doMC)
registerDoMC(12) #change to 4 if running on desktop 

##########################################################################
args = commandArgs(TRUE)
input_file = args[[1]]
output_file = args[[2]]
somatic_file = args[[3]]
TOBI_path = args[[4]]
vcftype = args[[5]]
#add check_missed optional flag
#add verbose flag
print(args)
time1 = Sys.time()
##########################################################################
cat("[Sourcing helper file] \n")
file_path = paste(TOBI_path, "/machine_learning/pp_scripts/pp_helpers.R", sep = "")
source(file_path)

##########################################################################
cat("[Reading mutations] \n")
mt = my_read_table(input_file)
if(vcftype == "default"){
  mt$chrom = sapply(mt$chrom, function(x) gsub("chr", "", x))
}


################################################################################
cat("[Adding somatic annotation] \n")
somatic = my_read_table(somatic_file)
a = somatic_annot(mt, somatic)
mt = a[[1]]
somatic = a[[2]]
rm(a)

#MORE DESCRIPTIVE OUTPUT?
table(mt$Y)
table(somatic$Not_Present)

################################################################################
cat("[Adding in gene_freq to the data] \n")
gene_freq_file = paste(TOBI_path, 
  "/machine_learning/pp_data/gene_count_acl_length_merged.txt", 
  sep = "")
gene_freq_length = my_read_table(gene_freq_file)

mt = gene_freq_fun(mt, gene_freq_length)

#MORE DESCRIPTIVE OUTPUT?
head(mt$gene_freq)

################################################################################
# optional
check_missed = FALSE
if (check_missed) {
  cat("[Checking which mutations in important genes are missed by filtering] \n")
  list_driver_genes = c("PIK3R1", "PTEN", "PIK3CA", "TP53", "EGFR", "IDH1", 
                        "BRAF", "RB1", "NF1", "PDGFRA", "LRP2", "PPP2R3A", 
                        "FAT2", "GPR116", "RIMS2", "NOS1", "PIK3C2B", "MDM4",
                        "MYCN", "KIT", "DDIT3", "TSHZ2", "LRP1B", "CTNND2", 
                        "HCN1", "PKHD1", "TEK", "PCNX", "HERC2", "LZTR1", 
                        "BCOR", "ATRX", "PCDH11X", "CDK4", "CDKN2A")
  
  # writing the mutations to be checked
  missed_driver = somatic[somatic$Not_Present == 1 
                          & somatic$gene_symbol %in% missed_driver_genes, ]
  missed_driver = missed_driver[, c("gene_symbol", "case_id", "chr", "start_position")]
  
  output_file = paste(TOBI_path, "/scratch/Results/Finding_lost_mutations/missed_driver_mutations_lib.txt", sep = "")
  my_write_table(missed_driver, output_file)
}


################################################################################
cat("[Adding total number of cosmic sample to the data] \n")
tot_cosm_file = paste(TOBI_path, 
                      "/machine_learning/pp_data/tot_num_samples_gene.txt", 
                      sep = "")
tot_cosm_samp = my_read_table(tot_cosm_file)
mt = tot_cosm_fun(mt, tot_cosm_samp)

#MORE DESCRIPTIVE OUTPUT?
head(mt$tot_cosm_samp)

################################################################################
cat("[Calculating mutation allele frequency] \n")

if(vcftype == "default"){
  mt$freq = (mt$dp4_3 + mt$dp4_4) / 
    (mt$dp4_1 + mt$dp4_2 + mt$dp4_3 + mt$dp4_4) * 100
  var_dp = round(mt$dp*mt$freq/100)
  mt$pval_som = dbinom(var_dp, mt$dp, 0.5)
}
if(vcftype == "TCGA"){
  mt$freq = 100 * mt$fa_2
  var_dp = round(mt$dp_2*mt$freq/100)
  mt$pval_som = dbinom(var_dp, mt$dp_2, 0.5)
  
}
# calculating the pvalue for being somatic based on dp and freq
# making id2
mt$id2 = id2(mt)

# Correcting factors
time1 = Sys.time()
mt = cleaning_features(mt)
Sys.time() - time1

# Removing empty columns
mt = removing_features(mt)

# Removing freq lower than 1
mt = mt[ mt$freq > 1, ]

# Removing mq lower than 40
#remove for TCGA
if(vcftype == "default"){
  mt = mt[ mt$mq >= 40, ]
}

################################################################################
cat("[Add mega, and filter indels] \n")
INDEL = FALSE
if (INDEL) {
  mega_file = paste(TOBI_path,
                    "/machine_learning/pp_data/219normals_chr_intervals.txt",
                    sep = "")
  mega = my_read_table(mega_file)
  mt = indel_filter(mt, mega)
} else {
  mt = mt[ mt$indel == ".", ]
}

#MORE DESCRIPTIVE OUTPUT?
table(mt$indel)

################################################################################
cat("[Add domain-specific cosmic mutations] \n")
pfam = FALSE
if (pfam) {
  pfam_cosmic_file = paste(TOBI_path, 
                           "/machine_learning/pp_data/pfam_cosmic.txt", 
                           sep = "")
  pfam_cosmic = my_read_table(pfam_cosmic_file)
  mt = pfam_cosmic_fun(mt, pfam_cosmic)
}

################################################################################
cat("[Writing the output] \n")

#fix for new versions of cosmic
if(!("cosmic_nsamp" %in% colnames(mt))){
  colnames(mt)[which(colnames(mt) == "cnt")] = "cosmic_nsamp"
}
my_write_table(mt, output_file)
print(Sys.time() - time1)
