#!/usr/bin/env Rscript
# runs pre-processing, by labelling variants as truly somatic or no
rm(list = ls())

library(doMC)
registerDoMC(12) #change to 4 if running on desktop 

##########################################################################
args <- commandArgs(TRUE)
input_file <- args[[1]]
output_file <- args[[2]]
somatic_file <- args[[3]]
TOBI_path <- args[[4]]
print(args)

time1 = Sys.time()
##########################################################################
# Sourcing the functions
file_path = paste(TOBI_path, "/machine_learning/pp_scripts/pp_helpers.R", sep = "")
source(file_path)

##########################################################################
# Reading mutations
mt = my_read_table(input_file)

################################################################################
# Adding somatic annotation
somatic = my_read_table(somatic_file)
a = somatic_annot(mt, somatic)
mt = a[[1]]
somatic = a[[2]]
rm(a)

table(mt$Y)
table(somatic$Not_Present)

################################################################################
# addin gene_freq to the data
gene_freq_file = paste(TOBI_path, 
                       "/machine_learning/pp_data/gene_count_acl_length_merged.txt", 
                       sep = "")
gene_freq_length = my_read_table(gene_freq_file)

mt = gene_freq_fun(mt, gene_freq_length)
head(mt$gene_freq)

################################################################################
# optional
# checking which mutations in important genes are missed by filtering
check_missed = FALSE
if (check_missed) {
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
# adding total number of cosmic sample to the data
tot_cosm_file = paste(TOBI_path, 
                      "/machine_learning/pp_data/tot_num_samples_gene.txt", 
                      sep = "")
tot_cosm_samp = my_read_table(tot_cosm_file)
mt = tot_cosm_fun(mt, tot_cosm_samp)
head(mt$tot_cosm_samp)

################################################################################
# Calculating mutation allele frequency
mt$freq = 100 * mt$fa_2 #"Fractions of reads (excluding MQ0 from both ref and alt) supporting each reported alternative allele, per sample"  #(mt$dp4_3 + mt$dp4_4) / 
    #(mt$dp4_1 + mt$dp4_2 + mt$dp4_3 + mt$dp4_4) * 100

# calculating the pvalue for being somatic based on dp and freq
var_dp = round(mt$dp_2*mt$freq/100)
mt$pval_som = dbinom(var_dp, mt$dp_2, 0.5)

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
#mt = mt[ mt$mq >= 40, ]  #cjm, field is absent

################################################################################
# Add mega, and filter indels
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

table(mt$indel)

################################################################################
# Add domain-specific cosmic mutations
pfam = FALSE
if (pfam) {
  pfam_cosmic_file = paste(TOBI_path, 
                           "/machine_learning/pp_data/pfam_cosmic.txt", 
                           sep = "")
  pfam_cosmic = my_read_table(pfam_cosmic_file)
  mt = pfam_cosmic_fun(mt, pfam_cosmic)
}

################################################################################
# Writing the output
my_write_table(mt, output_file)
print(Sys.time() - time1)
