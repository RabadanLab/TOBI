rm(list = ls())

library(doMC)
registerDoMC(cores = 4)

##########################################################################
# Sourcing the functions
loc = "/Volumes/ifs"
# loc = "/Volumes/ar3177"
file_path = paste(loc, "/bin/TOBI/machine_learning/pp_scripts", sep = "")
file.sources = list.files(path = file_path,
                          pattern="*.R",
                          full.names = TRUE)
for (item in file.sources) source(item)
rm(file.sources, item, file_path)

##########################################################################
# Reading mutations
file_name = "all_GBM_mutations_Oct-3-2014_filt_with_indel_techn_biol.txt"
# file_name = "all_GBM_mutations_104cases_filt_indel_techn_biol_pre_proc.txt"
mt_file = paste(loc, 
                "/scratch/Results/GBM/final_tables/", 
                file_name,
                sep = "")
mt = my_read_table(mt_file)

################################################################################
# Adding somatic annotation
somatic_file = paste(loc, 
                     "/scratch/Results/GBM/Somatic_Mutations/CBio/cbio_somatic_mutations_XY.txt", 
                     sep = "")
somatic = my_read_table(somatic_file)
a = somatic_annot(mt, somatic)
mt = a[[1]]
somatic = a[[2]]
rm(a)

table(mt$Y)
table(somatic$Not_Present)

################################################################################
# addin gene_freq to the data
gene_freq_file = paste(loc, 
                       "/scratch/meganormal/gene_count_acl_length_merged.txt", 
                       sep = "")
gene_freq_length = my_read_table(gene_freq_file)

mt = gene_freq_fun(mt, gene_freq_length)

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
  
  output_file = paste(loc, "/scratch/Results/Finding_lost_mutations/missed_driver_mutations_lib.txt", sep = "")
  my_write_table(missed_driver, output_file)
}


################################################################################
# adding total number of cosmic sample to the data
tot_cosm_file = paste(loc, 
                      "/scratch/cosmic/tot_num_samples_gene.txt", 
                      sep = "")
tot_cosm_samp = my_read_table(tot_cosm_file)
mt = tot_cosm_fun(mt, tot_cosm_samp)

################################################################################
# Calculating mutation allele frequency
mt$freq = (mt$dp4_3 + mt$dp4_4) / 
    (mt$dp4_1 + mt$dp4_2 + mt$dp4_3 + mt$dp4_4) * 100

# calculating the pvalue for being somatic based on dp and freq
var_dp = round(mt$dp*mt$freq/100)
mt$pval_som = dbinom(var_dp, mt$dp, 0.5)

# making id2
mt$id2 = id2(mt)

# Correcting factors
time1 = Sys.time()
mt = cleaning_features(mt)
Sys.time() - time1

# Removing empty columns
mt = removing_features(mt)


################################################################################
# Add mega, and filter indels
INDEL = FALSE
if (INDEL) {
  mega_file = paste(loc,
                    "/scratch/meganormal/219normals_chr_intervals.txt",
                    sep = "")
  mega = my_read_table(mega_file)
  mt = indel_filter(mt, mega)
} else {
  mt = mt[ mt$indel == ".", ]
}

table(mt$indel)

################################################################################
# Writing the output
output_file = paste(loc, "/scratch/Results/GBM/final_tables/all_GBM_mutations_104cases_filt_with_indel_techn_biol_pre_proc.txt", sep = "")
my_write_table(mt, output_file)

