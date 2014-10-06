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
### pre-processing
### adding the Y label
file_name = "all_GBM_mutations_Oct-3-2014_filt_with_indel_techn_biol.txt"
# file_name = "all_GBM_mutations_104cases_filt_indel_techn_biol_pre_proc.txt"
mt_file = paste(loc, 
                "/scratch/Results/GBM/final_tables/", 
                file_name,
                sep = "")

somatic_file = paste(loc, 
                     "/scratch/Results/GBM/Somatic_Mutations/CBio/cbio_somatic_mutations_XY.txt", 
                     sep = "")
gene_freq_file = paste(loc, 
                       "/scratch/meganormal/gene_count_acl_length_merged.txt", 
                       sep = "")
tot_cosm_file = paste(loc, 
                      "/scratch/cosmic/tot_num_samples_gene.txt", 
                      sep = "")
mega_file = paste(loc,
                  "/scratch/meganormal/219normals_chr_intervals.txt",
                  sep = "")

mt = my_read_table(mt_file)
somatic = my_read_table(somatic_file)
gene_freq_length = my_read_table(gene_freq_file)
tot_cosm_samp = my_read_table(tot_cosm_file)
mega = my_read_table(mega_file)

mt_list = unique(mt$case)
n = length(mt_list)
som_list = unique(somatic$case_id)
mt_list = mt_list[mt_list %in% som_list]

# Choose the corresponding cases from somatic
somatic = subset(somatic, somatic$case_id %in% mt_list)

somatic$ID = paste(somatic$case_id, somatic$chr, somatic$start_position, sep = ":")
mt$ID = paste(mt$case, mt$chrom, mt$pos, sep = ":")

# Finding the somatic mutations using CBio data
somatic$Not_Present = 1
mt$Y = "non_som"

time1 = Sys.time()
for (id in somatic$ID) {
    a = which(mt$ID == id)
    if (length(a)) {
        mt$Y[a] = "som"
        somatic$Not_Present[somatic$ID == id] = 0
    }
}
Sys.time() - time1

table(mt$Y)
table(somatic$Not_Present)

table(mt$indel)

# optional
# checking which mutations in important genes are missed by filtering
list_driver_genes = c("PIK3R1", "PTEN", "PIK3CA", "TP53", "EGFR", "IDH1", 
                      "BRAF", "RB1", "NF1", "PDGFRA", "LRP2", "PPP2R3A", 
                      "FAT2", "GPR116", "RIMS2", "NOS1", "PIK3C2B", "MDM4",
                      "MYCN", "KIT", "DDIT3", "TSHZ2", "LRP1B", "CTNND2", 
                      "HCN1", "PKHD1", "TEK", "PCNX", "HERC2", "LZTR1", 
                      "BCOR", "ATRX", "PCDH11X", "CDK4", "CDKN2A")

missed_somatic_genes = unique(somatic$gene_symbol[somatic$Not_Present == 1])
(missed_driver_genes = missed_somatic_genes[missed_somatic_genes %in% 
                                                list_driver_genes])
View(somatic[somatic$Not_Present == 1 & somatic$gene_symbol %in% missed_driver_genes, ])


# writing the mutations to be checked
# missed_driver = somatic[somatic$Not_Present == 1 & somatic$gene_symbol %in% missed_driver_genes, ]
# missed_driver = missed_driver[, c("gene_symbol", "case_id", "chr", "start_position")]
# output_file = paste(loc, "/scratch/Results/Finding_lost_mutations/missed_driver_mutations_lib.txt", sep = "")
# my_write_table(missed_driver, output_file)

# adding total number of cosmic sample to the data
l = dim(mt)[1]
time1 = Sys.time()
mt$tot_cosm_samp = foreach (i = 1:l, .combine = c) %dopar% {
    row_num = which(tot_cosm_samp$gene == mt$gene_name[i])
    if (length(row_num) > 0) {
        tot_cosm_samp$tot_nsample[row_num]
    } else {
        0
    }
}
Sys.time() - time1

# Calculating mutation allele frequency
mt$freq = (mt$dp4_3 + mt$dp4_4) / 
    (mt$dp4_1 + mt$dp4_2 + mt$dp4_3 + mt$dp4_4) * 100

# calculating the pvalue for being somatic based on dp and freq
var_dp = round(mt$dp*mt$freq/100)
mt$pval_som = dbinom(var_dp, mt$dp, 0.5)

# making id2
mt$id2 = id2(mt)

# Correcting factors
mt = cleaning_features(mt)

# Removing empty columns
mt = removing_features(mt)

# addin gene_freq to the data
gene_freq_length$freq = with(gene_freq_length, count / amino_acid_length)
gene_freq_length$freq = sapply(gene_freq_length$freq, function(x) ifelse(is.nan(x),0,x))

time1 = Sys.time()
tmp_output = foreach (i=1:dim(mt)[1], .combine = rbind) %dopar% {
  row_num = which(gene_freq_length$gene_name == mt$gene_name[i])
  if (length(row_num)) {
    c(gene_freq_length$freq[row_num], gene_freq_length$length[row_num])
  } else {
    c(NA, NA)
  }
}
mt$gene_freq = tmp_output[,1]
mt$gene_length = tmp_output[,2]
mt = mt[! is.na(mt$gene_freq), ]
rm(tmp_output)
Sys.time() - time1

INDEL = TRUE

if (INDEL) {
  # flag indels close to meganormals
  l = dim(mt)[1]
  time1 = Sys.time()
  mt$mega = foreach (i = 1:l, .combine = c) %dopar% {
    mega_fun(mt$chrom[i], mt$pos[i], mega)
  }
  Sys.time() - time1
  
  table(mt$mega[mt$indel == 1])
  table(mt$mega[mt$indel == 1 & mt$Y == "som"])
  
  mt$ref_len = nchar(mt$ref)
  mt$alt_len = nchar(mt$alt)
  
  mt$is_1 = as.numeric(mt$is_1)
  mt$is_2 = as.numeric(mt$is_2)
  
  # only deletions
  mt = mt[ mt$indel == "." | (mt$ref_len > mt$alt_len), ]
  
  # only Frame_shift
  mt = mt[ mt$indel == "." | (mt$effect == "FRAME_SHIFT"), ]
  
  # indels longer than 10 are excluded
  mt = mt[ mt$ref_len <= 10 & mt$alt_len <= 10, ]
  
  # remvoing indels close to meganormal
  mt = mt[ mt$indel == "." | mt$mega == 0, ]
  
  # removing indels with is_1 < 100
  mt = mt[ mt$indel == "." | mt$is_1 < 100, ]
  
  table(mt$Y[mt$indel == "1"])
  
  mt$ref_len = NULL
  mt$alt_len = NULL
}

# optional
# removing the abnormals
# num_cases = length(table(mt$case))
# y = data.frame(case = names(table(mt$case)), non_som = rep(0, num_cases), som = rep(0, num_cases))
# 
# for (case in unique(mt$case)) {
#     y$non_som[y$case == case] = table(mt$Y[mt$case == case])[1]
#     y$som[y$case == case] = table(mt$Y[mt$case == case])[2]
# }
# 
# exclude_list = as.character(y$case[is.na(y$som)])
# flag = ! mt$case %in% exclude_list

output_file = paste(loc, "/scratch/Results/GBM/final_tables/all_GBM_mutations_104cases_filt_with_indel_techn_biol_pre_proc.txt", sep = "")
my_write_table(mt, output_file)

