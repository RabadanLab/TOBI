rm(list = ls())

library(ggplot2)
library(doMC)
library(foreach)

registerDoMC(cores = 4)

GBM_Recur_Mutations = c("PIK3R1", "PTEN", "PIK3CA", "TP53", "EGFR", "IDH1", "BRAF", "RB1", "NF1", "PDGFRA", "LRP2", "PPP2R3A", "FAT2",
                        "GPR116", "RIMS2", "NOS1", "PIK3C2B", "MDM4", "MYCN", "KIT", "DDIT3", "TSHZ2", "LRP1B", "CTNND2", "HCN1", "PKHD1",
                        "TEK", "PCNX", "HERC2", "LZTR1", "BCOR", "ATRX", "PCDH11X", "CDK4", "CDKN2A")

mt = read.table("/Volumes/ifs/scratch/Results/GBM/final_tables/cleaned_GBM_84_testing_cases_filt_techn_biol_with_pred.txt", 
                   sep = "\t", 
                   header = TRUE, 
                   fill = TRUE, 
                   stringsAsFactors = FALSE)

mt$loc = paste(mt$chrom, mt$pos, sep = ":")

# to parallelize the for loop, a matrix x is filled, and then transformed into a data frame
time1 = Sys.time()
x = foreach(i=1:dim(genes)[1], .combine='rbind') %dopar% {
  y = rep(NA, 20)
  
  y[1] = mt$num_mut_gene[mt$gene_name == genes$name[i]][1]
  y[2] = max(mt$case_mut_gene[mt$gene_name == genes$name[i]])
  y[3] = max(table(mt$loc[mt$gene_name == genes$name[i]]))
  y[4] = length(unique(mt$case[mt$gene_name == genes$name[i]]))
  y[5] = mt$tot_cosm_samp[mt$gene_name == genes$name[i]][1]
  y[6] = sum(mt$effect == "NON_SYNONYMOUS_CODING" & mt$gene_name == genes$name[i]) / sum(mt$gene_name == genes$name[i]) * 100
  y[7] = sum(mt$effect == "STOP_GAINED" & mt$gene_name == genes$name[i]) / sum(mt$gene_name == genes$name[i]) * 100
  y[8] = sum(mt$effect == "SYNONYMOUS_CODING" & mt$gene_name == genes$name[i]) / sum(mt$gene_name == genes$name[i]) * 100
  y[9] = mt$gene_freq[mt$gene_name == genes$name[i]][1]
  y[10] = mt$gene_length[mt$gene_name == genes$name[i]][1]
  y[11] = mean(mt$effect_impact[mt$gene_name == genes$name[i]])
  y[12] = median(mt$dbNSFP_MutationAssessor_score[mt$gene_name == genes$name[i]])
  y[13] = median(mt$dbNSFP_MutationTaster_score[mt$gene_name == genes$name[i]])
  y[14] = sum((mt$id2 ==  3 | mt$id2 == 2) & mt$gene_name == genes$name[i]) / sum(mt$gene_name == genes$name[i]) * 100
  y[15] = sum(mt$id2 == 2 & mt$gene_name == genes$name[i]) / sum(mt$gene_name == genes$name[i]) * 100
  y[16] = sum(mt$id2 ==  1 & mt$gene_name == genes$name[i]) / sum(mt$gene_name == genes$name[i]) * 100
  y[17] = sum(mt$pred == "som" & mt$gene_name == genes$name[i]) / sum(mt$gene_name == genes$name[i]) * 100
  y[18] = sum(mt$pred == "non_som" & mt$gene_name == genes$name[i]) / sum(mt$gene_name == genes$name[i]) * 100
  y[19] = length(mt$case[mt$gene_name == genes$name[i] & mt$pred == "som"])
  y[20] = any(genes$name[i] == GBM_Recur_Mutations)
  y
}
Sys.time() - time1

genes2 = as.data.frame(unique(sort(mt$gene_name)))
genes2 = cbind(genes2, as.data.frame(x))
names(genes2) = c('name', 'tot_num_mut', 'max_num_mut_case', 'num_recur', 'num_cases', 'cosmic_nsamp', 'non_syn', 'stop_gained', 
                  'syn', 'gene_freq', 'gene_length', 'effect_impact', 'MutationAssessor', 'MutationTaster', 'per_mut_cosmic', 
                  'per_mut_dbsnp', 'per_mut_none', 'per_mut_som', 'per_mut_non_som', 'num_som_cases', 'Y')

table(genes$Y)

write.table(genes, file = "/Volumes/ifs/scratch/Results/GBM/final_tables/genes_GBM_84_testing_cases.txt", sep = "\t", quote = FALSE, na = ".", row.names = FALSE)

