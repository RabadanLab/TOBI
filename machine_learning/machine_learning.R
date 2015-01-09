rm(list = ls())

library(caret)
library(doMC)

registerDoMC(cores = 4)

input_file = "all_GBM_mutations_104cases_filt_indel_techn_biol_pre_proc.txt"
output_file = "cleaned_GBM_84_testing_cases_filt_techn_biol_with_pred_latest.txt"

################################################################################
# Sourcing the functions
loc = "/Volumes/ifs"
# loc = "/Volumes/ar3177"
file_path = paste(loc, "/bin/TOBI/machine_learning/ml_scripts", sep = "")
file.sources = list.files(path = file_path,
                          pattern="*.R",
                          full.names = TRUE)
for (item in file.sources) source(item)
rm(file.sources, item, file_path)
modelInfo = custom_gbm()
################################################################################
# Reading in the data and dividing it to training and testing
file_path = paste(loc, "/scratch/Results/GBM/final_tables/",
                  input_file, sep = "")
p = 20 # number of training cases
param = c("Y", "case", "chrom", "pos", "id", "ref", "alt", "gene_name", "dp", "effect",
          "amino_acid_change", "warnings", "id2", "qual", "freq", "amino_acid_length",
          "gene_freq", "gene_length", "tot_cosm_samp", "num_mut_gene", "case_mut_gene",
          "mq", "effect_impact", "cosmic_nsamp", "pval_som", "recur",
          "dbNSFP_1000Gp1", "dbNSFP_CADD_phred", "dbNSFP_MutationAssessor_score",
          "dbNSFP_MutationTaster_score", "dbNSFP_FATHMM_score",
          "dbNSFP_Polyphen2_HDIV_score", "dbNSFP_Polyphen2_HVAR_score",
          "dbNSFP_RadialSVM_score", "dbNSFP_SIFT_score")
a = create_training_testing(file_path, p, param, 
                            dbSNP_only = TRUE, indel = FALSE)
training = a[[1]]
testing = a[[2]]
rm(a, p)

################################################################################
# Formula
my_formula = Y ~ freq + pval_som + recur +
  # qual + mq + 
  # cosmic_nsamp + case_mut_gene
  gene_freq + amino_acid_length +
  num_mut_gene + id2 + 
  tot_cosm_samp + 
  dbNSFP_MutationTaster_score + dbNSFP_Polyphen2_HDIV_score + dbNSFP_SIFT_score +
  effect_impact + dbNSFP_MutationAssessor_score

# Boosting using Caret
time1 = Sys.time()
modboost = ML(training, testing, 
              modelInfo,
              #"gbm",
              custom = TRUE,
              my_formula,
              # optimization measure
              "Fscore", 
              # maximize
              TRUE, 
              #print table
              TRUE, 
              #plot cross-validation curves
              FALSE)
Sys.time() - time1

# Performance plot
# ****can only be done if the method is gbm and not modelInfo****
test_pr = predict(modboost, testing, type = "prob")$som
pr_df = performance_plot(test_pr, testing$Y, 
                         method = "Fscore", 
                         label = "gbm",
                         do_plot = TRUE)

################################################################################
# check the driver mutations
somatic_path = paste(loc, "/scratch/Results/GBM/Somatic_Mutations/CBio/cbio_somatic_SNPs_XY.txt", 
                     sep = "")

driver_mutations(testing, modboost,
                 somatic_path,
                 # output_plot
                 TRUE)

################################################################################
# checking the model with the normal data
normal_path = paste(loc, "/scratch/Results/GBM/final_tables/normal/cleaned_with_prediction_GBM_20_normals_filt_techn_biol_no_meg2.txt", sep = "")

check_normal(normal_path, modboost)

################################################################################
# Export cleaned data
mydata = testing
mydata$pred = predict(modboost, mydata[,-c(1:11)])

write.table(mydata, file = paste(loc, "/scratch/Results/GBM/final_tables/", output_file, sep = ""),
            sep = "\t", 
            quote = FALSE, 
            na = ".", 
            row.names = FALSE)
