rm(list = ls())

##################################################################################################################
# Sourcing the functions
file.sources = list.files(path = "/Volumes/ifs/bin/TOBI/machine_learning/script",
                          pattern="*.R",
                          full.names = TRUE)
for (item in file.sources) source(item)
rm(file.sources, item, file_path)
modelInfo = custom_gbm()
##################################################################################################################
library(caret)
library(doMC)

registerDoMC(cores = 4)

# Reading in the data and dividing it to training and testing
file_path = "/Volumes/ifs/scratch/Results/GBM/final_tables/all_GBM_mutations_104cases_filt_indel_techn_biol_pre_proc.txt"
p = 20 # number of training cases
a = create_training_testing(file_path, p)
training = a[[1]]
testing = a[[2]]
rm(a, p)

##################################################################################################################
# Formula
my_formula = Y ~ freq + id2 + gene_freq + amino_acid_length + num_mut_gene + case_mut_gene + pval_som + qual + tot_cosm_samp + mq + 
  effect_impact + cosmic_nsamp + dbNSFP_MutationAssessor_score + dbNSFP_MutationTaster_score + 
  dbNSFP_Polyphen2_HDIV_score + dbNSFP_SIFT_score

# Boosting using Caret
modboost = ML(training, testing, my_formula,
              # optimization measure
              "Fscore", 
              # maximize
              TRUE, 
              #print table
              TRUE, 
              #plot
              FALSE) 

##################################################################################################################
# check the driver mutations
driver_mutations(testing, modboost[[1]], 
                 # output_plot
                 TRUE)

##################################################################################################################
# checking the model with the normal data
normal_path = "/Volumes/ifs/scratch/Results/GBM/final_tables/normal/cleaned_with_prediction_GBM_20_normals_filt_techn_biol_no_meg2.txt"

check_normal(file_path, modboost[[1]])

##################################################################################################################
# Export cleaned data
mydata = testing
mydata$pred = predict(modboost[[1]], mydata[,-c(1:11)])

write.table(mydata, file = "/Volumes/ifs/scratch/Results/GBM/final_tables/cleaned_GBM_84_testing_cases_filt_techn_biol_with_pred.txt", 
            sep = "\t", 
            quote = FALSE, 
            na = ".", 
            row.names = FALSE)
