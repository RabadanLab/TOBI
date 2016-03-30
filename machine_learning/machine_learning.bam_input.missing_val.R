#!/usr/bin/env Rscript

# performs machine learning w/ gradient boosting, calculates performance statistics, saves top-predicted genes
rm(list = ls())

library(plyr); library(dplyr)
library(caret)
library(doMC)
library("session", lib.loc="~/Library/R/3.1/library")

registerDoMC(cores = 4)

args <- commandArgs(TRUE)
input_file <- args[[1]] #the pre-processed, aggregated mutation file
out_dir <- args[[2]] #directory where you want files to be saved
somatic_path <- args[[3]]
suffix <- args[[4]] #a label specific to this particular run (e.g. <date>_<disease>)
train_size <- args[[5]] # number of patients you want in the training set
TOBI_path <- args[[6]]
print(args)

################################################################################
# Sourcing the functions
#contains new na_minus1.R script that does mean and na_to_0 for converting non-numeric columns into 0
file_path = paste(TOBI_path, "/machine_learning/ml_scripts/", sep = "")
file.sources = list.files(path = file_path,
                          pattern="*.R",
                          full.names = TRUE)
for (item in file.sources) source(item)
rm(file.sources, item, file_path)
modelInfo = custom_gbm()

################################################################################
# Reading in the data and dividing it to training and testing
file_path = input_file
p = train_size #20 # number of training cases
param = c("Y", "case", "chrom", "pos", "id", "ref", "alt", "gene_name", "dp", "effect",
          "amino_acid_change", "warnings", "id2", "qual", "freq", "amino_acid_length",
          "gene_freq", "gene_length", "tot_cosm_samp", "num_mut_gene", "case_mut_gene",
          "mq", "effect_impact", "cosmic_nsamp", "pval_som", "dbNSFP_1000Gp1",
          "dbNSFP_CADD_phred", "dbNSFP_MutationAssessor_score",
          "dbNSFP_MutationTaster_score", "dbNSFP_FATHMM_score",
          "dbNSFP_Polyphen2_HDIV_score", "dbNSFP_Polyphen2_HVAR_score",
          "dbNSFP_RadialSVM_score", "dbNSFP_SIFT_score")
a = create_training_testing(file_path, p, param, 
                            dbSNP_only = TRUE, indel = FALSE)
training = a[[1]]
testing = a[[2]]
rm(a, p)

sort(unique(training$case))
dim(training)
sort(unique(testing$case))
dim(testing)

################################################################################
# Formula
my_formula = Y ~ freq + pval_som +
  gene_freq + amino_acid_length +
  num_mut_gene + case_mut_gene + #id2 +
  cosmic_nsamp + tot_cosm_samp +
  #dbNSFP_MutationTaster_score + dbNSFP_Polyphen2_HDIV_score + dbNSFP_SIFT_score +
  effect_impact + dbNSFP_CADD_phred  #dbNSFP_MutationAssessor_score
print(my_formula)

####### need to change file paths here, will eventually put this into a loop of 4 rounds so it runs in bg
perf_round <- list() #will score the performance tables here, than simply select the highest performer for machine learning 
fscor_vec <- vector()
for (i in 1:4){ 
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
                TRUE)
  Sys.time() - time1

  #testing data: prediction, performance stat, calculate top_pred_genes, save predictions
  testing$pred <- predict(modboost, testing[,-c(1:11)])
  #print(confusionMatrix(table(testing$pred, testing$Y, dnn=list('predicted','actual'))))
  print(perf_round[[i]] <- myPerf_y(testing)); fscor_vec[i] <- perf_round[[i]][3] ##storing+printing myPerf; saving fscore
  mydata <- testing
  print(pred_genes <- paste0(out_dir, "top_pred_genes.testing.", suffix, ".round",i,".txt"))
  driver_mutations_gene(mydata, somatic_path, pred_genes, output_plot = F)
  
  print(test_file <- paste0(out_dir, "cleaned.filt_indel_techn_biol.pre_proc.pred.testing.", suffix, ".round",i,".txt"))
  write.table(testing, file = test_file,
              sep = "\t", quote = FALSE, 
              na = ".", row.names = FALSE)
  
  #training data: prediction, performance stat, calculate top_pred_genes, save predictions
  training$pred <- predict(modboost, training[,-c(1:11)])
  print(confusionMatrix(table(training$pred, training$Y, dnn=list('predicted','actual'))))
  print(myPerf_y(training))
  mydata <- training
  rm(pred_genes); print(pred_genes <- paste0(out_dir, "top_pred_genes.training.", suffix, ".round",i,".txt")) 
  driver_mutations_gene(mydata, somatic_path, pred_genes, output_plot = F) #breakdown of what's filtered, top genes
  
  print(train_file <- paste0(out_dir, "cleaned.filt_indel_techn_biol.pre_proc.pred.training.", suffix, ".round",i,".txt"))
  write.table(training, file = train_file,
              sep = "\t", quote = FALSE, 
              na = ".", row.names = FALSE)
  print(param_pdf_file <- paste0(out_dir, "training.tuning_param.cases", train_size, suffix, ".",suffix, ".round",i,".pdf"))
  ggsave(param_pdf_file)

  #save session
  print(sesssion_path <- paste0(out_dir, suffix, ".ML.top_pred.round",i,".Rsession"))
  save.session(sesssion_path)
}

#print(which.max())

##########################pause !!!!! not relevant to my lgg

