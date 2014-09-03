rm(list = ls())

##################################################################################################################
# Functions
##################################################################################################################
mystep = function(x, t) {
  # x is a vector
  # t is threshold
  return(as.numeric(x >= t))
}

##################################################################################################################
nsample = function(v, m) {
  # this function randomly divides vector v into 'length(m)' groups with 'm[i]' sizes 
  # n is a vector of numbers
  # m is a vector where Sum(m) = length(v)
  y = list()
  for (i in 1:length(m)) {
    y[length(y) + 1] = list(sample(v, m[i]))
    v = v[ ! v %in% unlist(y[i]) ]
  }
  return(y)
}

##################################################################################################################
rocplot = function(x, y, method, do_plot = TRUE) {
  # x is the probablity vector for prediction of data points
  # y is the actual label of the data points
  
  n = 101
  cutoff = seq(0,1, length = n)
  roc.plot = data.frame(method = rep(method, n), 
                        cutoff = rep(NA, n), 
                        TPR = rep(NA, n), 
                        FPR = rep(NA, n), 
                        accuracy = rep(NA, n))
  
  for (i in 1:n) {
    x_bin = as.numeric(x >= cutoff[i])
    TP = sum(x_bin == 1 & y == "som")
    FN = sum(x_bin == 0 & y == "som")
    FP = sum(x_bin == 1 & y == "non_som")
    TN = sum(x_bin == 0 & y == "non_som")
    
    roc.plot$cutoff[i] = cutoff[i]
    roc.plot$TPR[i] = TP / (TP + FN)
    roc.plot$FPR[i] = 1 - TN / (TN + FP)
    roc.plot$accuracy[i] = (TP + TN) / (TP + FN + FP + TN)
  }
  
  AUC = 0
  for (i in 1:(n-1)) {
    AUC = with(roc.plot, AUC + (TPR[i] + TPR[i + 1])*(FPR[i] - FPR[i + 1])/2)
  }
  
  # ROC curve
  if (do_plot) {
    print(ggplot(data = roc.plot, aes(x = FPR, y = TPR)) + 
            geom_line(color = "red") + 
            annotate("text", x = 0.5, y = 0.5, label = paste("AUC =", round(AUC, 3))))
  } else {
    AUC
  }
  
  return(roc.plot)
}

##################################################################################################################
cleaning = function(x) {
  y = subset(x, select = c(Y, case, chrom, pos, id, ref, alt, gene_name, dp, effect, amino_acid_change, warnings,
                           id2, qual, freq, amino_acid_length, gene_freq, gene_length, cnv, tot_cosm_samp, num_mut_gene, case_mut_gene, mq,
                           effect_impact, cosmic_nsamp, pval_som, dbNSFP_1000Gp1, dbNSFP_CADD_phred,
                           dbNSFP_MutationAssessor_score, dbNSFP_MutationTaster_score,
                           dbNSFP_FATHMM_score, dbNSFP_Polyphen2_HDIV_score, dbNSFP_Polyphen2_HVAR_score,
                           dbNSFP_RadialSVM_score, dbNSFP_SIFT_score))
  
  y$Y[y$Y == 0] = "non_som"
  y$Y[y$Y == 1] = "som"
  y$Y = as.factor(y$Y)
  
  y$id2 = as.factor(y$id2)
  y$effect_impact = as.factor(y$effect_impact)
  
  y$cnv = as.numeric(y$cnv)
  y$cnv[is.na(y$cnv)] = 0
  
  y$warnings = as.numeric(y$warnings)
  y$qual = as.numeric(y$qual)
  
  y$amino_acid_length = as.numeric(y$amino_acid_length)
  y$amino_acid_length[is.na(y$amino_acid_length)] = -1
  
  y$dbNSFP_MutationAssessor_score = as.numeric(y$dbNSFP_MutationAssessor_score)
  y$dbNSFP_MutationAssessor_score[is.na(y$dbNSFP_MutationAssessor_score)] = -1
  
  y$dbNSFP_MutationTaster_score = as.numeric(y$dbNSFP_MutationTaster_score)
  y$dbNSFP_MutationTaster_score[is.na(y$dbNSFP_MutationTaster_score)] = -1
  
  y$cosmic_nsamp = as.numeric(y$cosmic_nsamp)
  y$cosmic_nsamp[is.na(y$cosmic_nsamp)] = -1
  
  y$dbNSFP_1000Gp1 = as.factor(y$dbNSFP_1000Gp1)

  y$dbNSFP_CADD_phred = as.numeric(y$dbNSFP_CADD_phred)
  y$dbNSFP_CADD_phred[is.na(y$dbNSFP_CADD_phred)] = -1
  
  y$dbNSFP_FATHMM_score = as.numeric(y$dbNSFP_FATHMM_score)
  y$dbNSFP_FATHMM_score[is.na(y$dbNSFP_FATHMM_score)] = -1
  
  y$dbNSFP_Polyphen2_HDIV_score = as.numeric(y$dbNSFP_Polyphen2_HDIV_score)
  y$dbNSFP_Polyphen2_HDIV_score[is.na(y$dbNSFP_Polyphen2_HDIV_score)] = -1
  
  y$dbNSFP_Polyphen2_HVAR_score = as.numeric(y$dbNSFP_Polyphen2_HVAR_score)
  y$dbNSFP_Polyphen2_HVAR_score[is.na(y$dbNSFP_Polyphen2_HVAR_score)] = -1
  
  y$dbNSFP_RadialSVM_score = as.numeric(y$dbNSFP_RadialSVM_score)
  y$dbNSFP_RadialSVM_score[is.na(y$dbNSFP_RadialSVM_score)] = -1
  
  y$dbNSFP_SIFT_score = as.numeric(y$dbNSFP_SIFT_score)
  y$dbNSFP_SIFT_score[is.na(y$dbNSFP_SIFT_score)] = -1

#   y$ = as.numeric(y$)
#   y$[is.na(y$)] = -1
  

  return(y)
}

################################################################################################################################################
is.letter <- function(x) grepl("[[:alpha:]]", x)

##################################################################################################################

library(caret)
# library(rpart)
# library(party)
# library(klaR)
# library(e1071)
# library(rattle)
# library(mclust)
library(doMC)
library(reshape)
# library(MASS)

registerDoMC(cores = 4)

normal = 0

mt = read.table("/Volumes/ifs/scratch/Results/GBM/final_tables/all_GBM_mutations_104cases_filt_indel_techn_biol_pre_proc.txt", 
                sep = "\t", 
                header = TRUE, 
                fill = TRUE, 
                stringsAsFactors = FALSE)

table(mt$Y)

somatic = read.table("/Volumes/ifs/scratch/Results/Somatic_Mutations/CBio/cbio_somatic_SNPs.txt", 
                     sep = "\t", 
                     header = TRUE, 
                     stringsAsFactors = FALSE)

# Choose the corresponding tumor from somatic
mt_list = unique(mt$case)
somatic = subset(somatic, somatic$case_id %in% mt_list)

mt = mt[mt$id2 != 0, ]

# num_mut_gene and case_mut_gene for mt
gene_mut_freq = as.data.frame(table(mt$gene_name))
mt$num_mut_gene = sapply(mt$gene_name, 
                               function(x) gene_mut_freq$Freq[gene_mut_freq$Var1 == x])
num_cases = length(unique(mt$case))
mt$num_mut_gene = mt$num_mut_gene / num_cases * 100

for (item in unique(mt$case)) {
  gene_mut_freq = as.data.frame(table(mt$gene_name[mt$case == item]))
  mt$case_mut_gene[mt$case == item] = sapply(mt$gene_name[mt$case == item], 
                                                         function(x) gene_mut_freq$Freq[gene_mut_freq$Var1 == x])
}
mt$case_mut_gene = mt$case_mut_gene / num_cases * 100

# cleaning up the columns
mt_subset = cleaning(mt)

##################################################################################################################
# Formula
my_formula = Y ~ freq + id2 + gene_freq + gene_length + num_mut_gene + case_mut_gene + pval_som + qual + cnv + tot_cosm_samp + mq + 
  effect_impact + cosmic_nsamp + dbNSFP_MutationAssessor_score + dbNSFP_MutationTaster_score + 
  dbNSFP_Polyphen2_HDIV_score + dbNSFP_SIFT_score

my_formula = Y ~ freq + id2 + gene_freq + gene_length + num_mut_gene + case_mut_gene + pval_som + qual + cnv + tot_cosm_samp + mq + 
  effect_impact + cosmic_nsamp

##################################################################################################################
# Boosting using Caret

fitControl = trainControl(method = "repeatedcv",
                          number = 10,
                          repeats = 10,
                          allowParallel = TRUE)

# fitControl = trainControl(allowParallel = TRUE)

gbmGrid <- expand.grid(.interaction.depth = (1:5),
                       .n.trees = (1:4)*50, 
                       .shrinkage = .1)

set.seed(2)
modboost = train(my_formula,
                 data = mt_subset,
                 trControl = fitControl,
                 method = "gbm",
                 verbose = FALSE,
                 metric = "Kappa",
                 tuneGrid = gbmGrid)
modboost

ggplot(modboost)
summary(modboost)

confusionMatrix(table(predict(modboost, mt_subset[,-c(1:11)]), mt_subset[,1], dnn=list('predicted','actual')))

##################################################################################################################
# validation

validation = read.table("/Volumes/ifs/scratch/Results/GBM/final_tables/cleaned_GBM_53cases_filt_techn_biol.txt", 
                        sep = "\t", 
                        header = TRUE, 
                        fill = TRUE, 
                        stringsAsFactors = FALSE)

validation$id2 = as.factor(validation$id2)
validation$effect_impact = as.factor(validation$effect_impact)
validation$Y = as.factor(validation$Y)

confusionMatrix(table(predict(modboost, validation[,-c(1:11)]), validation[,1], dnn=list('predicted','actual')))

##################################################################################################################
# ROC curve
mydata = validation

pred_glm = predict(modglm, mydata[,-c(1:11)], type = "response")
pred_tree = predict(modtree, mydata[,-c(1:11)], type="prob")[,2]
pred_boost = predict(modboost, mydata[,-c(1:11)], type="prob")[,2]
pred_bayes = predict(modbayes, mydata[,-c(1:11)], type="raw")[,2]
# pred_glm = predict(modglm, mydata[,-c(1:11)], type="prob")[,2]
 
act = mydata[,1]

roc_glm = rocplot(pred_glm, act, "glm")
roc_tree = rocplot(pred_tree, act, "tree")
roc_boost = rocplot(pred_boost, act, "boost")
roc_bayes = rocplot(pred_bayes, act, "bayes")

# plot all ROC curves
all_rocs = rbind(roc_glm, roc_tree, roc_boost)
ggplot(data = all_rocs, aes(x = FPR, y = TPR, color = method)) + geom_line()

# sen and spe curves
ggplot(data = roc_boost, aes(x = cutoff)) + 
  geom_line(aes(y = TPR), color = "blue") +
  geom_line(aes(y = FPR), color = "red")


##################################################################################################################
# Sanity check
mydata = mt_subset
mydata = validation

mydata$pred = predict(modboost, mydata[,-c(1:11)])

# mydata$pred = predict(modbayes, mydata[,-c(1:11)])
# mydata$pred = predict(modtree, mydata[,-c(1:11)]) 
# mydata$pred = predict(modglm, mydata[,-c(1:11)], )

# Checking which genes are predicted correctly as somatic using the CBio list
# GBM_Recur_Mutations = read.table("/Volumes/ifs/scratch/Results/GBM/GBM_Recur_Mutations.txt", sep = "\t", header = TRUE, fill = TRUE, stringsAsFactors = FALSE)

GBM_Recur_Mutations = c("PIK3R1", "PTEN", "PIK3CA", "TP53", "EGFR", "IDH1", "BRAF", "RB1", "NF1", "PDGFRA", "LRP2", "PPP2R3A", "FAT2",
                      "GPR116", "RIMS2", "NOS1", "PIK3C2B", "MDM4", "MYCN", "KIT", "DDIT3", "TSHZ2", "LRP1B", "CTNND2", "HCN1", "PKHD1",
                      "TEK", "PCNX", "HERC2", "LZTR1", "BCOR", "ATRX", "PCDH11X", "CDK4", "CDKN2A")

# (not_predicted = GBM_Recur_Mutations[! GBM_Recur_Mutations %in% mydata$gene_name[mydata$pred == "som"]])
# (not_mutated = GBM_Recur_Mutations[! GBM_Recur_Mutations %in% mydata$gene_name[mydata$Y == "som"]])
# (missed = not_predicted[! not_predicted %in% not_mutated])
# (mistaken = not_mutated[! not_mutated %in% not_predicted])

# checking which mutations in important genes are missed by Machine Learning
# missed_mt = mydata[mydata$Y == "som" & mydata$pred == "non_som", ]
# imp_missed_mt = missed_mt[missed_mt$gene_name %in% GBM_Recur_Mutations, ]

# case = "EGFR"
# gene = mydata[mydata$gene_name == case & mydata$Y == "som", ]
# View(gene[order(gene$pred),])

# correcting somatic
somatic = read.table("/Volumes/ifs/scratch/Results/Somatic_Mutations/CBio/cbio_somatic_SNPs.txt", 
                     sep = "\t", 
                     header = TRUE, 
                     stringsAsFactors = FALSE)
# Choose the corresponding tumor from somatic
mydata_list = unique(mydata$case)
somatic = subset(somatic, somatic$case_id %in% mydata_list)

l = length(GBM_Recur_Mutations)
imp_mt = data.frame(gene = GBM_Recur_Mutations, 
                    tot_num_mut = rep(0, l), 
                    num_mut_filt = rep(0, l), 
                    TP = rep(0, l), 
                    FP = rep(0, l), 
                    FN = rep(0, l))
for (i in 1:l) {
  imp_mt$tot_num_mut[i] = sum(somatic$gene_symbol == imp_mt$gene[i])
  imp_mt$num_mut_filt[i] = sum(mydata$gene_name == imp_mt$gene[i] & mydata$Y == "som")
  imp_mt$TP[i] = sum(mydata$gene_name == imp_mt$gene[i] & mydata$Y == "som" & mydata$pred == "som")
  imp_mt$FP[i] = sum(mydata$gene_name == imp_mt$gene[i] & mydata$Y == "non_som" & mydata$pred == "som")
  imp_mt$FN[i] = sum(mydata$gene_name == imp_mt$gene[i] & mydata$Y == "som" & mydata$pred == "non_som")
}
imp_mt = imp_mt[imp_mt$tot_num_mut != 0, ]
sum(imp_mt$tot_num_mut)
sum(imp_mt$num_mut_filt)
sum(imp_mt$TP)


imp_mt_melt = melt(imp_mt, id = c("gene"))

# cbPalette <- c("black", "red")
q = qplot(x = gene, y = value, fill=variable,
      data = imp_mt_melt, geom="bar", stat="identity",
      position="dodge") 

q + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          text = element_text(size = 18),
          legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="bottom") + coord_fixed(ratio = 0.6)

##################################################################################################################
# checking the model with the normal data
if (normal == 1) {
  mt_normal = read.table("/Volumes/ifs/scratch/Results/GBM/final_tables/normal/cleaned_with_prediction_GBM_20_normals_filt_techn_biol_no_meg2.txt", 
                         sep = "\t", header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
  mt_normal = cleaning(mt_normal)
  mt_normal$cnv = 1
  
  mt_normal$pred = predict(modboost, mt_normal[,-c(1:11)])
  table(predict(modboost, mt_normal[,-c(1:11)]), mt_normal[,1], dnn=list('predicted','actual'))
}

##################################################################################################################
# Export
write.table(mt_subset, file = "/Volumes/ifs/scratch/Results/GBM/final_tables/cleaned_GBM_104cases_filt_techn_biol.txt", 
            sep = "\t", 
            quote = FALSE, 
            na = ".", 
            row.names = FALSE)

