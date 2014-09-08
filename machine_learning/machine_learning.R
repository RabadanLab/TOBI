rm(list = ls())

##################################################################################################################
# Functions
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
  y$Y = factor(y$Y, levels = c("som", "non_som"))
  
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
  return(y)
}

##################################################################################################################
num_mut_gene = function(x) {
  gene_mut_freq = as.data.frame(table(x$gene_name))
  y = sapply(x$gene_name, 
             function(x) gene_mut_freq$Freq[gene_mut_freq$Var1 == x])
  num_cases = length(unique(x$case))
  y / num_cases * 100
}

##################################################################################################################
case_mut_gene = function(x) {
  x$case_mut_gene = NA
  for (item in unique(x$case)) {
    gene_mut_freq = as.data.frame(table(x$gene_name[x$case == item]))
    x$case_mut_gene[x$case == item] = sapply(x$gene_name[x$case == item], 
                                             function(x) gene_mut_freq$Freq[gene_mut_freq$Var1 == x])
  }
  num_cases = length(unique(x$case))
  x$case_mut_gene / num_cases * 100
}

##################################################################################################################
myPerf = function(data, lev = NULL, model = NULL) {
  TP = sum(data[, "pred"] == "som" & data[, "obs"] == "som")
  FP = sum(data[, "pred"] == "som" & data[, "obs"] == "non_som")
  FN = sum(data[, "pred"] == "non_som" & data[, "obs"] == "som")
  sens = TP / (TP + FN)
  prec = TP / (TP + FP)
  fscore = 2*TP / (2*TP + FP + FN)
  out = c(sens, prec, fscore)
  names(out) = c("Sens", "Prec", "Fscore")
  out
}

##################################################################################################################
library(caret)
library(doMC)
library(reshape)

registerDoMC(cores = 4)

normal = 0

mt = read.table("/Volumes/ifs/scratch/Results/GBM/final_tables/all_GBM_mutations_104cases_filt_indel_techn_biol_pre_proc.txt", 
                sep = "\t", 
                header = TRUE, 
                fill = TRUE, 
                stringsAsFactors = FALSE)

# remove the dbSNP only
mt = mt[mt$id2 != 0, ]

table(mt$Y)

somatic = read.table("/Volumes/ifs/scratch/Results/Somatic_Mutations/CBio/cbio_somatic_SNPs.txt", 
                     sep = "\t", 
                     header = TRUE, 
                     stringsAsFactors = FALSE)

# divide the samples into training and testing sets
mt_list = unique(mt$case)

p = 20
set.seed(23)
training_list = sample(mt_list, p)
testing_list = mt_list[! mt_list %in% training_list]

training = mt[ mt$case %in% training_list, ]
testing = mt[ mt$case %in% testing_list, ]

# choose the corresponding tumor from somatic
somatic = subset(somatic, somatic$case_id %in% mt_list)

# num_mut_gene and case_mut_gene
training$num_mut_gene = num_mut_gene(training)
training$case_mut_gene = case_mut_gene(training)

testing$num_mut_gene = num_mut_gene(testing)
testing$case_mut_gene = case_mut_gene(testing)


# cleaning up the columns
training = cleaning(training)
testing = cleaning(testing)

##################################################################################################################
# Custom gbm model
modelInfo <- list(label = "Stochastic Gradient Boosting with Threshold Selection",
                  library = c("gbm", "plyr"),
                  type = c("Classification"),
                  parameters = data.frame(parameter = c('n.trees', 'interaction.depth', 
                                                        'shrinkage', 'threshold'),
                                          class = rep("numeric", 4),
                                          label = c('# Boosting Iterations', 'Max Tree Depth', 
                                                    'Shrinkage', "Probability Cutoff")),
                  grid = function(x, y, len = NULL) 
                    expand.grid(interaction.depth = seq(1, len),
                                n.trees = floor((1:len) * 50),
                                shrinkage = .1,
                                threshold = seq(0.1, 0.99, length.out = 5)),
                  loop = function(grid) {
                    ## Create a set of unique depth and shrinkage values
                    ## that will be used as a loop
                    loop <- ddply(grid,  c('interaction.depth', 'shrinkage'), 
                                  function(x) c(n.trees = max(x$n.trees),
                                                threshold = max(x$threshold)))
                    
                    ## Now, within each row of `loop`, determine the associated
                    ## submodels that are not the largest values for the number
                    ## of iterations and the prob threshold
                    submodels <- vector(mode = "list", length = nrow(loop))
                    
                    for(i in seq(along = submodels)) {
                      submod <- subset(grid, 
                                       interaction.depth == loop$interaction.depth[i] & 
                                         shrinkage == loop$shrinkage[i] & 
                                         (n.trees < loop$n.trees[i] | 
                                            threshold < loop$threshold[i]))
                      submod <- submod[order(-submod$n.trees, -submod$threshold),]
                      submodels[[i]] <- submod[, c("n.trees", "threshold")]
                    } 
                    list(loop = loop, submodels = submodels)           
                  },
                  fit = function(x, y, wts, param, lev, last, classProbs, ...) { 
                    ## train will figure out whether we are doing classification or reggression
                    ## from the class of the outcome and automatically specify the value of
                    ## 'distribution' in the control file. If the user wants to over-ride this,
                    ## this next bit will allow this.
                    theDots <- list(...)
                    
                    ## check to see if weights were passed in (and availible)
                    if(!is.null(wts)) theDots$w <- wts     
                    if(is.factor(y) && length(lev) == 2) y <- ifelse(y == lev[1], 1, 0)
                    
                    modArgs <- list(x = x,
                                    y = y,
                                    interaction.depth = param$interaction.depth,
                                    n.trees = param$n.trees,
                                    shrinkage = param$shrinkage, 
                                    distribution = "bernoulli")
                    
                    if(length(theDots) > 0) modArgs <- c(modArgs, theDots)
                    
                    do.call("gbm.fit", modArgs)
                  },
                  predict = function(modelFit, newdata, submodels = NULL) {
                    out <- predict(modelFit, newdata, type = "response",
                                   n.trees = modelFit$tuneValue$n.trees)
                    out[is.nan(out)] <- NA
                    out <- ifelse(out >= modelFit$tuneValue$threshold, 
                                  modelFit$obsLevels[1], 
                                  modelFit$obsLevels[2])
                    if(!is.null(submodels)) {
                      tmp <- vector(mode = "list", length = nrow(submodels) + 1)
                      tmp[[1]] <- out
                      for(index in 1:nrow(submodels)) {
                        tmp_pred <- predict(modelFit, newdata, type = "response", n.trees = submodels$n.trees[index])
                        tmp[[index+1]] <- ifelse(tmp_pred >= submodels$threshold[index], 
                                                 modelFit$obsLevels[1], 
                                                 modelFit$obsLevels[2])
                      }
                      out <- tmp
                    }
                    out  
                  },
                  prob = function(modelFit, newdata, submodels = NULL) {
                    out <- predict(modelFit, newdata, type = "response",
                                   n.trees = modelFit$tuneValue$n.trees)
                    out[is.nan(out)] <- NA
                    
                    out <- ifelse(out >= modelFit$tuneValue$threshold, 
                                  modelFit$obsLevels[1], 
                                  modelFit$obsLevels[2])     
                    if(!is.null(submodels)) {
                      tmp <- vector(mode = "list", length = nrow(submodels) + 1)
                      tmp[[1]] <- out
                      for(index in 1:nrow(submodels)) {
                        tmp_pred <- predict(modelFit, newdata, type = "response", n.trees = submodels$n.trees[index])
                        tmp_pred <- cbind(tmp_pred, 1 - tmp_pred)
                        colnames(tmp_pred) <- modelFit$obsLevels
                        tmp[[index+1]] <- tmp_pred
                      }
                      out <- tmp
                    }                    
                    out  
                  },
                  predictors = function(x, ...) {
                    vi <- relative.influence(x, n.trees = x$tuneValue$n.trees)
                    names(vi)[vi > 0]
                  },
                  varImp = function(object, numTrees = NULL, ...) {
                    if(is.null(numTrees)) numTrees <- object$tuneValue$n.trees
                    varImp <- relative.influence(object, n.trees = numTrees)
                    out <- data.frame(varImp)
                    colnames(out) <- "Overall"
                    rownames(out) <- object$var.names
                    out   
                  },
                  levels = function(x) {
                    if(x$distribution$name %in% c("gaussian", "laplace", "tdist")) 
                      return(NULL)
                    if(is.null(x$classes)) {
                      out <- if(any(names(x) == "obsLevels")) x$obsLevels else NULL
                    } else {
                      out <- x$classes
                    }
                    out
                  },
                  tags = c("Tree-Based Model", "Boosting", "Ensemble Model", "Implicit Feature Selection"),
                  sort = function(x) {
                    # This is a toss-up, but the # trees probably adds
                    # complexity faster than number of splits
                    x[order(x$n.trees, x$interaction.depth, x$shrinkage),] 
                  })

##################################################################################################################
# Formula
my_formula = Y ~ freq + id2 + gene_freq + amino_acid_length + num_mut_gene + case_mut_gene + pval_som + qual + tot_cosm_samp + mq + 
  effect_impact + cosmic_nsamp + dbNSFP_MutationAssessor_score + dbNSFP_MutationTaster_score + 
  dbNSFP_Polyphen2_HDIV_score + dbNSFP_SIFT_score

##################################################################################################################
# Boosting using Caret
fitControl = trainControl(method = "repeatedcv",
                          number = 5,
                          repeats = 5,
                          allowParallel = TRUE,
                          summaryFunction = myPerf,
                          classProbs = TRUE)

gbmGrid = expand.grid(n.trees = (2:4)*50,
                      interaction.depth = 3:7,
                      shrinkage = 0.1,
                      threshold = seq(0.1, 0.99, length.out = 25))

# for all training samples
time1 = Sys.time()
set.seed(123)
modboost = train(my_formula,
                     data = training,
                     trControl = fitControl,
                     method = modelInfo,
                     verbose = FALSE,
                     metric = "Fscore",
                     maximize = TRUE,
                     tuneGrid = gbmGrid)
cat(Sys.time() - time1)

# modboost
# ggplot(modboost)
# summary(modboost)

print(confusionMatrix(table(predict(modboost, training[,-c(1:11)]), training[,1], dnn=list('predicted','actual'))))
print(confusionMatrix(table(predict(modboost, testing[,-c(1:11)]), testing[,1], dnn=list('predicted','actual'))))

##################################################################################################################
# Sanity check
mydata = testing

mydata$pred = predict(modboost, mydata[,-c(1:11)])

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
# Export cleaned data

mydata = testing

mydata$pred = predict(modboost, mydata[,-c(1:11)])


write.table(mydata, file = "/Volumes/ifs/scratch/Results/GBM/final_tables/cleaned_GBM_84_testing_cases_filt_techn_biol_with_pred.txt", 
            sep = "\t", 
            quote = FALSE, 
            na = ".", 
            row.names = FALSE)




##################################################################################################################
Fscore = function (data, lev = NULL, model = NULL) 
{
  require(pROC)
  if (!all(levels(data[, "pred"]) == levels(data[, "obs"]))) 
    stop("levels of observed and predicted data do not match")
  rocObject <- try(pROC::roc(data$obs, data[, lev[1]]), silent = TRUE)
  rocAUC <- if (class(rocObject)[1] == "try-error") 
    NA
  else rocObject$auc
  sens = sensitivity(data[, "pred"], data[, "obs"],lev[1])
  spec = specificity(data[, "pred"], data[, "obs"], lev[2])
  prec = precision(data[, "pred"], data[, "obs"])
  fscore = 2 * sens * prec / (sens + prec)
  out <- c(rocAUC, 
           sens, 
           spec,
           prec,
           fscore)
  names(out) <- c("ROC", "Sens", "Spec", "Prec", "Fscore")
  out
}
####################################################
precision = function(predictions, labels) {
  TP = sum(predictions == "som" & labels == "som")
  FP = sum(predictions == "som" & labels == "non_som")
  TP / (TP + FP)
}