case_mut_gene = function(x) {
  # this function calculates the total number of mutations in a gene for every case, 
  # and adds it to the mutation as a feature
  # x is the dataframe
  # output is a vector
  # suggested usage: x$case_mut_gene = case_mut_gene(x)
  x$case_mut_gene = NA
  for (item in unique(x$case)) {
    gene_mut_freq = as.data.frame(table(x$gene_name[x$case == item]))
    x$case_mut_gene[x$case == item] = sapply(x$gene_name[x$case == item], 
                                             function(x) gene_mut_freq$Freq[gene_mut_freq$Var1 == x])
  }
  num_cases = length(unique(x$case))
  x$case_mut_gene / num_cases * 100
}

check_normal = function(file_path, modboost, param) {
  # file_path is the normal data after pre-processing
  cat("- Reading in the normal data\n")
  mt = read.table(file_path, 
                  sep = "\t", 
                  header = TRUE, 
                  fill = TRUE, 
                  stringsAsFactors = FALSE)
  # remove the dbSNP only
  mt = mt[mt$id2 != 0, ]
  
  cat("- Creating num_mut_gene and case_mut_gene\n")
  mt$num_mut_gene = num_mut_gene(mt)
  mt$case_mut_gene = case_mut_gene(mt)
  
  cat("- Cleaning the data\n")
  mt = cleaning(mt, param)
  
  cat("- Confusion table\n")
  print(table(predict(modboost, mt[,-c(1:11)]), mt[,1], dnn=list('predicted','actual')))
}

cleaning = function(x, param, indel = FALSE) {
  # this function gets the mutation list AFTER pre-processing step, 
  # and "cleans" it for machine learning step
  # x is the dataframe
  # the output is the cleaned dataframe
  
  # add indel related features, if true
  if (indel) {
    param = c(param, "indel", "is_1", "is_2", "mega")
  }
  # test if parameters match with column names
  if (any(! param %in% colnames(x))) {
    not_found = param[! param %in% colnames(x)]
    stop(paste(not_found, collapse = " & "), " not found. Check if all 
         parameters (param) are present as columns in mutation df.")
  }
  
  y = subset(x, select = param)
  
  y$Y = as.factor(y$Y)
  y$Y = factor(y$Y, levels = c("som", "non_som"))
  
  y$id2 = as.factor(y$id2)
  y$effect_impact = as.factor(y$effect_impact)
  
  y$warnings = as.numeric(y$warnings)
  if ("qual" %in% colnames(y)) y$qual = as.numeric(y$qual)
  
  columns = c("amino_acid_length", "cosmic_nsamp")
  
  if (indel) {
    columns = c(columns, "indel", "is_1", "is_2")
    for (column in columns) y[, column] = na_to_0(y[, column])
  } else {
    for (column in columns) y[, column] = na_to_0(y[, column])
  }
  
  dbNSFP = c("dbNSFP_CADD_phred", "dbNSFP_RadialSVM_score", 
             "dbNSFP_MutationAssessor_score", "dbNSFP_MutationTaster_score", 
             "dbNSFP_Polyphen2_HDIV_score", "dbNSFP_Polyphen2_HVAR_score", 
             "dbNSFP_SIFT_score", "dbNSFP_FATHMM_score")
  for (column in dbNSFP) y[, column] = na_to_mean(y[, column]) 
  y
}

create_training_testing = function(file_path, p, param, 
                                   dbSNP_only = TRUE, indel = FALSE) {
  cat("- Reading file\n")
  y = read.table(file_path, 
                 sep = "\t", 
                 header = TRUE, 
                 fill = TRUE, 
                 stringsAsFactors = FALSE)
  # remove the dbSNP only
  if (dbSNP_only) y = y[y$id2 != 0, ]
  
  # divide the samples into training and testing sets
  y_list = unique(y$case)
  
  no_testing = FALSE
  if (length(y_list) == p) no_testing = TRUE
  
  cat("- Training\n")
  set.seed(23)
  training_list = sample(y_list, p)
  training = y[ y$case %in% training_list, ]
  training$num_mut_gene = num_mut_gene(training)
  training$case_mut_gene = case_mut_gene(training)
  training$recur = recurrent(training)
  training = cleaning(training, param, indel)
  
  if (! no_testing) {
    cat("- Testing\n")
    testing_list = y_list[! y_list %in% training_list]
    testing = y[ y$case %in% testing_list, ]
    testing$num_mut_gene = num_mut_gene(testing)
    testing$case_mut_gene = case_mut_gene(testing)
    testing$recur = recurrent(testing)
    testing = cleaning(testing, param, indel)
    
    y = list(training, testing)
  } else {
    y = list(training, FALSE)
  }
  
  return(y)
}

## Get the model code for the original bayesGlm method:

thresh_code <- getModelInfo("bayesglm", regex = FALSE)[[1]]
thresh_code$type <- c("Classification")
## Add the threshold as another tuning parameter
thresh_code$parameters <- data.frame(parameter = c("parameter", "threshold"),
                                     class = c("character", "numeric"),
                                     label = c("parameter",
                                               "Probability Cutoff"))
## The default tuning grid code:
thresh_code$grid <- function(x, y, len = NULL) {
  expand.grid(parameter = "none",
              threshold = seq(.01, .99, length = len))
}

## Here we fit a single random forest model (with a fixed mtry)
## and loop over the threshold values to get predictions from the same
## randomForest model.
thresh_code$loop = function(grid) {
  library(plyr)
  loop <- ddply(grid, c("parameter"),
                function(x) c(threshold = max(x$threshold)))
  submodels <- vector(mode = "list", length = nrow(loop))
  for(i in seq(along = loop$threshold)) {
    index <- which(grid$parameter == loop$parameter[i])
    cuts <- grid[index, "threshold"]
    submodels[[i]] <- data.frame(threshold = cuts[cuts != loop$threshold[i]])
  }
  list(loop = loop, submodels = submodels)
}

thresh_code$fit = function(x, y, wts, param, lev, last, classProbs, ...) {
  dat <- if(is.data.frame(x)) x else as.data.frame(x)
  dat$.outcome <- y
  theDots <- list(...)
  if(!any(names(theDots) == "family"))
    theDots$family <- if(is.factor(dat$.outcome)) binomial() else gaussian()              
  
  ## pass in any model weights
  if(!is.null(wts)) theDots$weights <- wts
  
  modelArgs <- c(list(formula = as.formula(".outcome ~ ."),
                      data = dat),
                 theDots)
  
  out <- do.call("bayesglm", modelArgs)
  out$call <- NULL
  out
}

## Now get a probability prediction and use different thresholds to
## get the predicted class
thresh_code$predict = function(modelFit, newdata, submodels = NULL) {
  class1Prob <- predict(modelFit,
                        newdata,
                        type = "response")[, modelFit$obsLevels[1]]
  ## Raise the threshold for class #1 and a higher level of
  ## evidence is needed to call it class 1 so it should 
  ## decrease sensitivity and increase specificity
  out <- ifelse(class1Prob >= modelFit$tuneValue$threshold,
                modelFit$obsLevels[1],
                modelFit$obsLevels[2])
  if(!is.null(submodels)) {
    tmp2 <- out
    out <- vector(mode = "list", length = length(submodels$threshold))
    out[[1]] <- tmp2
    for(i in seq(along = submodels$threshold)) {
      out[[i+1]] <- ifelse(class1Prob >= submodels$threshold[[i]],
                           modelFit$obsLevels[1],
                           modelFit$obsLevels[2])
    }
  }
  out
}

## The probabilities are always the same but we have to create
## mulitple versions of the probs to evaluate the data across
## thresholds
thresh_code$prob = function(modelFit, newdata, submodels = NULL) {
  out <- as.data.frame(predict(modelFit, newdata, type = "response"))
  if(!is.null(submodels)) {
    probs <- out
    out <- vector(mode = "list", length = length(submodels$threshold)+1)
    out <- lapply(out, function(x) probs)
  }
  out
}

custom_gbm = function() {
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
  modelInfo
}

driver_mutations = function(mydata, modboost, somatic_path, output_plot = TRUE) {
  # this function finds how many driver mutation are left after filtering and ML,
  # and plots the outcome per gene
  # mydata is a dataframe
  
  require(ggplot2)
  require(reshape2)
  
  mydata$pred = predict(modboost, mydata[,-c(1:11)])
  
  GBM_Recur_Mutations = c("PIK3R1", "PTEN", "PIK3CA", "TP53", "EGFR", "IDH1", 
                          "BRAF", "RB1", "NF1", "PDGFRA", "LRP2", "PPP2R3A", 
                          "FAT2", "GPR116", "RIMS2", "NOS1", "PIK3C2B", "MDM4", 
                          "MYCN", "KIT", "DDIT3", "TSHZ2", "LRP1B", "CTNND2", 
                          "HCN1", "PKHD1", "TEK", "PCNX", "HERC2", "LZTR1", 
                          "BCOR", "ATRX", "PCDH11X", "CDK4", "CDKN2A")
  
  GBM_Recur_Mutations = GBM_Recur_Mutations[GBM_Recur_Mutations %in% 
                                              unique(mydata$gene_name)]
  
  # somatic mutations
  somatic = read.table(somatic_path, 
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
    imp_mt$num_mut_filt[i] = sum(mydata$gene_name == imp_mt$gene[i] & 
                                   mydata$Y == "som")
    imp_mt$TP[i] = sum(mydata$gene_name == imp_mt$gene[i] & 
                         mydata$Y == "som" & 
                         mydata$pred == "som")
    imp_mt$FP[i] = sum(mydata$gene_name == imp_mt$gene[i] & 
                         mydata$Y == "non_som" & 
                         mydata$pred == "som")
    imp_mt$FN[i] = sum(mydata$gene_name == imp_mt$gene[i] & 
                         mydata$Y == "som" & 
                         mydata$pred == "non_som")
  }
  imp_mt = imp_mt[imp_mt$tot_num_mut != 0, ]
  cat("Total number of driver mutations:", sum(imp_mt$tot_num_mut), "\n")
  cat("After filter:", sum(imp_mt$num_mut_filt), "\n")
  cat("After ML:", sum(imp_mt$TP), "\n")
  
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
  }
  
  if (output_plot) {
    imp_mt_melt = melt(imp_mt, id = c("gene"))
    
    q = ggplot(imp_mt_melt, aes(gene, value, fill = variable)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = gg_color_hue(5),
                        labels = c("Before filter", "After filter", "TP", "FP", "FN")) +
      xlab("") +
      ylab("# Mutations") +
      theme_bw(base_size = 14) +
      theme(legend.position=c(0.5, 0.95), 
            legend.title=element_blank(), 
            legend.direction = "horizontal") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank()) +
      guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))
    
    print(q)
  }
}

ML = function(training, testing, model, custom, my_formula, my_metric = "Fscore", 
              max_min = TRUE, do_table = TRUE, do_plot = FALSE) {
  # this function runs the machine learning
  # training can be a list. if not, it will be converted to one
  # training should be a dataframe
  # for Fscore, use max_min=TRUE
  # for Dist, use max_min=FALSE
  
  require(caret)
  require(ggplot2)
  
  fitControl = trainControl(method = "repeatedcv",
                            number = 5,
                            repeats = 5,
                            allowParallel = TRUE,
                            summaryFunction = myPerf,
                            classProbs = TRUE)
  
  if (custom) {
    gbmGrid = expand.grid(n.trees = (2:4)*50,
                          interaction.depth = 3:7,
                          shrinkage = 0.1,
                          threshold = seq(0.05, 0.95, length.out = 25))
  } else {
    gbmGrid = expand.grid(n.trees = (2:4)*50,
                          interaction.depth = 3:7,
                          shrinkage = 0.1)
  }
  
  # check if "training" is a list, and convert it if not
  #   if (is.data.frame(training)) {
  #     a = list()
  #     a[[1]] = training
  #     training = a
  #     rm(a)
  #   }
  #   n = length(training)
  
  #   modboost = list()
  #   set.seed(123)
  #   for (i in 1:n) {
  modboost = train(my_formula,
                   data = training,
                   trControl = fitControl,
                   method = model,
                   verbose = FALSE,
                   metric = my_metric,
                   maximize = max_min,
                   tuneGrid = gbmGrid)
  if (do_plot) {
    print(ggplot(modboost) + 
            theme_bw(base_size = 16) +
            theme(legend.position="top"))
  }
  
  
  # library(reshape2)
  # metrics <- modboost$results[, c(4:7)]
  # metrics <- melt(metrics, id.vars = "threshold",
  #                 variable.name = "Resampled",
  #                 value.name = "Data")
  # 
  # ggplot(metrics, aes(x = threshold, y = Data, color = Resampled)) +
  #   geom_point() +
  #   ylab("") + xlab("Probability Cutoff") +
  #   theme(legend.position = "top") + theme_gray(base_size = 18)
  
  if (do_table) {
    cat("- Testing confusion table\n")
    if (testing == FALSE) {
      testing = training
      cat("- Testing dataset not available. Training dataset is used.\n")
    }
    print(confusionMatrix(table(predict(modboost, testing[,-c(1:11)]), testing[,1], dnn=list('predicted','actual'))))
  }
  
  return(modboost)
}

myPerf = function(data, lev = NULL, model = NULL) {
  # this functuon is used in trainControl to give different measures for optimization
  # data has two columns: "pred" and "obs"
  TP = sum(data[, "pred"] == "som" & data[, "obs"] == "som")
  TN = sum(data[, "pred"] == "non_som" & data[, "obs"] == "non_som")
  FP = sum(data[, "pred"] == "som" & data[, "obs"] == "non_som")
  FN = sum(data[, "pred"] == "non_som" & data[, "obs"] == "som")
  sens = TP / (TP + FN)
  prec = TP / (TP + FP)
  spec = TN / (TN + FP)
  fscore = 2*TP / (2*TP + FP + FN)
  dist = sqrt((1 - sens)^2 + (1 - spec)^2)
  accuracy = (TP + TN) / (TP + TN + FP + FN)
  out = c(sens, prec, fscore, dist, accuracy)
  names(out) = c("Sens", "Prec", "Fscore", "Dist", "Accuracy")
  out
}

myPerf_y = function(data, lev = NULL, model = NULL) {
  # this functuon is used in trainControl to give different measures for optimization
  TP = sum(data[, "pred"] == "som" & data[, "Y"] == "som")
  TN = sum(data[, "pred"] == "non_som" & data[, "Y"] == "non_som")
  FP = sum(data[, "pred"] == "som" & data[, "Y"] == "non_som")
  FN = sum(data[, "pred"] == "non_som" & data[, "Y"] == "som")
  sens = TP / (TP + FN)
  prec = TP / (TP + FP)
  spec = TN / (TN + FP)
  fscore = 2*TP / (2*TP + FP + FN)
  dist = sqrt((1 - sens)^2 + (1 - spec)^2)
  accuracy = (TP + TN) / (TP + TN + FP + FN)
  out = c(sens, prec, fscore, dist, accuracy)
  names(out) = c("Sens", "Prec", "Fscore", "Dist", "Accuracy")
  out
}

na_to_0 = function(x) {
  # This function converts a vector into numbers, and the resulting NAs into -1
  x = as.numeric(x)
  x[is.na(x)] = 0 
  x
}

na_to_mean = function(x) {
  # This function converts a vector into numbers, 
  # and the resulting NAs into the mean of the present values 
  x = as.numeric(x)
  y = as.numeric(x[x!=-1])
  new = mean(y, na.rm = TRUE)
  x[x==-1] = new #-1
  x
}

num_mut_gene = function(x) {
  # this function calculates the total number of mutation in a gene for all samples, 
  # and adds it to the mutation as a feature
  # x is the dataframe
  # output is a vector
  # suggested usage: x$num_mut_gene = num_mut_gene(x)
  gene_mut_freq = as.data.frame(table(x$gene_name))
  y = sapply(x$gene_name, 
             function(x) gene_mut_freq$Freq[gene_mut_freq$Var1 == x])
  num_cases = length(unique(x$case))
  y / num_cases * 100
}

performance_plot = function(x, y, method = "ROC", label = "gbm", do_plot = TRUE) {
  # x is the probablity vector for prediction of data points
  # y is the actual label of the data points
  # method is "ROC", "PR" or "Fscore"
  # label is the name of machine learning model, e.g. gbm
  require(ggplot2)
  
  n = 101
  cutoff = seq(0,1, length = n)
  roc.plot = data.frame(method = rep(label, n), 
                        cutoff = rep(NA, n), 
                        TPR = rep(NA, n), 
                        FPR = rep(NA, n), 
                        accuracy = rep(NA, n),
                        Precision = rep(NA, n),
                        Recall = rep(NA, n),
                        fscore = rep(NA, n))
  
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
    roc.plot$Precision[i] = TP / (TP + FP)
    roc.plot$Recall[i] = roc.plot$TPR[i]
    roc.plot$fscore[i] = 2*TP / (2*TP + FP + FN)
  }
  
  AUC = 0
  for (i in 1:(n-1)) {
    AUC = with(roc.plot, AUC + (TPR[i] + TPR[i + 1])*(FPR[i] - FPR[i + 1])/2)
  }
  
  # ROC curve
  if (do_plot) {
    if (method == "ROC") {
      cat("AUC is: ", AUC, "\n")
      print(ggplot(data = roc.plot, aes(x = FPR, y = TPR)) + 
              geom_line(color = "red") + 
              annotate("text", x = 0.5, y = 0.5, label = paste("AUC =", round(AUC, 3))) +
              theme_bw(base_size = 16) +
              theme(legend.position="top"))
      
    } else if (method == "PR") {
      print(ggplot(data = roc.plot, aes(x = Recall, y = Precision)) + 
              geom_line(color = "red") +
              theme_bw(base_size = 16) +
              theme(legend.position="top"))
    } else if (method == "Fscore") {
      print(ggplot(data = roc.plot, aes(x = cutoff, y = fscore)) + 
              geom_line(color = "red") +
              theme_bw(base_size = 16) +
              theme(legend.position="top"))
    } else {
      cat("Method is wring. Use ROC or PR.")
    }
  }
  
  return(roc.plot)
}

recurrent = function(y) {
  require(dplyr)
  
  y %>% 
    mutate(uniq = paste(chrom, pos, ref, alt, sep = ":")) -> y
  
  y %>%
    group_by(uniq) %>%
    summarize(recur = n()) %>%
    inner_join(y, ., by = "uniq") %>%
    select(- uniq) -> y
  
  return(y$recur)
}

driver_mutations_gene = function(mydata, somatic_path, out_path, output_plot = FALSE) {
  # this function counts genes that are predicted somatic, and gives the numbers in true somatic set, then TP, FP, FN,
  # and plots the outcome per gene
  # mydata is a dataframe
  library(dplyr)
  require(ggplot2)
  require(reshape2)
  
  #  args <- commandArgs(TRUE)
  #  variant_path <- args[[1]] # all variants after TOBI filters and predictions
  #  somatic_path <-  args[[2]] # path to true somatic calls
  #  out_path <- args[[3]] # path to final file
  #  print(args)
  #  
  #  mydata <- read.table(variant_path, header = T, stringsAsFactors = F)
  
  # somatic mutations
  somatic = read.table(somatic_path, 
                       sep = "\t", 
                       header = TRUE, 
                       stringsAsFactors = FALSE)
  
  # Choose the corresponding tumor from somatic
  mydata_list = unique(mydata$case)
  somatic = subset(somatic, somatic$case_id %in% mydata_list)
  
  # generate performance values
  time1 = Sys.time()
  #from true somatic calls
  som_notes <- somatic %>% group_by(gene_symbol) %>% mutate(tot_true_som= n()) %>% 
    rename(gene_name=gene_symbol) %>% select(gene_name, tot_true_som) %>% distinct()
  #true calls after TOBI filtering
  true_after_filt <- mydata %>% filter(Y== "som") %>% group_by(gene_name, Y) %>% mutate(true_after_filt= n()) %>%  
    ungroup() %>%  select(gene_name, true_after_filt) %>% distinct()
  #predicted as somatic from TOBI
  tot_pred_som <- mydata %>% filter(pred== "som") %>% group_by(gene_name, pred) %>% mutate(tot_pred_som= n()) %>%  
    ungroup() %>%  select(gene_name, tot_pred_som) %>% distinct()
  #number of unique cases with true somatic mutations in the gene
  true_case <- mydata %>% filter(Y== "som") %>% select(gene_name, case, Y) %>% distinct() %>% 
    group_by(gene_name) %>% mutate(true_case= n()) %>%  
    ungroup() %>%  select(gene_name, true_case) %>% distinct()
  #number of unique cases with predicted somatic mutations in the gene
  pred_case <- mydata %>% filter(pred== "som") %>% select(gene_name, case, pred) %>% distinct() %>% 
    group_by(gene_name) %>% mutate(pred_case= n()) %>%  
    ungroup() %>%  select(gene_name, pred_case) %>% distinct()
  TP <- mydata %>% filter(Y== "som" & pred== "som") %>% group_by(gene_name) %>% mutate(TP= n()) %>% 
    select(gene_name, TP) %>% distinct()
  FP <- mydata %>% filter(Y== "non_som" & pred== "som") %>% group_by(gene_name) %>% mutate(FP= n()) %>% 
    select(gene_name, FP) %>% distinct()
  FN <- mydata %>% filter(Y== "som" & pred== "non_som") %>% group_by(gene_name) %>% mutate(FN= n()) %>% 
    select(gene_name, FN) %>% distinct()
  
  #select longest transcript for protein length since I've lost transcript labelling data
  protein_length <- mydata %>% rename(protein_length=amino_acid_length) %>% 
    select(gene_name, protein_length) %>% distinct() %>% group_by(gene_name) %>% top_n(1,)  
  
  # fraction of synonymous mutations over all
  syn<- mydata %>% filter(effect == "SYNONYMOUS_CODING" | effect == "SYNONYMOUS_STOP") %>% group_by(gene_name) %>% 
    mutate(syn=n()) %>% select(gene_name,syn) %>% distinct()
  tot_length <- mydata %>% group_by(gene_name) %>% mutate(tot_length=n()) %>% select(gene_name,tot_length) %>% distinct()
  syn_ratio <- full_join(syn, tot_length, "gene_name") %>% mutate_each(funs(replace(., which(is.na(.)), 0))) %>% 
    mutate(syn_ratio=syn/tot_length) %>% select(gene_name,syn_ratio) %>% distinct()
  
  # all cases with mutations after filtering
  case_gene <- mydata %>% select(gene_name,case) %>% distinct() %>% group_by(gene_name) %>% 
    mutate(all_cases=n()) %>% select(gene_name,all_cases) %>% distinct()
  
  # join all above categories 
  imp_mt <-  right_join(full_join(full_join(full_join(full_join(full_join(full_join(full_join(full_join(full_join(som_notes, true_after_filt, "gene_name"),
                                                                                                        tot_pred_som, "gene_name"), 
                                                                                              true_case, "gene_name"), 
                                                                                    pred_case, "gene_name"),
                                                                          TP, "gene_name"),
                                                                FP, "gene_name"), 
                                                      FN, "gene_name"), 
                                            protein_length, "gene_name"),
                                  syn_ratio, "gene_name"), 
                        case_gene, "gene_name")  %>%  
    mutate_each(funs(replace(., which(is.na(.)), 0)))
  print(  Sys.time() - time1)
  
  cat("Total number of predicted somatically altered variants:", sum(imp_mt$tot_pred_som), "\n")
  cat("Total true somatically altered variants:", sum(imp_mt$tot_true_som), "\n")
  cat("After filtering and variants:", sum(imp_mt$TP), "\n")
  
  write.table(imp_mt, out_path, 
              sep = "\t", 
              quote = FALSE, 
              na = ".", 
              row.names = FALSE)
  
  if (output_plot) {
    imp_mt_melt = melt(imp_mt, id = c("gene"))
    q = qplot(x = gene, y = value , fill=variable,
              data = imp_mt_melt, geom="bar", stat="identity",
              position="dodge")
    #   q
    q = q + 
      theme_bw() + 
      coord_fixed(ratio = 0.6) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), 
            text = element_text(size = 18),
            legend.title = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position="bottom")
    # q
    print(q)
  }
}