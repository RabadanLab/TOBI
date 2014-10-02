ML = function(training, testing, my_formula, my_metric = "Fscore", max_min = TRUE, do_table = TRUE, do_plot = FALSE) {
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
  
  gbmGrid = expand.grid(n.trees = (2:4)*50,
                        interaction.depth = 3:7,
                        shrinkage = 0.1,
                        threshold = seq(0.05, 0.95, length.out = 25))
  
  # check if "training" is a list, and convert it if not
  if (is.data.frame(training)) {
    a = list()
    a[[1]] = training
    training = a
    rm(a)
  }
  n = length(training)
  
  modboost = list()
#   set.seed(123)
  for (i in 1:n) {
    modboost[[i]] = train(my_formula,
                          data = training[[1]],
                          trControl = fitControl,
                          method = modelInfo,
                          verbose = FALSE,
                          metric = my_metric,
                          maximize = max_min,
                          tuneGrid = gbmGrid)
  }
  if (do_plot) {
    print(ggplot(modboost[[i]]) + theme_grey(base_size = 18))
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
  
  if (do_table & (n == 1)) {
    cat("- Testing confusion table\n")
    if (testing == FALSE) {
      testing = training[[1]]
      cat("- Testing dataset not available. Training dataset is used.\n")
    }
    print(confusionMatrix(table(unlist(predict(modboost[[i]], testing[,-c(1:11)])), testing[,1], dnn=list('predicted','actual'))))
  }
  
  return(modboost)
}