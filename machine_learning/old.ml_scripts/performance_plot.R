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