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