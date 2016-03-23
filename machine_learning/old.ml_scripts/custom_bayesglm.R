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