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