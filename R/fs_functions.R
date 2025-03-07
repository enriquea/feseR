
## Correlation filter function

#' filter.corr
#'
#' Filter features based on its correletion with a predictor variable.
#'
#' @param features A numeric matrix as input.
#' @param class Response variable as numeric vector.
#' @param mincorr Correlation coeficient cutoff used for filtering.
#'
#' @return A matrix filtered.
#'
#' @export

filter.corr <- function(features, class, mincorr = 0.3){

  # Matrix input validation
  valid.matrix(mx = features)

  if(!is.vector(class) & (nrow(features) != length(class))){
     stop('Error processing input data...')
  }

  # hold original input
  input <- features

  # split large dataframe in smaller n-columns dfs (faster implementation)
  list.temp <- lapply(split(as.list(as.data.frame(features)), cut(1:ncol(features), 10, labels = FALSE)), FUN = as.data.frame, check.names = FALSE)

  # apply correlation function
  list.corr <- lapply(list.temp, function(x) FSelector::linear.correlation(class ~., x))

  # unlist and rename
  weights <- plyr::rbind.fill(list.corr)
  rownames(weights) <- unlist(lapply(list.corr, row.names))

  # filter out features with low correlation
  attr_subset <- subset(weights, attr_importance >= mincorr)
  features <- features[, colnames(features) %in% rownames(attr_subset)]

  # compare input/output
  compare.matrix(input.matrix = input, output.matrix = features, description = 'Number of removed features (Univariate correlation filter):')

  # return filtered features (matrix)
  return(as.matrix(features))
}



## Gain information filter

#' filter.gain.inf
#'
#' Gain information filter implementation.
#'
#' @param features A numeric matrix as input.
#' @param class Response variable as numeric vector.
#' @param n.percent The percent of features with more gain to be returned.
#' @param zero.gain.out Is TRUE (default), zero-gain features will be filtered out (n.percent will be ignored).
#'
#' @return A matrix filtered.
#'
#' @export

filter.gain.inf <- function(features, class, n.percent = 0.75, zero.gain.out = TRUE) {

  # Matrix input validation
  valid.matrix(mx = features)

  if(!is.vector(class) & (nrow(features) != length(class))){
    stop('Error processing input data...')
  }

  # hold original input
  input <- features

  # split large dataframe in smaller (10-folds) ones
  list.temp <- lapply(split(as.list(as.data.frame(features)),
                            cut(1:ncol(features), 10, labels = FALSE)),
                            FUN = as.data.frame, check.names = FALSE)

  # apply correlation function
  list.gains <- lapply(list.temp, function(x) FSelector::information.gain(class ~., x))

  # unlist and rename
  weights <- plyr::rbind.fill(list.gains)
  rownames(weights) <- unlist(lapply(list.gains, row.names))

  # filter out features with zero-gain
  if(isTRUE(zero.gain.out)){
    attr_subset <- subset(weights, attr_importance > 0)
  } else {
    attr_subset <- cutoff.k.percent(weights, n.percent)
  }

  # filter out features with low gain information
  features <- features[, colnames(features) %in% rownames(attr_subset)]

  # compare input/output
  compare.matrix(input.matrix = input, output.matrix = features, description = 'Number of removed features (Gain Information filter):')

  # return filtered features (matrix)
  return(as.matrix(features))
}


## filter matrix correlation


#' filter.matrix.corr
#'
#' Compute matrix correlation between features and filter using a threshold.
#'
#' @param features A numeric matrix as input.
#' @param maxcorr Correlation coeficient cutoff used for filtering.
#'
#' @return A matrix filtered.
#'
#' @export


filter.matrix.corr <- function(features, maxcorr = 0.75){

  # Matrix input validation
  valid.matrix(mx = features)

  # hold original input
  input <- features

  # computed correlation
  descrCorr <- cor(features)

  # find correlated features
  highCorr <- caret::findCorrelation(descrCorr, maxcorr)

  # filter correlated features
  if(length(highCorr)) {
    features <- features[, -highCorr]
  }

  # compare input/output
  compare.matrix(input.matrix = input, output.matrix = features, description = 'Number of removed features (Matrix correlation filter)')

  return (as.matrix(features))
}


## filter PCA

#' filter.pca
#'
#' Compute Principal Components from input matrix. Control how many components must be returned.
#'
#' @param features A numeric matrix as input.
#' @param center Must the features be centered (default TRUE).
#' @param scale  Must the features be scaled (default TRUE).
#' @param cum.var.cutoff Cumulative variance proportion cutoff. Used to return the Principal Componets (PCs) which satisfies that value. If 1 (default), all PCs will be returned.
#'
#' @return A matrix of Principal components.
#' 
#' @export


filter.pca <- function(features, center = TRUE, scale = TRUE, cum.var.cutoff = 1){

     # Matrix input validation
     valid.matrix(mx = features)

     # hold original input
     input <- features

     # compute principal componets
     features.pca <- prcomp(features, center = center, scale. = scale)

     # getting componet variance distribution
     varianceMatrix <- summary(features.pca)$importance

     # getting cumulative variance proportion
     cvp <- varianceMatrix[c('Cumulative Proportion'),]

     # removing low-variance componets
     highVariancePCs <- cvp[cvp <= cum.var.cutoff]

     # filtering PCs
     pcs <- features.pca$x[,names(highVariancePCs)]

     # print number of final components
     message(paste('Features reduced to:', ncol(pcs), 'Principal components', sep = ' '))

     # return features
     return(as.matrix(pcs))
}



#' rfeRF
#'
#' Recursive feature elimination (RFE) method wrapped with a Random Forest (RF) algorithm for feature importance evaluation.
#'
#' @param features A numeric matrix as input.
#' @param class Response variable as numeric vector. It will be coerced to factor.
#' @param number.cv Number of cross-validation folds (10 default). Used during training phase.
#' @param group.sizes  A numeric vector of integers corresponding to the number of features that should be retained.
#' @param metric Metric to evaluate performance ('Accuracy' (default), 'Kappa' or 'ROC').
#' @param tolerance Allow tolerance for evaluation metric (Default zero).
#' @param verbose Make the output verbose. 
#'
#' @return A list the elements. See \code{\link[caret]{rfe}} for more details.
#'
#'
#' @export
#' 
rfeRF = function(features, class, number.cv = 10, 
                 group.sizes = c(1:10, seq(15, 100, 5)), 
                 metric = "Accuracy", verbose = TRUE, tolerance = 0) {
    
    # Matrix input validation
    valid.matrix(mx = features)
  
    if (!is.vector(class) & (nrow(features) != length(class))) {
        stop("Error processing input data...")
    }
  
    funcs = rfFuncs2
    if (metric == "ROC") {
        funcs$summary = twoClassSummary  
    }
    if(tolerance != 0) {
      funcs$selectSize = function (x, metric, tol = tolerance, maximize) {
        if (!maximize) {
          best <- min(x[, metric])
          perf <- (x[, metric] - best)/best * 100
          flag <- perf <= tol
        }
        else {
          best <- max(x[, metric])
          perf <- (best - x[, metric])/best * 100
          flag <- perf <= tol
        }
        min(x[flag, "Variables"])
      }
    }
    
    #### recursive feature elimination-random forest
    rfProfile <- caret::rfe(x = features, 
                     y = as.factor(class), 
                     maximize = TRUE, 
                     metric = metric, 
                     sizes = group.sizes, 
                     rfeControl = rfeControl(functions = funcs, 
                                             method = "cv", 
                                             number = number.cv,
                                             allowParallel = TRUE,
                                             verbose = verbose))
    return(rfProfile) # return RF profile
}                       
                       


## combine Feature selection methods

#' combineFS
#'
#' The main function controlling the Feature Selection workflow. This function combines sequencially different FS step
#' keeping the following structure: Univariate filter -> Multivariate filter -> Wrapper method. It includes an internal
#' cross-validation step used in the last FS step. In addition, it is posible to set up an external loop which operates
#' over randomized and class-balanced test data. The process returns information from both training and testing phases.
#'
#' @param features A numeric matrix as input.
#' @param class Response variable as numeric vector. It will be coerced to factor.
#'
#' @param univariate Descrition of the Univariate filter to be used. Set 'corr' (default) for correlation filter, 'gain' for gain information or none.
#' @param mincorr The threshold controling the Univariate correlation filter.
#' @param n.percent If 'gain' is selected, this parameter controls the percent of features (with higher) to be returned.
#' @param zero.gain.out Is TRUE (default), zero-gain features will be filtered out (n.percent will be ignored).
#'
#' @param multivariate Multivariate filter to be used. Set 'mcorr' (default) for correlation filter, 'pca' for Principal Component Analysis or none.
#' @param maxcorr The threshold controling the matrix correlation filter (default value 0.75).
#' @param cum.var.cutoff If 'pca' is selected, this parameter controls the PCA process. See function \code{\link{filter.pca}}.
#'
#' @param wrapper Wrapper method to be used. Set 'rfe.rf' (default) for recursive feature elimination wrapped with random forest.
#' @param number.cv See \code{\link{rfeRF}} for description.
#' @param group.sizes See \code{\link{rfeRF}} for description.
#' @param metric Metric to evaluate performance ('Accuracy' (default), 'Kappa' or 'ROC').
#' @param tolerance Allow tolerance for evaluation metric (Default zero).
#' @param verbose Make the output verbose. 
#'
#' @param extfolds Number of times (default 10) to repeat the entire FS process randomizing the dataset (test/training).
#' @param partition Parameter controling the data partition in test and training dataset. It generates random and class-balanced dataset.
#'
#' @return A list with the following elements.
#'
#' opt.variables vector with optimal (final) features names.
#' training dataframe with the metrics from the training phase.
#' testing dataframe with the metrics from the testing phase.
#' best.model the best model obtained (max. accuracy and min. number of final features).
#' runtime the workflow runtime (secs).
#'
#'
#' @export
#' 
combineFS = function(features, class, univariate = "corr", mincorr = 0.3, 
                     n.percent = 0.75, zero.gain.out = TRUE, multivariate = "mcorr", 
                     maxcorr = 0.75, cum.var.cutoff = 1, wrapper = "rfe.rf", 
                     number.cv = 10, group.sizes = c(1:10, seq(15, 100, 5)), 
                     extfolds = 10, partition = 2/3, metric = "Accuracy", tolerance = 0, verbose = TRUE) 
{
  
  set.seed(123)
  time_start <- proc.time()
  list.process <- list()
  
  nVars <- vector()
  accv <- vector()
  # profile = list()
  test_stats <- NULL
  
  results <- foreach(i = seq(1, extfolds, 1), 
                     .packages = c("caret", "ROCR", "doParallel"), 
                     .export = c("rfeRF", "valid.matrix", "rfFuncs2"),
                     .combine = 'comb',
                     .multicombine = TRUE) %dopar% {
                       
                       message("Repeat process number ", i)
                       
                       inTrain <- createDataPartition(as.factor(class), p = partition, list = FALSE)
                       trainDescr <- features[inTrain, ] # training set
                       testDescr <- features[-inTrain, ] # test set
                       trainClass <- class[inTrain] # response training
                       testClass <- class[-inTrain] # response test
                       
                       message(paste0("Step 1/3: Applying univariate filter: ", univariate))
                       
                       if (univariate == "corr") {
                         trainDescr <- filter.corr(trainDescr, trainClass, mincorr = mincorr)
                       }
                       else if (univariate == "gain") {
                         trainDescr <- filter.gain.inf(trainDescr, trainClass, zero.gain.out = zero.gain.out)
                       }
                       else if (univariate == 'none') {
                         message("Univariate filter skipped...")
                       } else {
                         stop("Undefined univariate filter...")
                       }
                       
                       message(paste0("Step 2/3: Applying multivariate filter: ", multivariate))
                       
                       if (multivariate == "mcorr") {
                         trainDescr <- filter.matrix.corr(trainDescr, maxcorr = maxcorr)
                       }
                       else if (multivariate == "pca") {
                         # TODO 
                         # trainDescr <- filter.pca(trainDescr, cum.var.cutoff = cum.var.cutoff)
                         # if PCA is used as multivariate filter, feautures in test data need to be
                         # also converted to componentes (apply with out any filtering).
                         # testDescr <- filter.pca(testDescr, cum.var.cutoff = 1) # no filtering (return all PCs)
                         warning("PCA currently disabled in this setting. PCA filter skipped...")
                       }
                       else if (multivariate == 'none') {
                         message("Multivariate filter skipped...")
                       }
                       else {
                         stop("Undefined multivariate filter...")
                       }
                       
                       if (wrapper == "rfe.rf") {
                         
                         message("Step 3/3: Recursive feature elimination wrapped with random forest model.  It could take a while...")
                         
                         profile = rfeRF(features = trainDescr, class = trainClass, 
                                         number.cv = number.cv, group.sizes = group.sizes,
                                         metric = metric, verbose = verbose, tolerance = tolerance)
                       }
                       else if (wrapper == "ga.rf") {
                         profile = gaRF(features = trainDescr, class = trainClass)
                       }
                       else {
                         stop("Undefined wrapper parameter...")
                       }
                       predictedClass <- predict(profile, newdata = testDescr)
                       if (metric %in% c("Accuracy", "Kappa")) {
                         confMatrix <- confusionMatrix(data = factor(predictedClass$pred), 
                                                       reference = factor(testClass))
                         res <-  confMatrix$overall
                         accv <- confMatrix$overall[metric]
                       } else if (metric == "ROC") {
                         predictions = as.vector(predictedClass[,1])
                         pred = prediction(as.numeric(predictions), testClass)
                         
                         perf_AUC = performance(pred, "auc") #Calculate the AUC value
                         AUC = perf_AUC@y.values[[1]]
                         
                         ss = performance(pred, "sens", "spec")
                         ss = c(ss@x.values[[1]][2], ss@y.values[[1]][2])
                         
                         accv <- AUC
                         res <- as.data.frame(t(c(AUC, ss)))
                         colnames(res) <- c("ROC", "Sens", "Spec")
                         
                         
                       }
                       
                       nVars <- profile$optsize
                       return(list(accv, nVars, res, profile))
                     }
  
  accv = unlist(results[,1])
  nVars = unlist(results[,2])
  test_stats = do.call(rbind, results[,3])
  profile = do.call(list,results[,4])
  for (i in seq(1, extfolds, 1)) {
    if (accv[i] >= max(accv)*(1 - tolerance/100) && nVars[i] <= min(nVars[which(accv >= 
                                                                                max(accv)*(1 - tolerance/100))])) {
      bestModel <- profile[[i]]
    }
  }
  
  time_end <- (proc.time() - time_start)["elapsed"]
  #names(test_stats) <- names(results[, 3]) #names(confMatrix$overall)
  test_stats <- cbind(Run = c(1:extfolds), Variables = nVars, 
                      test_stats)
  test_stats <- test_stats[, colnames(test_stats) %in% 
                             c("Run", "Variables", "Accuracy", "Kappa", "AccuracyPValue", "ROC", "Sens", "Spec")]
  results_training <- bestModel$results
  opt.variables <- bestModel$optVariables
  list.process <- list(opt.variables = opt.variables, 
                       training = results_training, 
                       testing = test_stats, 
                       best.model = bestModel, 
                       tolerance = paste0(tolerance, "%"), 
                       runtime = time_end)
  message("Process finalized!!!")
  return(list.process)
}  


# Helper function to combine results in foreach                       
comb <- function(...) {
    mapply('list', ..., SIMPLIFY=TRUE)
}

                      
# Fix for error in rfFuncs() in caret package when rerank=TRUE
rfFuncs2 <- caret::rfFuncs
rfFuncs2$fit <- function(x, y, first, last, ...) {
    loadNamespace("randomForest")
    randomForest::randomForest(x, y, importance = T,...)
}

rfFuncs2$functions$rank <- function(object, x, y) {
    vimp <- varImp(object)
    if (is.factor(y)) {
        if (all(levels(y) %in% colnames(vimp))) {
            avImp <- apply(vimp[, levels(y), drop = FALSE], 1, mean)
            vimp$Overall <- avImp
        }
    }
    vimp <- vimp[order(vimp$Overall, decreasing = TRUE), , drop = FALSE]
    vimp$var <- rownames(vimp)
    vimp
}

rfFuncs2$rank <- function(object, x, y) {
    vimp <- varImp(object)
    if (is.factor(y)) {
        if (all(levels(y) %in% colnames(vimp))) {
            avImp <- apply(vimp[, levels(y), drop = FALSE], 1, 
                           mean)
            vimp$Overall <- avImp
        }
    }
    vimp <- vimp[order(vimp$Overall, decreasing = TRUE), , drop = FALSE]
    vimp$var <- rownames(vimp)
    vimp
}

  
