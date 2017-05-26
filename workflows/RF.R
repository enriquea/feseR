library(randomForest)
library(caret)

source('R/visualization.R') # load viz function

# fs_worflow_description
fs_workflow_description <- 'Naive Random Forest (RF) implementation for variable
                            importance evaluation.'

# assing workflow label (used to save plot, tables)
wf_label <- "RF"

# read original data
load('data/TNBC.rda')

# data
dat <- TNBC

# assing label to data (used to save plot, tables)
dat_label <- "TNBC"

# brief data description
dat_description <- comment(dat)

# retrive features
class(dat) <- "numeric"
features <- dat[,-ncol(dat)]

# retrive class variables (expected last column)
class <- dat[,ncol(dat)]

# getting only those features (i.e. gene/protein expression values) with values for all instances (i.e. samples)
# This step could remove some features!!!
features <- features[ , colSums(is.na(features)) == 0]

# Scale data features
# These transformations coerce the original predictors to have zero mean and standard deviation equal one.
features <- scale(features, center=TRUE, scale=TRUE)

# define workflow metrics
ext_folds <- 20 # folds to repit the entire workflow
nVars <- vector() # number of features in the final model
accv <- vector() # prediction accuracy
test_stats <- data.frame() # test stats summary per runs

set.seed(1)

time_start <- proc.time() # timing the workflow

for(i in seq(1,ext_folds,1)){
  
  # Divide the dataset in train and test sets
  inTrain <- createDataPartition(as.factor(class), p = 2/3, list = FALSE)
  
  # Create the Training Dataset for Descriptors 
  trainDescr <- features[inTrain,]
  
  # Create the Testing dataset for Descriptors
  testDescr <- features[-inTrain,]
  
  trainClass <- class[inTrain]
  testClass <-  class[-inTrain]
  
  ## generating RF model
  rfProfile <- randomForest(x = trainDescr, 
                          y = as.factor(trainClass), 
                          importance = TRUE,
                          keep.forest = TRUE)
  
  ## predict variable response (class) for test sammples with the new model
  predictedClass <- predict(rfProfile, newdata = testDescr)
  
  ## computing confussion matrix
  confMatrix <- confusionMatrix(data = predictedClass, reference = testClass)
  
  ## building data frame with test phase stats
  test_stats <- rbind(test_stats, confMatrix$overall)
  
  ## getting accuracy from confusion matrix
  accv[i] <- confMatrix$overall['Accuracy']
  
  ## retrive number of variables used in the model
  nVars[i] <- length(varUsed(rfProfile, count = FALSE))
  
  ## Keep best model. It combines two conditions: (Maximal accuracy) AND (Minimal variables)
  if (accv[i] >= max(accv) && nVars[i] <= min(nVars[which(accv == max(accv))])){
    bestModel <- rfProfile # replace model only if improve the metrics (accuracy and number of variables)
    best_test_run <- i # used later for highlighting
  }
}

# summarizing workflow metrics
time_end <- (proc.time() - time_start)['elapsed'] # runtime (sec.)
nVar_mean <- mean(nVars)
acc_mean <- mean(accv) # accuracy mean on testing
acc_sd <- sd(accv)     # accuracy satandard deviation on testing

# dataframe with metrics (i.e. Accuracy) information from training phase (only for the best model)
results_training <- bestModel$results

# formatting table with test metrics
names(test_stats) <- names(confMatrix$overall)
test_stats <- cbind(run = c(1:ext_folds), Variables = nVars, test_stats)
test_stats[,c('AccuracyLower','AccuracyUpper','McnemarPValue','AccuracyNull')] <- list(NULL) # remove some metrics

# save best model from training
save(bestModel, file = paste('models/','model_',dat_label,'_',wf_label,'_',Sys.Date(),'.rda',sep = ""))

## plot importance variables
png(file = paste('figures/','plot_variables_',dat_label,'_',wf_label,'_',Sys.Date(),'.png',sep = ""), width = 800, height = 800)
plotImportance <- varImpPlot(bestModel)
dev.off()

# Visualizing samples on principal component before feature selection process
features <- dat[,-ncol(dat)]
class <-    dat[,ncol(dat)]
features <- features[,colSums(is.na(features)) == 0]
# compute PCs using inbuilt function R 'prcomp', 
# set scale=FALSE since expression values are already normalized.
features.pca <- prcomp(features, center = TRUE, scale. = TRUE)
p1 <- plotPCVariances(prcomp = features.pca, groups = class)
png(paste('figures/','PCA_initial_',dat_label,'_',wf_label,'_',Sys.Date(),'.png',sep = ""), width = 1200, height = 1200)
multiplot(plotlist = p1, cols=2)
dev.off()

# Visualizing samples on principal component after feature selection process
setVars <- varUsed(bestModel, count = FALSE)
features <- features[, setVars]
# compute PCs using inbuilt function R 'prcomp', 
# set scale=FALSE since expression values are already normalized.
features.pca <- prcomp(features, center = TRUE, scale. = TRUE)
p2 <- plotPCVariances(prcomp = features.pca, groups = class)
png(paste('figures/','PCA_end_',dat_label,'_',wf_label,'_',Sys.Date(),'.png',sep = ""), width = 1200, height = 1200)
multiplot(plotlist = p2, cols=2)
dev.off()

## Render a final report using pre-defined template
rmarkdown::render('analysis/report_template.Rmd',
                  output_dir = 'analysis',
                  output_file = paste('report','_',dat_label,'_', wf_label,'_',Sys.Date(), '.pdf', sep = ''), 
                  output_format = 'pdf_document')

