library(caret)
library(randomForest)
library(FSelector)
library(plyr)

source('R/visualization.R') # load viz function

# Brief description about the Feature Selection workflow
fs_workflow_description <- 'Gain Information filter (GI)
                            with Multivariate Correlation filter (CM)
                            follow by Recursive Feature Elimination (RFE) 
                            wrapped with Random Forest (RF).'

# assing workflow label (used to save plot, tables)
wf_label <- "GI_CM_RFE_RF"

# read original data
load('data/GSE6919_GPL93.rda')

# data
dat <- GSE6919_GPL93

# assing label to data (used to save plot, tables)
dat_label <- "GSE6919_GPL93"

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

## Gain Information filter implementation

# split large dataframe in smaller (10-folds) ones
list.temp <- lapply(split(as.list(as.data.frame(features)), 
                          cut(1:ncol(features), 10, labels = FALSE)), 
                    FUN = as.data.frame, check.names = FALSE)

# apply correlation function
list.gains <- lapply(list.temp, function(x) information.gain(class ~., x))

# unlist and rename
weights <- rbind.fill(list.gains)
rownames(weights) <- unlist(lapply(list.gains, row.names))

# filter out features with zero-gain
#attr_subset <- cutoff.k.percent(weights, .75)
attr_subset <- subset(weights, attr_importance > 0)

#attr_subset <- subset(weights, attr_importance >= 0.3)
features <- features[, colnames(features) %in% rownames(attr_subset)]


#### Multivariate Correlation filter
descrCorr <- cor(features)
highCorr <- findCorrelation(descrCorr, 0.75)
features <- features[, -highCorr]


# define workflow metrics
ext_folds <- 20 # folds to repit the entire workflow
nVars <- vector() # number of features in the final model
accv <- vector() # prediction accuracy
test_stats <- data.frame() # test stats summary per runs

set.seed(1)

time_start <- proc.time() # timing the workflow

# training process and prediction on test data
for(i in seq(1,ext_folds,1)){
  
  # Divide the dataset in train and test sets
  # This process is random and class-balanced
  inTrain <- createDataPartition(as.factor(class), p = 2/3, list = FALSE)
  
  # Create the training dataset
  trainDescr <- features[inTrain,]
  
  # Create the testing dataset
  testDescr <- features[-inTrain,]
  
  # create the training class subset
  trainClass <- class[inTrain]
  
  # create the testing class subset
  testClass <- class[-inTrain]
  
  # recursive feature elimination wrapped with random forest
  rfProfile <- rfe(x = trainDescr, 
                   y = as.factor(trainClass),
                   maximize = TRUE,
                   metric = 'Accuracy',
                   sizes = c(1:10,15,20,25,30,35,40,45,50,60,70,80,90,100), 
                   rfeControl = rfeControl(functions = rfFuncs, 
                                           method = "cv",
                                           number = 10,
                                           verbose = TRUE))
  
  ## predict variable response (class) for test sammples with the new model
  predictedClass <- predict(rfProfile, newdata = testDescr)
  
  ## computing confussion matrix
  confMatrix <- confusionMatrix(data = predictedClass$pred, reference = testClass)
  
  ## building data frame with test phase stats
  test_stats <- rbind(test_stats, confMatrix$overall)
  
  ## getting accuracy from confusion matrix
  accv[i] <- confMatrix$overall['Accuracy']
  
  ## retrive number of variables used in the model
  nVars[i] <- rfProfile$optsize
  
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

# save accuracy table from training
write.table(results_training, file = paste('analysis/','acc_table_',dat_label,'_',wf_label,'_',Sys.Date(),'.txt',sep = ""), 
            sep = "\t", 
            row.names = FALSE)

# save best model from training
save(bestModel, file = paste('models/','model_',dat_label,'_',wf_label,'_',Sys.Date(),'.rda',sep = ""))

# plot variables subset vs. accuracy
png(file = paste('figures/','plot_variables_',dat_label,'_',wf_label,'_',Sys.Date(),'.png',sep = ""), width = 800, height = 800)
plot.rfe(x = bestModel, xlim = c(-5,105))
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
setVars <- bestModel$optVariables
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

## clening up environment
# rm(list = ls())
