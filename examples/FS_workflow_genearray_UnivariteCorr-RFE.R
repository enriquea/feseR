library(caret)
library(randomForest)
library(FSelector)

# original data (gene array expression dataset)
load('data/GSE5325_genearray.rda')
class(GSE5325) <- "numeric"
features <- GSE5325[,-ncol(GSE5325)]

# getting only those genes with expression level for all samples (~8000)
features <- features[ , colSums(is.na(features)) == 0]

# getting response variable (estrogen receptor)
class <- as.matrix(GSE5325[,ncol(GSE5325)])

# Scale data features
features <- scale(features, center=TRUE, scale=FALSE)

#### univariate correlation filter implementation

# build input dataframe to FSelector
merge <- as.data.frame(cbind(features, class))
# compute correlations
weights <- linear.correlation(class ~., merge)
# filter out features with low correlation (cutoff 0.3)
attr_subset <- subset(weights, attr_importance >= 0.3)
features <- features[, colnames(features) %in% rownames(attr_subset)]

# define workflow metrics
n <- 20 # folds to repit the entire workflow
tv <- vector() # workflow time
nVars <- vector() # number of features in the final model
accv <- vector() # prediction accuracy


for(i in seq(1,n,1)){
  
  t1 <- proc.time()  
  
  # Divide the dataset in train and test sets
  inTrain <- createDataPartition(as.factor(class), p = 2/3, list = FALSE)
  
  # Create the training dataset
  trainDescr <- features[inTrain,]
  
  # Create the testing dataset
  testDescr <- features[-inTrain,]
  
  # create the training class subset
  trainClass <- class[inTrain]
  
  # create the testing class subset
  testClass <- class[-inTrain]
  
  #### recursive feature elimination-random forest
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
  
  ## compute prediction accuracy
  accTable <- postResample(predictedClass, as.factor(testClass))
  accv[i] <- accTable[['Accuracy']] 
  
  ## keep best model
  if (accv[length(accv)] >= max(accv)){
    bestModel <- rfProfile
  }
  
  ## retrive number of variable used
  nVar <- rfProfile$optsize
  nVars[i] <- nVar
  
  ## timing workflow
  t2 <- proc.time() - t1
  tv[i] <- t2[['elapsed']]
  
}

## summary metrics
nVar_mean <- mean(nVars)
time_mean <- mean(tv)
acc_mean <- mean(accv)
acc_sd <- sd(accv)

# Dataframe with variables metric (i.e. Accuracy) information
results <- bestModel$results

# save results
write.table(results, file = "analysis/metrics_genearray_univariateCorr_RFE-RF.txt", sep = "\t", row.names = FALSE)

# save comple model profile
save(bestModel, file = "profiles/univariate_rfe_rf.rda")

# plotting rfe object
png("figures/GSE5325_univariate03_RFE-RF.png", width = 800, height = 800)
plot.rfe(x = bestModel, xlim = c(-5,105))
dev.off()
