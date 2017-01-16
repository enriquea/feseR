library(randomForest)
library(caret)

# Variable selection with random forest

# original data (array expression dataset)
class(GSE5325) <- "numeric"

features <- GSE5325[,-ncol(GSE5325)]

# response variable (estrogen receptor)
class <- GSE5325[,ncol(GSE5325)]

# getting only those genes with expression level for all samples
features <- features[ , colSums(is.na(features)) == 0]

# Scale data features
features <- scale(features, center=TRUE, scale=FALSE)

# define workflow metrics
n <- 20 # folds to repit the entire workflow
tv <- vector() # workflow time
nVars <- vector() # number of features in the final model
accv <- vector() # prediction accuracy



for(i in seq(1,n,1)){
  
  t1 <- proc.time()  
  
  # Divide the dataset in train and test sets
  inTrain <- createDataPartition(as.factor(class), p = 2/3, list = FALSE)
  
  # Create the Training Dataset for Descriptors 
  trainDescr <- features[inTrain,]
  
  # Create the Testing dataset for Descriptors
  testDescr <- features[-inTrain,]
  
  trainClass <- class[inTrain]
  testClass <-  class[-inTrain]
  
  ## generating RF model
  rfModel <- randomForest(x = trainDescr, 
                          y = as.factor(trainClass), 
                          importance = TRUE,
                          keep.forest = TRUE)
  
  ## variable importance rank
  # list.var <- importance(rfModel)
  
  ## predict variable response (class) for all sammples with the new model
  predictedClass <- predict(rfModel, newdata = testDescr)
  
  ## compute prediction accuracy
  accTable <- postResample(predictedClass, as.factor(testClass))
  accv[i] <- accTable[['Accuracy']] 
  
  ## keep best model
  if (accv[length(accv)] >= max(accv)){
    bestModel <- rfModel
  }
  
  ## retrive number of variable used
  nVar <- length(varUsed(rfModel, count = FALSE))
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


## plot importance variables
png(filename = 'figures/gene_expression_randomforest.png', width = 800, height = 800)
plotImportance <- varImpPlot(bestModel)
dev.off()
