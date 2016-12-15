library(caret)
library(randomForest)

## Principal Component Analysis 
## See https://github.com/chinadd/PCA for data description

# original data (array expression dataset)
class(GSE5325) <- "numeric"
data_features <- GSE5325[,-ncol(GSE5325)]

# getting only those genes with expression level for all samples
data_features <- data_features[ , colSums(is.na(data_features)) == 0]

# response variable (estrogen receptor)
data_class <- as.matrix(GSE5325[,ncol(GSE5325)])

# Scale data features
data_features <- scale(data_features, center=TRUE, scale=TRUE)

# Divide the dataset in train and test sets
inTrain <- createDataPartition(data_class, p = 3/4, list = FALSE)

# Create the Training Dataset for Descriptors 
trainDescr <- data_features[inTrain,]

# Create the Testing dataset for Descriptors
testDescr <- data_features[-inTrain,]

trainClass <- data_class[inTrain]
testClass <- data_class[-inTrain]

# Here, we can included a correlation matrix analysis to remove 
# the redundant features before to perform the backwards selection approach.
descrCorr <- cor(trainDescr)
highCorr <- findCorrelation(descrCorr, 0.5)
trainDescr <- trainDescr[, -highCorr]
testDescr <- testDescr[, -highCorr]

#### Example 1. RFE-SVM

# caret function: the rfe is the backwards selection, 
# c is the possible sizes of the features sets, 
# and the optimization method is a support vector machine.

# Warning (If regression): You are trying to do regression and your outcome only has two possible values Are you trying to do classification? 
# If so, use a 2 level factor as your outcome column.
svmProfile <- rfe(x = trainDescr, 
                  y = as.factor(trainClass), 
                  sizes = c(1:5,10,20,30,50,100,200,300,400,500), 
                  rfeControl = rfeControl(functions = caretFuncs, number = 10), 
                  method = "svmRadial", 
                  fit = FALSE)

# Dataframe with variables metric (i.e. RMSE) information
resultTable <- svmProfile$results

# plotting rfe object
png("figures/GSE5325_importance_variables.png", width = 800, height = 800)
plot.rfe(x = svmProfile, metric = 'RMSE', output = 'ggplot')
dev.off()


#### Example 2. RFE-RandomForest

# caret function: the rfe is the backwards selection, 
# c is the possible sizes of the features sets, 
# and the optimization method is random forest.
# Note that response variable must be converted to factor since
# this is a classification problem (non-regression)

rfProfile <- rfe(x = trainDescr, 
                 y = as.factor(trainClass), 
                 sizes = c(1:5,10,20,30,50,100,200,300,400,500), 
                 rfeControl = rfeControl(functions = rfFuncs, method="cv", number=10))

# Dataframe with variables metric (i.e. Accuracy) information
results <- rfProfile$results

# plotting rfe object
png("figures/GSE5325_importance_variables_random_forest.png", width = 800, height = 800)
plot.rfe(x = rfProfile)
dev.off()


#### Example 3. RFE-LDA

## Selection by filtering - Linear Discriminat analyis 
sbfProfile <- sbf(x = trainDescr,
                  y = as.factor(trainClass), 
                  sbfControl = sbfControl(functions = ldaSBF, method = "repeatedcv", repeats = 5, saveDetails = TRUE))

summary(sbfProfile)

png("figures/GSE5325_hist_selection_by_filtering_LDA.png", width = 800, height = 800)
resampleHist(object = sbfProfile, type = 'histogram')
dev.off()

