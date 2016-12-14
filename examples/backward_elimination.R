library(caret)

# data features (example: Branca dataset with AAindex descriptors estimated)
# use computing_all_aaindex_features.R to create this dataset.
data_features <- brancaFeatures

# Load data classes (example: Response variable, Branca peptides experimental variable)
data_class <- as.matrix(brancaDataSet)[,2]

# Scale data features
data_features<- scale(data_features, center=TRUE, scale=TRUE)

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

# caret function: the rfe is the backwards selection, 
# c is the possible sizes of the features sets, 
# and method the optimization method is a support vector machine.
svmProfile <- rfe(x=trainDescr, 
                  y=trainClass, 
                  sizes=c(1:5), 
                  rfeControl=rfeControl(functions = caretFuncs, number = 2), 
                  method = "svmRadial", 
                  fit = FALSE)
