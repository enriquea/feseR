library(corrplot)
library(caret)
library(FSelector)
library(doMC)

# load features dataset. Last column is the class variable (experimental pI)
load("~/peptideData.rda")

peptideFeatures <- peptideData[,-ncol(peptideData)]
peptideClass    <- as.numeric(peptideData[,ncol(peptideData)])


#### overrall correlation

#scale all the features 
feature.scale <- scale(peptideFeatures, center=TRUE, scale=TRUE)

#compute the correlation matrix for all features
featureCorr <- cor(feature.scale)



#### univariate correlation filter implementation

# build input dataframe to FSelector
merge <- as.data.frame(cbind(feature.scale , peptideClass))

# compute correlations
weights <- linear.correlation(peptideClass ~., merge)

# filter out features with low correlation (cutoff 0.3)
attr_subset <- subset(weights, attr_importance >= 0.3)
filteredLow <- feature.scale[, colnames(feature.scale) %in% rownames(attr_subset)]

# compute the correlation matrix for filtered (low) peptides features
featureCorrLow <- cor(filteredLow)



#### matrix correlation filter implementation

# getting highly correlated features and pply correlation filter at 0.70 (arbitrary value)
highlyCor <- findCorrelation(featureCorr , 0.70)

#then we remove all the variable correlated with more 0.7.
filteredHigh <- feature.scale[,-highlyCor]

# computed correlation on filtered data
featureCorrHigh <- cor(filteredHigh)



#### merged correlation filters (low and high)

# getting highly correlated features on low-filtered dataset
highlyCor <- findCorrelation(featureCorrLow, 0.70)
filteredFull <- filteredLow[,-highlyCor]
featureCorrFull <- cor(filteredFull)


#### RFE-SVM

features <- filteredHigh
class <- peptideClass

# Divide the dataset in train and test sets
inTrain <- createDataPartition(class, p = 3/4, list = FALSE)

# Create the Training Dataset
trainDescr <- features[inTrain,]

# Create the Testing dataset for Descriptors
testDescr <- features[-inTrain,]

# training class subset
trainClass <- class[inTrain]

# testing class subset
testClass <-  class[-inTrain]

# running in multicore (n*2)
registerDoMC(cores = 2)

# caret function: the rfe is the backwards selection, 
# c is the possible sizes of the features sets, 
# and method the optimization method is a support vector machine.
svmProfile <- rfe(x= trainDescr, 
                  y= as.numeric(trainClass), 
                  sizes = c(1:5, 10, 15, 20, 25), 
                  rfeControl = rfeControl(functions = caretFuncs, number = 2, verbose = TRUE), 
                  method = "svmRadial", 
                  fit = FALSE)

# Dataframe with variables metric (i.e. RMSE) information
resultsSVM <- svmProfile$results

# plotting rfe object
png("figures/peptides_SVM_importance_variables.png", width = 800, height = 800)
plot.rfe(x = svmProfile, metric = 'RMSE', output = 'ggplot')
dev.off()
