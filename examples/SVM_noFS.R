library(e1071)
library(hydroGOF)
library(pROC)
library(pIR)

load("data/peptideData.rda")
start.time <- Sys.time()
ptm <- proc.time()

peptideFeatures <- peptideData[,-ncol(peptideData)]
feature.scale <- scale(peptideFeatures, center=TRUE, scale=TRUE)

peptideClass    <- as.numeric(peptideData[,ncol(peptideData)])

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

feature.scale <- filteredFull


svm_model_after_tune <- svm(peptideClass ~ ., data=feature.scale, kernel="radial", cost=1, gamma=0.5)
summary(svm_model_after_tune)
pred <- predict(svm_model_after_tune,feature.scale)

rmse(pred, peptideClass)
cor(pred, peptideClass)
proc.time() - ptm
time.taken <- Sys.time() - start.time
time.taken
