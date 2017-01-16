library(corrplot)
library(caret)
library(FSelector)

# original data (array expression dataset)
class(GSE5325) <- "numeric"
data_features <- GSE5325[,-ncol(GSE5325)]

# getting only those genes with expression level for all samples
data_features <- data_features[ , colSums(is.na(data_features)) == 0]

# response variable (estrogen receptor)
data_class <- as.matrix(GSE5325[,ncol(GSE5325)])

# Scale data features
data_features <- scale(data_features, center=TRUE, scale=TRUE)


#compute the correlation matrix for all features
featuresCorr <- cor(data_features)

#visualize and save the matrix, clustering features by correlation index.
png("figures/genes_corrplot_wo_filtering.png", width = 1800, height = 1800)
p1 <- corrplot(featuresCorr, order = "hclust", tl.pos = 'n', tl.cex = 0.2) # Make plot
dev.off()


#### univariate correlation filter implementation


# build input dataframe to FSelector
merge <- as.data.frame(cbind(data_features, data_class))

# compute correlations
weights <- linear.correlation(data_class ~., merge)

# filter out features with low correlation (cutoff 0.3)
attr_subset <- subset(weights, attr_importance >= 0.3)
filteredLow <- data_features[, colnames(data_features) %in% rownames(attr_subset)]

# compute the correlation matrix for filtered (low) peptides features
featuresCorr1 <- cor(filteredLow)

#visualize and save the matrix, clustering features by correlation index.
png("figures/genes_corrplot_univariate_low03.png", width = 20000, height = 20000)
p2 <- corrplot(featuresCorr1, order = "hclust", tl.pos = 'n', tl.cex = 0.8) # Make plot
dev.off()


#### matrix correlation filter implementation


# getting highly correlated features and apply correlation filter at 0.70 (arbitrary value)
highlyCor <- findCorrelation(featuresCorr, 0.70)

#then we remove all the variable correlated with more 0.7.
filteredHigh <- data_features[,-highlyCor]

# computed correlation on filtered data
featuresCorr2 <- cor(filteredHigh)

#visualize and save new filtered data.
png("figures/genes_corrplot_high07.png", width = 20000, height = 20000)
p3 <- corrplot(featuresCorr2, order = "hclust", tl.pos = 'n', tl.cex = 0.8) # Make plot
dev.off()


#### merged correlation filters (low and high)

# getting highly correlated features on low-filtered dataset
highlyCor <- findCorrelation(featuresCorr1, 0.70)
filteredFull <- filteredLow[,-highlyCor]
featuresCorr3 <- cor(filteredFull)
png("figures/genes_corrplot_full_filtered.png", width = 10000, height = 10000)
p4 <- corrplot(featuresCorr3, order = "hclust", tl.pos = 'n', tl.cex = 0.8) # Make plot
dev.off()
