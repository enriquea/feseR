
library(corrplot)
library(caret)
library(FSelector)

# Branca features dataset
peptideFeatures <- peptideFeatures

# peptide class
peptideClass <- as.vector(peptideDataSet)[,2]

#scale all the features 
peptides.scale <- scale(peptideFeatures, center=TRUE, scale=TRUE)

#compute the correlation matrix for all features
corpeptides <- cor(peptides.scale)

#visualize and save the matrix, clustering features by correlation index.
png("figures/peptides_corrplot_wo_filtering.png", width = 1800, height = 1800)
p1 <- corrplot(corpeptides, order = "hclust", tl.pos = 'n', tl.cex = 0.2) # Make plot
dev.off()


#### univariate correlation filter implementation


# build input dataframe to FSelector
merge <- as.data.frame(cbind(peptides.scale, peptideClass))

# compute correlations
weights <- linear.correlation(peptideClass ~., merge)

# filter out features with low correlation (cutoff 0.3)
attr_subset <- subset(weights, attr_importance >= 0.3)
low.filtered.peptides <- peptides.scale[, colnames(peptides.scale) %in% rownames(attr_subset)]

# compute the correlation matrix for filtered (low) peptides features
corpeptides1 <- cor(low.filtered.peptides)

#visualize and save the matrix, clustering features by correlation index.
png("figures/peptides_corrplot_univariate_low03.png.png", width = 1800, height = 1800)
p1 <- corrplot(corpeptides1, order = "hclust", tl.cex = 0.8) # Make plot
dev.off()


#### matrix correlation filter implementation


# getting highly correlated features and pply correlation filter at 0.70 (arbitrary value)
highlyCor <- findCorrelation(corpeptides, 0.70)

#then we remove all the variable correlated with more 0.7.
peptides.filtered <- peptides.scale[,-highlyCor]

# computed correlation on filtered data
corpeptides2 <- cor(peptides.filtered)

#visualize and save new filtered data.
png("figures/peptides_corrplot_high07.png", width = 1800, height = 1800)
p2 <- corrplot(corpeptides2, order = "hclust", tl.cex = 0.8) # Make plot
dev.off()


#### merged correlation filters (low and high)

# getting highly correlated features on low-filtered dataset
highlyCor <- findCorrelation(corpeptides1, 0.70)
peptides.filtered <- low.filtered.peptides[,-highlyCor]
corpeptides3 <- cor(peptides.filtered)
png("figures/peptides_corrplot_full_filtered.png", width = 1800, height = 1800)
p2 <- corrplot(corpeptides3, order = "hclust", tl.cex = 0.8) # Make plot
dev.off()