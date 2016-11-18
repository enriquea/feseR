#corrplot: the library to compute correlation matrix.
library(corrplot)
library(caret)

# Branca features dataset
branca <- brancaFeatures

#scale all the features 
branca.scale <- scale(branca, center=TRUE, scale=TRUE)

#compute the correlation matrix
corbranca <- cor(branca.scale)

#visualize and save the matrix, clustering features by correlation index.
png("figures/correlation_plot_wo_filtering.png", width = 1800, height = 1800)
p1 <- corrplot(corbranca, order = "hclust", tl.pos = 'n', tl.cex = 0.2) # Make plot
dev.off()

## getting highly correlated features and pply correlation filter at 0.70 (arbitrary value)
highlyCor <- findCorrelation(corbranca, 0.70)

#then we remove all the variable correlated with more 0.7.
branca.filtered <- branca.scale[,-highlyCor]

# computed correlation on filtered data
corbranca2 <- cor(branca.filtered)

#visualize and save new filtered data.
png("figures/correlation_plot_filtered.png", width = 1800, height = 1800)
p2 <- corrplot(corbranca2, order = "hclust", tl.cex = 0.8) # Make plot
dev.off()
