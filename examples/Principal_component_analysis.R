
## Principal Component Analysis 
## See https://github.com/chinadd/PCA for data description

# original data (array expression dataset)
class(GSE5325) <- "numeric"
features <- GSE5325[,-ncol(GSE5325)]

# response variable (estrogen receptor)
class <- GSE5325[,ncol(GSE5325)]

# getting only those genes with expression level for all samples
m <- features[ , colSums(is.na(features)) == 0]

# compute PCs using inbuilt function R 'prcomp', set scale=FALSE since expression values are already normalized.
m.pca <- prcomp(m, center = TRUE, scale. = FALSE) 

#print(m.pca)
#plot(m.pca, type = "l")

# get summary pca.comp object
summ_stat <- summary(m.pca)

# getting dataframe with PCs variance information
df_variances <- as.data.frame(t(summ_stat$importance))
#add Principal Component index
df_variances <- cbind(PCs=c(1:nrow(df_variances)), df_variances) 
colnames(df_variances) <- c('PCs','Standard_desviation','Proportion_of_variance','Accumulative_variance')

# plotting histogram from PCs variance
plots_variance <- plotPCVariances(dat = df_variances)
png('figures/Principal_component_variance.png', width = 800, height = 800)
multiplot(plotlist = plots_variance, cols=2)
dev.off()
