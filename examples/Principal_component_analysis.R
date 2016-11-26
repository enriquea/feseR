
## Principal Component Analysis 
## See https://github.com/chinadd/PCA for data description

# load array expression dataset
load(file = 'data/GSE5325_array_expression_filtered.rda', envir = .GlobalEnv)

# get gene index from original dataset
rnames <- genearray$ID_REF

# get sample names from original dataset and exclude name from first column (gene reference index)
cnames <- names(genearray)[-1]

# convert to matrix
m <- as.matrix(sapply(genearray[,-1], as.numeric))

# rename colums matrix
colnames(m) <- cnames

# rename rows matrix
rownames(m) <- rnames

#replacing 
#m[is.na(m)] <- 0

# getting only those genes with expression level for all samples
m1 <- m[complete.cases(m),]

#apply transpose to get the data in this way: rows=samples, columns=variables(genes).
m2 <- t(m1)

# compute PCs using inbuilt function R 'prcomp', set scale=FALSE since expression values are already normalized.
m.pca <- prcomp(m2, center = TRUE, scale. = FALSE) 

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
