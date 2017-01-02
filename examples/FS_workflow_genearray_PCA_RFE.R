

# original data (array expression dataset)
class(GSE5325) <- "numeric"
features <- GSE5325[,-ncol(GSE5325)]

# response variable (estrogen receptor)
class <- GSE5325[,ncol(GSE5325)]

# getting only those genes with expression level for all samples
m <- features[ , colSums(is.na(features)) == 0]

# compute PCs using inbuilt function R 'prcomp', set scale=FALSE since expression values are already normalized.
m.pca <- prcomp(m, center = TRUE, scale. = FALSE) 

# getting PCs from prcomp object
trainDescr <- m.pca$x
trainClass <- class

rfProfile <- rfe(x = trainDescr, 
                 y = as.factor(trainClass), 
                 sizes = c(1:20,25,30,40,50), 
                 rfeControl = rfeControl(functions = rfFuncs, method="cv", number=10))

# Dataframe with variables metric (i.e. Accuracy) information
results <- rfProfile$results

# plotting rfe object
png("figures/GSE5325_PCA_RFE-RandomForest.png", width = 800, height = 800)
plot.rfe(x = rfProfile)
dev.off()


## Selection by filtering - Linear Discriminat analyis 
sbfProfile <- sbf(x = trainDescr,
                  y = as.factor(trainClass), 
                  sbfControl = sbfControl(functions = ldaSBF, method = "repeatedcv", repeats = 5, saveDetails = TRUE))

summary(sbfProfile)

png("figures/GSE5325_PCA_RFE-LDA.png", width = 800, height = 800)
resampleHist(object = sbfProfile, type = 'histogram')
dev.off()
