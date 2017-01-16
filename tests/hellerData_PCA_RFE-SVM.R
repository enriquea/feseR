library(caret)
library(doMC)

# retrieving aminoacid sequences
sequence.list <- as.list(hellerDataSet$sequence)

# computing AAindex features
temp.list <- lapply(sequence.list, function(x) as.vector(computeAllAAindexValues(x)))

# convert to matrix and rename rows and cols
peptidesFeatures <- do.call(rbind, temp.list)
rownames(peptidesFeatures) <- unlist(sequence.list)
colnames(peptidesFeatures) <- rownames(aaindexMatrix)

# data features (example: Branca dataset with AAindex descriptors estimated)
# use computing_all_aaindex_features.R to create this dataset.
data_features <- peptidesFeatures

# Load data classes (example: Response variable, Branca peptides experimental variable)
data_class <- as.matrix(hellerDataSet)[,2]

# compute PCs using inbuilt function R 'prcomp', set scale=FALSE since expression values are already normalized.
peptides.pca <- prcomp(data_features, center = TRUE, scale. = TRUE) 

# get summary pca.comp object
summ_stat <- summary(peptides.pca)

# getting dataframe with PCs variance information
df_variances <- as.data.frame(t(summ_stat$importance))
#add Principal Component index
df_variances <- cbind(PCs=c(1:nrow(df_variances)), df_variances) 
colnames(df_variances) <- c('PCs','Standard_desviation','Proportion_of_variance','Accumulative_variance')

# holds only PCs until reach variance=1
pcs_subset <- subset(df_variances, Accumulative_variance < 1L)

# plotting histogram from PCs variance
plots_variance <- plotPCVariances(dat = pcs_subset)
png('Principal_component_variance_peptides.png', width = 800, height = 800)
multiplot(plotlist = plots_variance, cols=2)
dev.off()


# Divide the dataset in train and test sets
inTrain <- createDataPartition(data_class, p = 3/4, list = FALSE)

# Create the Training Dataset for Principal components
components <- peptides.pca$x
# getting relevant components
components <- components[,1:20]
trainDescr <- components[inTrain,]

# Create the Testing dataset for Descriptors
testDescr <- components[-inTrain,]

trainClass <- data_class[inTrain]
testClass <- data_class[-inTrain]


registerDoMC(cores = 2)
# caret function: the rfe is the backwards selection, 
# c is the possible sizes of the features sets, 
# and method the optimization method is a support vector machine.
svmProfile <- rfe(x= trainDescr, 
                  y= as.numeric(trainClass), 
                  sizes = c(1:5,10,15), 
                  rfeControl = rfeControl(functions = caretFuncs, number = 2, verbose = TRUE), 
                  method = "svmRadial", 
                  fit = FALSE)

# Dataframe with variables metric (i.e. RMSE) information
resultsSVM <- svmProfile$results

# plotting rfe object
png("figures/peptides_heller_PCA_RFE-SVM.png", width = 800, height = 800)
plot.rfe(x = svmProfile, metric = 'RMSE', output = 'ggplot')
dev.off()