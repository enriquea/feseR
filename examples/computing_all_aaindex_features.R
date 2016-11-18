
## load functions and data
source('R/computeAAindex.R')
source('R/sequenceProcessing.R')

## Computing AAindex features

# loading peptides file
brancaDataSet <- read.table(file = 'data/branca_peptides.csv', header = TRUE, sep = ',', stringsAsFactors = FALSE)

# retrieving aminoacid sequences
sequence.list <- as.list(brancaDataSet$sequence)

# computing AAindex features
temp.list <- lapply(sequence.list, function(x) as.vector(computeAllAAindexValues(x)))

# convert to matrix and rename rows and cols
brancaFeatures <- do.call(rbind, temp.list)
rownames(brancaFeatures) <- unlist(sequence.list)
colnames(brancaFeatures) <- rownames(aaindexMatrix)

# save data with all AAindex pre-computed descriptors (mean)
save(brancaFeatures, file = 'data/brancaAAindexFeatures.rda')
