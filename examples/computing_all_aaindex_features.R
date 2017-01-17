## Computing AAindex features

# loading peptides file
filepath <- system.file("extdata", 'branca_peptides.csv', package = "feseR")
peptideDataSet <- read.table(file = filepath, header = TRUE, sep = ',', stringsAsFactors = FALSE)

# retrieving aminoacid sequences
sequence.list <- as.list(peptideDataSet$sequence)

# computing AAindex features
temp.list <- lapply(sequence.list, function(x) as.vector(computeAllAAindexValues(x)))

# convert to matrix and rename rows and cols
peptideFeatures <- do.call(rbind, temp.list)
rownames(peptideFeatures) <- unlist(sequence.list)
colnames(peptideFeatures) <- rownames(aaindexMatrix)

# save data with all AAindex pre-computed descriptors (mean)
save(peptideFeatures, file = 'data/brancaAAindexFeatures.rda')
