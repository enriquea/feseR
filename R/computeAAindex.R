## load libraries
require(plyr)

## This function compute amino acid sequence AAindex property given an accession code.
computeAAindexValueByAccession <- function(seq, accession){
  
      if (!accession %in% names(aaindex)){
         stop("Not valid accession code...")
      }
     aaindexTable <- getAAindexTable(accession = accession)
  
     sequence <- toupper(seq)
     sequence <- reformatSequence(seq = sequence)
     lev <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
     seqTable <- table(factor(prot <- strsplit(sequence, "")[[1]], levels = lev))
  
     temp <- seqTable*aaindexTable
     
     meanIndex <- sum(as.vector(temp))/length(prot)
     return(meanIndex)
}


## This function compute all AAindex values given a sequence

computeAllAAindexValues <- function(seq){
   #loadAAindex()
      sequence <- toupper(seq)
      sequence <- reformatSequence(seq = sequence)
      lev <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
      seqTable <- table(factor(naa <- strsplit(sequence, "")[[1]], levels = lev))
      # AAindex matrix
      M <- aaindexMatrix
      # sequence vector
      P <- as.vector(seqTable)
      
      features <- M %*% P
      features <- features*(1/length(naa))     
   return(features)
}


## Retrieving AAindex table from accesion code
getAAindexTable <- function(accession){
  #retrive table using acc code
  aaindexTable <- as.table(aaindex[[accession]]$I)
  names(aaindexTable) <- mapAA(names(aaindexTable))
  return(aaindexTable)
}


## Reducing AAindexDataBase and getting simplified output like matrix
getAAindexMatrix <- function(){
   #loadAAindexMatrix() 
   accs <- names(aaindex)
   aaindexVectors <- lapply(accs, function(x) as.vector(aaindex[[x]]$I))
   aaindexMatrix <- do.call(rbind, aaindexVectors)
   rnames <- accs
   cnames <- names(aaindex[[1]]$I)
   rownames(aaindexMatrix) <- rnames
   colnames(aaindexMatrix) <- mapAA(cnames)
   # replace NA with zeros
   aaindexMatrix[is.na(aaindexMatrix)] <- 0
   return(aaindexMatrix)
}


## This function load AAindex database
loadAAindexDataBase <- function(){
   load(file = 'data/aaindex.rda', envir = .GlobalEnv)
}

## This function load AAindexMatrix to fast computation.
loadAAindexMatrix <- function(){
   load(file = 'data/aaindexMatrix.rda', envir = .GlobalEnv)
}

## loading aaindex database.
loadAAindexDataBase()

## loading aaindex simplified matrix.
loadAAindexMatrix()