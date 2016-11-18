## load libraries
require(plyr)


## This function compute amino acid sequence AAindex property given an accession code.
computeAAindexValue <- function(seq, accession){
    
   loadAAindex()
      if (!accession %in% names(aaindex)){
         stop("Not valid accession code...")
      }
     aaindexTable <- getAAindexTable(accession = accession)
  
     sequence <- toupper(seq)
     sequence <- reformatSequence(seq = sequence)
     lev <- c("A",  "L",  "R", "K", "N", "M", "D", "F", "C", "P", "Q", "S", "E", "T", "G", "W", "H", "Y", "I", "V")
     seqTable <- table(factor(prot <- strsplit(sequence, "")[[1]], levels = lev))
  
     temp <- seqTable*aaindexTable
     
     meanIndex <- sum(as.vector(temp))/length(prot)
     return(meanIndex)
}


## This function compute all AAindex values given a sequence

computeAllAAindexValues <- function(seq){
           loadAAindex()
           acc.list <- names(aaindex)
           aaindex.vector <- lapply(acc.list, function(x) computeAAindexValue(seq, x))
           return(aaindex.vector)
}


## Retrieving AAindex table from accesion code
getAAindexTable <- function(accession){
  #load AAindex data
  loadAAindex()
  #retrive table using acc code
  aaindexTable <- as.table(aaindex[[accession]]$I)
  names(aaindexTable) <- mapAA(names(aaindexTable))
  return(aaindexTable)
}


## This function load AAindex database
loadAAindex <- function(){
   load(file = 'data/aaindex.rda', envir = .GlobalEnv)
}
