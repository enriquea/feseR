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


## renaming amino acid from 3-letter to 1-letter code

mapAA <- function(seqvec){
  seqvec <- revalue(seqvec, c("Ala"="A", "Arg"="R", "Asn"="N", "Asp"="D", "Cys"="C", 
                              "Gln"="Q", "Glu"="E", "Gly"="G", "His"="H", "Ile"="I", 
                              "Leu"="L", "Lys"="K", "Met"="M", "Phe"="F", "Pro"="P",
                              "Ser"="S", "Thr"="T", "Trp"="W", "Tyr"="Y", "Val"="V"))
  return(seqvec)
}


## This function load AAindex database
loadAAindex <- function(){
   load(file = 'data/aaindex.rda', envir = .GlobalEnv)
}

#' reformatSequence
#'
#' This function reformat the sequence to remove inconsistencies
#'
#' @param seq sequence
#'


reformatSequence <- function(seq){
  seq <- gsub("X", "", seq)
  seq <- gsub("B", "", seq)
  seq <- gsub("J", "", seq)
  seq <- gsub("Z", "", seq)
  seq <- gsub("U", "", seq)
  return (seq)
}

#' removePTM
#'
#' This function reformat the sequence to remove all PTM labers.
#' Useful for example to compute AAindex Isoelectric point
#'
#' @param seq sequence
#'

removePTM <- function(seq){
  seq <- gsub("o", "", seq)
  seq <- gsub("m", "", seq)
  seq <- gsub("n", "", seq)
  seq <- gsub("p", "", seq)
  return (seq)
}


#' processTerminalSequence
#'
#' This function reformat the sequence removing Search Engine notation.
#' 
#' Example:
#' 
#' sequence <- "K.SDFGHQASSR.L"
#' s <- processTerminalSequence(seq = sequence)
#' 
#' The result will be:
#' s = SDFGHQASSR
#' 
#' @param seq sequence
#' 
processTerminalSequence <- function (seq){
  
  before_dot <- 3
  after_dot   <- nchar(seq) - 2
  seq <- substr(seq, start = before_dot, stop = after_dot)
  return (seq)
}