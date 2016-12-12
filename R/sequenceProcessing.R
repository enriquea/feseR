
#' mapAA
#'
#' This function renames amino acid from 3-letter to 1-letter cod
#'
#' @param seqvec vector of amino acid (3-letter code)
#'
mapAA <- function(seqvec){
  seqvec <- revalue(seqvec, c("Ala"="A", "Arg"="R", "Asn"="N", "Asp"="D", "Cys"="C", 
                              "Gln"="Q", "Glu"="E", "Gly"="G", "His"="H", "Ile"="I", 
                              "Leu"="L", "Lys"="K", "Met"="M", "Phe"="F", "Pro"="P",
                              "Ser"="S", "Thr"="T", "Trp"="W", "Tyr"="Y", "Val"="V"))
  return(seqvec)
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