
## Utils functions

#' valid.matrix
#'
#' Check if input is a valid matrix
#'
#' @param mx A numeric matrix as input.
#'
#' @return A matrix.
#'

valid.matrix <- function(mx) {
  if(!is.matrix(mx)) {
    stop('Expected a matrix as input...')
  } else if (any(is.na(mx))) {
    stop('Matrix with some values missing...')
  } else if (length(rownames(mx)) == 0) {
    stop('Row names are missing...')
  } else if (length(colnames(mx)) == 0) {
    stop('Column names are missing...')
  } else {
    return(mx)
  }
}



#' compare.matrix 
#'
#' @param input.matrix A numeric matrix as input.
#' @param output.matrix A numeric matrix as input.
#'
#' @return A message object.
#'

compare.matrix <- function(input.matrix, output.matrix) {

  if(!is.matrix(input.matrix) & !is.matrix(output.matrix)) {
    stop('Expected a matrix as input...')
  }

  ncol1 <- ncol(input.matrix)
  ncol2 <- ncol(output.matrix)

  if (ncol2 < 2){
    stop('Removed all features from input matrix...')
  } else {
    pct_removed = round(((ncol1 - ncol2)/ncol1*100), 2)
    message(paste0(
      'Kept ', ncol2, ' features out of ', ncol1, '...',
      ' Number of removed features: ', pct_removed, '%'))
  }
}


#' computeMissingRate
#'
#' @param v Input vector
#'
#' @return Missing rate value
#' @export
#'
#' @examples
computeMissingRate <- function(v){
   v <- as.numeric(v)
   mr <- sum(is.na(v)) / length(v)
   return(round(mr, 2))
}



#' getMissingnessRateByFeatures
#'
#' Given an input matrix, compute the missingness rate for every feature (expected in cols)
#'
#' @param mx A numeric matrix as input.
#'
#' @return A dataframe with features and missingness rate.
#' @export
#'
#' @examples
getMissingnessRateByFeatures <- function(mx){
  if(length(colnames(mx)) == 0) stop('Column (feature) names are missing...')
  missingness_rate <- apply(mx, 2, computeMissingRate)
  df <- data.frame(features=names(missingness_rate),
                   missingness=missingness_rate)
  return(df)
}


#' impute
#'
#' @param v Input numeric vector
#' @param method One of 'mean' of 'median'
#'
#' @return A vector with missing values imputed.
#' @export
#'
#' @examples
impute <- function(v, method='mean'){
   if(!method %in% c('mean', 'median')){
     stop('method must be one of mean or median...')
     }
   v <- as.numeric(v)
   if(method=='mean'){
     v[is.na(v)] <- mean(v, na.rm = T)
   } 
   if(method=='median'){
     v[is.na(v)] <- median(v, na.rm = T)
   }
   return(v)
}


#' imputeMatrix
#' 
#' Impute missing values by computing the mean or median across samples (rows)
#' for every feature (cols).
#'
#' @param mx  A numeric matrix as input.
#' @param method Imputation method for missing values (mean or median).
#'
#' @return Imputed matrix
#' @export
#'
#' @examples
imputeMatrix <- function(mx, method='mean'){
   message(paste0("Imputing matrix with method: ", method))
   mx.imputed <- apply(mx, 2, impute, method = method)
   return(mx.imputed)
}


#' filterMissingnessRate
#'
#' @param mx A numeric matrix as input.
#' @param max_missing_rate Maximal missing rate allow for a feature (default: 0.25).
#'
#' @return Filtered matrix with features (col) passing the max missing rate threshold removed.
#' @export
#'
#' @examples
filterMissingnessRate <- function(mx, max_missing_rate=0.25){
  if(!is.matrix(mx)) {
    stop('Expected a matrix as input...')
  }
  if (length(rownames(mx)) == 0) {
    stop('Row names are missing...')
  } 
  if (length(colnames(mx)) == 0) {
    stop('Column names are missing...')
  } 
  
  df = getMissingnessRateByFeatures(mx)
  
  features_to_keep = subset(df, df[,'missingness'] < 
                              max_missing_rate)[,'features']
  message(
    paste0("Found ", length(features_to_keep),
           ' features with max missing rate ', max_missing_rate))
  
  # keep features passing the missing rate cutoff
  mx.filtered = mx[,colnames(mx) %in% features_to_keep]
  
  compare.matrix(input.matrix = mx,
                 output.matrix = mx.filtered)
  
  return(mx.filtered)
}