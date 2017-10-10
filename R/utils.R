
## Useful functions for checking/comparing matrices

## valid matrix

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
    mx
  }
}


## compare matrix

compare.matrix <- function(input.matrix, output.matrix, description = 'Number of removed features:') {

  if(!is.matrix(input.matrix) & !is.matrix(output.matrix)) {
    stop('Expected a matrix as input...')
  }

  ncol1 <- ncol(input.matrix)
  ncol2 <- ncol(output.matrix)

  if (ncol2 < 2){
    stop('Removed all features from input matrix...')
  } else {
    message(paste(description, round(ncol2/ncol1*100, 2), '%', sep = ' '))
  }
}
