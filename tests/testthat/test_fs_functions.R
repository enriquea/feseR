
library(feseR)

context("Feature selection evals")

# read original data
data(TNBC)

# retrive features
features <- TNBC[,-ncol(TNBC)]

# retrive class variables (expected last column)
class <- TNBC[,ncol(TNBC)]

# getting only those features (i.e. gene/protein expression values) with values for all instances (i.e. samples)
features <- features[ , colSums(is.na(features)) == 0]

# Scale data features
# These transformations coerce the original predictors to have zero mean and standard deviation equal one.
features <- scale(features, center=TRUE, scale=TRUE)


test_that("Test filter correlation function", {
  output <- filter.corr(features = features, class = class, mincorr = 0.3)
  expect_equal(ncol(output), 1944)
})

test_that("Test Gain Information filter function", {
  output <- filter.gain.inf(features = features, class = class, zero.gain.out = TRUE)
  expect_equal(ncol(output), 615)
})


test_that("Test Matrix correlation filter function", {
  output <- filter.matrix.corr(features = features, maxcorr = 0.75)
  expect_equal(ncol(output), 1531)
})


test_that("Test pca function for data reduction", {
  output <- filter.pca(features = features, cum.var.cutoff = .95)
  expect_equal(ncol(output), 31)
})

test_that("Test combineFS function (FS workflow: X2-MC-RFE-RF)", {
  results <- combineFS(features = features, class = class,
                       univariate = 'corr', mincorr = 0.3,
                       multivariate = 'mcorr', maxcorr = .75,
                       wrapper = 'rfe.rf',
                       number.cv = 10,
                       extfolds = 10,
                       verbose = FALSE)

  max.acc.testing <- max(as.data.frame(results$testing)$Accuracy)
  expect_equal(max.acc.testing, 1)
})

test_that("Test combineFS function (FS workflow: GI-PCA-RFE-RF)", {
  results <- combineFS(features = features, class = class,
                       univariate = 'gain', zero.gain.out = TRUE,
                       multivariate = 'pca', cum.var.cutoff = 1,
                       wrapper = 'rfe.rf',
                       number.cv = 10,
                       extfolds = 10,
                       verbose = FALSE)

  max.acc.testing <- max(as.data.frame(results$testing)$Accuracy)
  expect_more_than(max.acc.testing, 0.90)

})

