library(testthat)

test_that("bic_score_bn returns numeric for Gaussian", {
  data(iris)
  X <- iris[, 1:4]
  adj <- matrix(0, ncol(X), ncol(X))  # empty DAG
  score <- bic_score_bn(adj, X, "gaussian")
  expect_true(is.numeric(score))
})

test_that("bic_score_bn returns numeric for Multinomial", {
  Titanic_df <- as.data.frame(Titanic)
  Titanic_full <- Titanic_df[rep(seq_len(nrow(Titanic_df)), Titanic_df$Freq), 1:4]
  adj <- matrix(0, ncol(Titanic_full), ncol(Titanic_full))
  score <- bic_score_bn(adj, Titanic_full, "multinomial")
  expect_true(is.numeric(score))
})

