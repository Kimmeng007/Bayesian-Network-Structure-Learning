library(testthat)

test_that("naive_bn_structure_learning works for Gaussian", {
  data(iris)
  X <- iris[, 1:4]
  result <- naive_bn_structure_learning(X, "gaussian")
  expect_true(is.matrix(result$best_adj))
  expect_equal(ncol(result$best_adj), ncol(X))
  expect_true(is.numeric(result$best_score))
})

test_that("naive_bn_structure_learning works for Multinomial", {
  Titanic_df <- as.data.frame(Titanic)
  Titanic_full <- Titanic_df[rep(seq_len(nrow(Titanic_df)), Titanic_df$Freq), 1:4]
  result <- naive_bn_structure_learning(Titanic_full, "multinomial")
  expect_true(is.matrix(result$best_adj))
  expect_equal(ncol(result$best_adj), ncol(Titanic_full))
  expect_true(is.numeric(result$best_score))
})
#To run all tests
#testthat::test_dir("tests/testthat")
