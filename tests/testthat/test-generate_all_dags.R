library(testthat)

test_that("generate_all_dags produces only acyclic DAGs", {
  dags <- generate_all_dags(3)
  expect_true(length(dags) > 0)
  for (adj in dags) {
    expect_true(is_acyclic(adj))
  }
})
