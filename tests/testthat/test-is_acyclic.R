test_that("is_acyclic correctly identifies DAGs", {
  # Acyclic graph
  adj_acyclic <- matrix(c(0, 0, 0,
                          1, 0, 0,
                          0, 1, 0), nrow = 3, byrow = TRUE)
  expect_true(is_acyclic(adj_acyclic))

  # Cyclic graph
  adj_cyclic <- matrix(c(0, 1, 0,
                         0, 0, 1,
                         1, 0, 0), nrow = 3, byrow = TRUE)
  expect_false(is_acyclic(adj_cyclic))
})
