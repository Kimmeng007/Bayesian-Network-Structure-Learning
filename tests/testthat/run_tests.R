# run_tests.R - Test runner for Bayesian Network functions
library(testthat)

source("R/naive.R")

cat("Running Bayesian Network Structure Learning tests...\n")
cat("================================================\n\n")


# Option 2: Or run all test files in the directory
test_results <- test_dir(".", pattern = "^test_.*\\.R$", reporter = "progress")

# Print detailed results
print(test_results)

cat("\n================================================\n")
cat("Tests completed.\n")
if (all(test_results$passed)) {
  cat("✓ All tests passed!\n")
} else {
  cat("✗ Some tests failed.\n")
}
