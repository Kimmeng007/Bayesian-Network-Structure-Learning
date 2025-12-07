#naive structure learning
naive_bn_structure_learning <- function(data, distribution = c("gaussian","multinomial")) {
  distribution <- match.arg(distribution)
  p <- ncol(data)

  cat("Generating all possible DAGs...\n")
  start_time <- Sys.time()  # Start timing
  dags <- generate_all_dags(p)
  cat(length(dags), "acyclic DAGs generated.\n")

  best_score <- -Inf
  best_adj <- NULL

  for (adj in dags) {
    sc <- bic_score_bn(adj, data, distribution)
    if (sc > best_score) {
      best_score <- sc
      best_adj <- adj
    }
  }
  end_time <- Sys.time()  # End timing
  exec_time <- end_time - start_time

  cat("Execution time:", exec_time, "\n")
  list(best_adj = best_adj, best_score = best_score, time = exec_time)
}
