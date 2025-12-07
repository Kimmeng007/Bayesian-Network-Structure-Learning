source("R/CommonFunction.R")

#generate all DAGs
generate_all_dags <- function(p) {
  n_edges <- p * (p - 1)
  all_bin <- expand.grid(rep(list(c(0,1)), n_edges))

  dags <- list()
  idx <- 1

  for (i in 1:nrow(all_bin)) {
    adj <- matrix(0, p, p)
    v <- all_bin[i, ]
    adj[t(matrix(1:p, p, p)) != matrix(1:p, p, p)] <- as.numeric(v)

    if (is_acyclic(adj)) {
      dags[[idx]] <- adj
      idx <- idx + 1
    }
  }

  dags
}

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
