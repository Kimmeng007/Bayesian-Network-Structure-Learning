#Check whether an adjacency matrix represents an acyclic graph
is_acyclic <- function(adj) {
  p <- nrow(adj)
  visited <- rep(FALSE, p)
  rec     <- rep(FALSE, p)

  has_cycle <- function(v) {
    visited[v] <<- TRUE
    rec[v] <<- TRUE

    for (u in which(adj[v, ] == 1)) {
      if (!visited[u] && has_cycle(u)) return(TRUE)
      if (rec[u]) return(TRUE)
    }

    rec[v] <<- FALSE
    FALSE
  }

  !any(sapply(1:p, has_cycle))
}

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
