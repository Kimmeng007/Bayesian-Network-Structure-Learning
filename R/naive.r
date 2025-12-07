# Check if an adjacency matrix is acyclic
is_acyclic <- function(adj) {
  p <- nrow(adj)
  visited  <- integer(p)
  in_stack <- integer(p)

  dfs <- function(u) {
    visited[u]  <<- 1L
    in_stack[u] <<- 1L
    children <- which(adj[u, ] != 0L)

    for (v in children) {
      if (!visited[v] && dfs(v)) return(TRUE)
      if (in_stack[v]) return(TRUE)
    }

    in_stack[u] <<- 0L
    FALSE
  }

  for (u in seq_len(p)) {
    if (!visited[u] && dfs(u)) return(FALSE)
  }

  TRUE
}


# Local BIC score for a Gaussian Bayesian network node
local_bic_gaussian <- function(j, parents, X) {
  X <- as.matrix(X)
  y <- X[, j]
  N <- nrow(X)

  if (length(parents) == 0L) {
    mu <- mean(y)
    rss <- sum((y - mu)^2)
    sigma2_hat <- rss / N
    loglik <- -0.5 * N * (log(2 * pi * sigma2_hat) + 1)
    k <- 2L  # intercept + variance
    return(loglik - 0.5 * log(N) * k)
  }

  Xp <- X[, parents, drop = FALSE]
  X_design <- cbind(1, Xp)

  XtX <- crossprod(X_design)
  XtY <- crossprod(X_design, y)
  beta_hat <- solve(XtX, XtY)

  resid <- y - X_design %*% beta_hat
  rss <- sum(resid^2)
  sigma2_hat <- rss / N

  loglik <- -0.5 * N * (log(2 * pi * sigma2_hat) + 1)
  k <- length(parents) + 2L  # betas + intercept + variance
  loglik - 0.5 * log(N) * k
}


# Local BIC score for a multinomial (discrete) Bayesian network node
local_bic_multinomial <- function(j, parents, X) {
  # X: data.frame or matrix, assumed no NA
  df <- as.data.frame(X, stringsAsFactors = FALSE)

  df[[j]] <- as.factor(df[[j]])
  if (length(parents) > 0L) {
    for (p in parents) df[[p]] <- as.factor(df[[p]])
  }

  N <- nrow(df)
  if (N == 0L) return(-Inf)

  y <- df[[j]]

  # no parents case
  if (length(parents) == 0L) {
    tab <- table(y)
    pos <- which(tab > 0)
    if (length(pos) == 0L) return(-Inf)
    probs_hat <- tab[pos] / sum(tab[pos])
    loglik <- sum(tab[pos] * log(probs_hat))
    r_j <- length(tab)
    k <- r_j - 1L
    return(as.numeric(loglik - 0.5 * log(N) * k))
  }

  # parents present: build contingency table
  parent_names <- names(df)[parents]
  tmp <- df[, c(parent_names, names(df)[j]), drop = FALSE]
  tab_df <- as.data.frame(table(tmp), stringsAsFactors = FALSE)
  names(tab_df)[ncol(tab_df)] <- "Freq"

  # key for parent configurations (interaction), drop unused levels
  key <- do.call("interaction", c(tab_df[parent_names], drop = TRUE))
  if (all(is.na(key))) return(-Inf)

  # compute N_pa for each parent config
  N_pa <- tapply(tab_df$Freq, key, sum)
  keys_unique <- names(N_pa)

  loglik <- 0
  for (kname in keys_unique) {
    rows_k <- which(key == kname)
    freqs <- tab_df$Freq[rows_k]
    pos <- which(freqs > 0)
    if (length(pos) == 0) next
    probs_hat <- freqs[pos] / sum(freqs[pos])
    # safe: probs_hat > 0 and freqs[pos] > 0
    loglik <- loglik + sum(freqs[pos] * log(probs_hat))
  }

  r_j <- nlevels(y)
  q_j <- length(N_pa)
  k <- (r_j - 1L) * q_j
  return(as.numeric(loglik - 0.5 * log(N) * k))
}

# Global BIC score for a Bayesian network
bic_score_bn <- function(adj, X, distribution = c("gaussian", "multinomial")) {
  distribution <- match.arg(distribution)

  # explicit NA check (user requested): fail fast if any NA present
  if (anyNA(X)) stop("bic_score_bn: data X contains NA. Please handle or impute before calling.")

  if (!is_acyclic(adj)) return(-Inf)

  p <- ncol(X)
  total <- 0

  for (j in seq_len(p)) {
    parents <- which(adj[, j] != 0L)

    if (distribution == "gaussian") {
      total <- total + local_bic_gaussian(j, parents, X)
    } else if (distribution == "multinomial") {
      total <- total + local_bic_multinomial(j, parents, X)
    } else {
      stop("Error: Please enter one of the distribution assumption for this model : ('gaussian' or 'multinomial')")
    }
  }

  total
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
