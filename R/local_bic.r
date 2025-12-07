#Local Gaussian BIC
local_bic_gaussian <- function(j, parents, X) {
  X <- as.matrix(X)
  y <- X[, j]
  N <- nrow(X)

  if (length(parents) == 0L) {
    mu <- mean(y)
    rss <- sum((y - mu)^2)
    sigma2_hat <- rss / N
    loglik <- -0.5 * N * (log(2 * pi * sigma2_hat) + 1)
    k <- 2L
    return(loglik - 0.5 * log(N) * k)
  }

  # with parents
  Xp <- X[, parents, drop = FALSE]
  X_design <- cbind(1, Xp)
  XtX <- crossprod(X_design)
  XtY <- crossprod(X_design, y)

  beta_hat <- solve(XtX, XtY)
  resid <- y - X_design %*% beta_hat
  rss <- sum(resid^2)
  sigma2_hat <- rss / N

  loglik <- -0.5 * N * (log(2 * pi * sigma2_hat) + 1)
  k <- length(parents) + 2L
  loglik - 0.5 * log(N) * k
}

#Local Multinomial BIC
local_bic_multinomial <- function(j, parents, X) {
  df <- as.data.frame(X)
  df[] <- lapply(df, factor)

  y <- df[[j]]
  N <- nrow(df)

  if (length(parents) == 0L) {
    tab <- table(y)
    probs_hat <- tab / sum(tab)
    loglik <- sum(tab * log(probs_hat))
    r_j <- length(tab)
    k <- r_j - 1L
    return(loglik - 0.5 * log(N) * k)
  }

  parent_names <- names(df)[parents]
  tmp <- df[, c(parent_names, names(df)[j]), drop = FALSE]

  tab <- as.data.frame(table(tmp))

  key <- interaction(tab[, parent_names], drop = TRUE)
  N_pa <- tapply(tab$Freq, key, sum)
  N_pa_each <- N_pa[key]

  # ---- LAPLACE SMOOTHING ----
  alpha <- 1
  r_j <- nlevels(y)
  N_pa_each_smoothed <- N_pa_each + alpha * r_j
  probs_hat <- (tab$Freq + alpha) / N_pa_each_smoothed

  loglik <- sum(tab$Freq * log(probs_hat))

  q_j <- length(N_pa)
  k <- (r_j - 1L) * q_j

  loglik - 0.5 * log(N) * k
}

#Global BIC Score
bic_score_bn <- function(adj, X, distribution = c("gaussian", "multinomial")) {
  distribution <- match.arg(distribution)

  if (!is_acyclic(adj)) return(-Inf)

  p <- ncol(X)
  total <- 0

  for (j in seq_len(p)) {
    parents <- which(adj[, j] == 1L)

    if (distribution == "gaussian") {
      total <- total + local_bic_gaussian(j, parents, X)
    } else {
      total <- total + local_bic_multinomial(j, parents, X)
    }
  }

  total
}
