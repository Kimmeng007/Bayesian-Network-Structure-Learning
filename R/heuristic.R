############################################################
# heuristic.R
############################################################

#------------------ R IMPLEMENTATIONS ----------------------#

#' Greedy hill-climbing structure learning for Bayesian networks (R version)
#'
#' @param X Data frame or matrix of size \code{N x p} containing the data.
#'   Must be numeric for \code{distribution = "gaussian"} and categorical
#'   (or coercible to factors) for \code{distribution = "multinomial"}.
#' @param distribution Character string, either \code{"gaussian"} or
#'   \code{"multinomial"}.
#' @param max_iter Integer, maximum number of hill-climbing iterations.
#' @param verbose Logical; if \code{TRUE}, prints progress information.
#'
#' @return A list with components:
#'   \item{adj}{Final adjacency matrix (p x p) of the learned DAG.}
#'   \item{score}{BIC score of the final DAG.}
#'   \item{history}{Numeric vector of BIC scores across iterations.}
#'   \item{iterations}{Number of iterations actually performed.}
#'   \item{distribution}{The distribution assumption used.}
#'
#' @export
hill_climbing_bn <- function(X,
                             distribution = c("gaussian", "multinomial"),
                             max_iter = 100,
                             verbose = FALSE) {

  distribution <- match.arg(distribution)
  X <- as.data.frame(X)
  p <- ncol(X)

  # Initialize with empty graph
  adj <- matrix(0L, nrow = p, ncol = p)
  colnames(adj) <- rownames(adj) <- colnames(X)

  best_score <- bic_score_bn(adj, X, distribution)
  score_history <- c(best_score)

  if (verbose) {
    cat("Hill-Climbing (R): initial score =", best_score, "\n")
  }

  for (iter in seq_len(max_iter)) {
    best_move <- NULL
    best_move_score <- best_score

    # Explore all possible local moves (add, remove, reverse)
    for (i in seq_len(p)) {
      for (j in seq_len(p)) {
        if (i == j) next

        # ADD edge i -> j
        if (adj[i, j] == 0L) {
          new_adj <- adj
          new_adj[i, j] <- 1L
          if (is_acyclic(new_adj)) {
            sc <- bic_score_bn(new_adj, X, distribution)
            if (sc > best_move_score) {
              best_move_score <- sc
              best_move <- list(type = "add", from = i, to = j, adj = new_adj)
            }
          }

        } else {
          # REMOVE edge i -> j
          new_adj <- adj
          new_adj[i, j] <- 0L
          sc <- bic_score_bn(new_adj, X, distribution)
          if (sc > best_move_score) {
            best_move_score <- sc
            best_move <- list(type = "remove", from = i, to = j, adj = new_adj)
          }

          # REVERSE edge i -> j to j -> i
          if (adj[j, i] == 0L) {
            new_adj <- adj
            new_adj[i, j] <- 0L
            new_adj[j, i] <- 1L
            if (is_acyclic(new_adj)) {
              sc <- bic_score_bn(new_adj, X, distribution)
              if (sc > best_move_score) {
                best_move_score <- sc
                best_move <- list(type = "reverse", from = i, to = j, adj = new_adj)
              }
            }
          }
        }
      }
    }

    # If no improving move, stop
    if (is.null(best_move) || best_move_score <= best_score) {
      if (verbose) {
        cat("Hill-Climbing (R): no improvement at iteration", iter, "\n")
      }
      break
    }

    # Apply best move
    adj <- best_move$adj
    best_score <- best_move_score
    score_history <- c(score_history, best_score)

    if (verbose) {
      cat(
        "Hill-Climbing (R): iter", iter,
        "- move", best_move$type, best_move$from, "->", best_move$to,
        "- score =", best_score, "\n"
      )
    }
  }

  list(
    adj = adj,
    score = best_score,
    history = score_history,
    iterations = length(score_history) - 1L,
    distribution = distribution
  )
}


#' Tabu search structure learning for Bayesian networks (R version)
#'
#' @param X Data frame or matrix of size \code{N x p} containing the data.
#' @param distribution "gaussian" or "multinomial".
#' @param max_iter Integer, maximum number of iterations.
#' @param tabu_tenure Integer, number of iterations a move remains tabu.
#' @param verbose Logical; if \code{TRUE}, prints progress information.
#'
#' @return A list with components:
#'   \item{adj}{Final adjacency matrix of the last visited DAG.}
#'   \item{score}{BIC score of the final DAG.}
#'   \item{best_adj}{Adjacency matrix of the best-scoring DAG encountered.}
#'   \item{best_score}{Best BIC score encountered.}
#'   \item{history}{Scores across iterations.}
#'   \item{iterations}{Number of iterations performed.}
#'   \item{distribution}{Distribution used.}
#'   \item{tabu_tenure}{Tabu tenure.}
#'
#' @export
tabu_search_bn <- function(X,
                           distribution = c("gaussian", "multinomial"),
                           max_iter = 200,
                           tabu_tenure = 10,
                           verbose = FALSE) {

  distribution <- match.arg(distribution)
  X <- as.data.frame(X)
  p <- ncol(X)

  # Initialize with empty graph
  adj <- matrix(0L, nrow = p, ncol = p)
  colnames(adj) <- rownames(adj) <- colnames(X)

  current_score <- bic_score_bn(adj, X, distribution)
  best_score <- current_score
  best_adj <- adj

  score_history <- c(current_score)

  if (verbose) {
    cat("Tabu (R): initial score =", current_score, "\n")
  }

  tabu_list <- character(0)

  for (iter in seq_len(max_iter)) {
    best_candidate <- NULL
    best_candidate_score <- -Inf

    for (i in seq_len(p)) {
      for (j in seq_len(p)) {
        if (i == j) next
        edge_key <- paste(i, j, sep = "_")

        # ADD
        if (adj[i, j] == 0L) {
          new_adj <- adj
          new_adj[i, j] <- 1L
          if (is_acyclic(new_adj)) {
            sc <- bic_score_bn(new_adj, X, distribution)
            is_tabu <- edge_key %in% tabu_list
            is_asp <- sc > best_score
            if (!is_tabu || is_asp) {
              if (sc > best_candidate_score) {
                best_candidate_score <- sc
                best_candidate <- list(
                  type = "add", from = i, to = j,
                  adj = new_adj, edge_key = edge_key
                )
              }
            }
          }

        } else {
          # REMOVE
          new_adj <- adj
          new_adj[i, j] <- 0L
          sc <- bic_score_bn(new_adj, X, distribution)
          is_tabu <- edge_key %in% tabu_list
          is_asp <- sc > best_score
          if (!is_tabu || is_asp) {
            if (sc > best_candidate_score) {
              best_candidate_score <- sc
              best_candidate <- list(
                type = "remove", from = i, to = j,
                adj = new_adj, edge_key = edge_key
              )
            }
          }

          # REVERSE i -> j to j -> i
          if (adj[j, i] == 0L) {
            new_adj <- adj
            new_adj[i, j] <- 0L
            new_adj[j, i] <- 1L
            if (is_acyclic(new_adj)) {
              sc <- bic_score_bn(new_adj, X, distribution)
              edge_key_rev <- paste(j, i, sep = "_")
              is_tabu <- edge_key_rev %in% tabu_list
              is_asp <- sc > best_score
              if (!is_tabu || is_asp) {
                if (sc > best_candidate_score) {
                  best_candidate_score <- sc
                  best_candidate <- list(
                    type = "reverse", from = i, to = j,
                    adj = new_adj, edge_key = edge_key_rev
                  )
                }
              }
            }
          }
        }
      }
    }

    if (is.null(best_candidate)) {
      if (verbose) cat("Tabu (R): no admissible move at iteration", iter, "\n")
      break
    }

    # Apply move
    adj <- best_candidate$adj
    current_score <- best_candidate_score
    score_history <- c(score_history, current_score)

    if (current_score > best_score) {
      best_score <- current_score
      best_adj <- adj
    }

    tabu_list <- c(tabu_list, best_candidate$edge_key)
    if (length(tabu_list) > tabu_tenure) {
      tabu_list <- tail(tabu_list, tabu_tenure)
    }

    if (verbose) {
      cat(
        "Tabu (R): iter", iter,
        "- move", best_candidate$type, best_candidate$from, "->", best_candidate$to,
        "- score =", current_score,
        "- best_score =", best_score, "\n"
      )
    }
  }

  list(
    adj = adj,
    score = current_score,
    best_adj = best_adj,
    best_score = best_score,
    history = score_history,
    iterations = length(score_history) - 1L,
    distribution = distribution,
    tabu_tenure = tabu_tenure
  )
}

#------------------ C++ WRAPPERS ---------------------------#

.prepare_discrete_cpp <- function(X) {
  df <- as.data.frame(X)
  facs <- lapply(df, factor)
  X_int <- sapply(facs, function(f) as.integer(f) - 1L)
  X_int <- as.matrix(X_int)
  n_levels <- vapply(facs, nlevels, integer(1))
  list(X_int = X_int, n_levels = n_levels)
}

#' Hill-climbing using C++ (Gaussian or multinomial)
#'
#' @param X Data frame or matrix. Numeric for \code{distribution = "gaussian"};
#'   categorical / factors for \code{distribution = "multinomial"}.
#' @param distribution "gaussian" or "multinomial".
#' @param max_iter Maximum number of iterations.
#' @param verbose Logical; print progress if \code{TRUE}.
#'
#' @return List with fields \code{adj}, \code{score}, \code{history},
#'   \code{iterations}, and \code{distribution}.
#'
#' @export
hill_climbing_bn_cpp <- function(X,
                                 distribution = c("gaussian", "multinomial"),
                                 max_iter = 100,
                                 verbose = FALSE) {
  distribution <- match.arg(distribution)
  p <- ncol(X)

  if (distribution == "gaussian") {
    Xc <- as.matrix(X)
    if (!is.numeric(Xc)) {
      stop("hill_climbing_bn_cpp: X must be numeric for distribution = 'gaussian'.")
    }
    X_disc <- matrix(0L, nrow = nrow(Xc), ncol = ncol(Xc))
    n_levels <- rep(1L, ncol(Xc))

    res <- hill_climbing_bn_cpp_internal(
      X = Xc,
      X_disc = X_disc,
      n_levels = n_levels,
      max_iter = as.integer(max_iter),
      verbose = verbose,
      score_type = 0L
    )
    res$distribution <- "gaussian"
    return(res)
  }

  # multinomial
  prep <- .prepare_discrete_cpp(X)
  X_disc <- prep$X_int
  n_levels <- prep$n_levels
  Xc <- matrix(0.0, nrow = nrow(X_disc), ncol = ncol(X_disc))

  res <- hill_climbing_bn_cpp_internal(
    X = Xc,
    X_disc = X_disc,
    n_levels = n_levels,
    max_iter = as.integer(max_iter),
    verbose = verbose,
    score_type = 1L
  )
  res$distribution <- "multinomial"
  res
}

#' Tabu search using C++ (Gaussian or multinomial)
#'
#' @param X Data frame or matrix. Numeric for \code{distribution = "gaussian"};
#'   categorical / factors for \code{distribution = "multinomial"}.
#' @param distribution "gaussian" or "multinomial".
#' @param max_iter Maximum number of iterations.
#' @param tabu_tenure Integer, tabu tenure (length of tabu list).
#' @param verbose Logical; print progress if \code{TRUE}.
#'
#' @return List with fields \code{adj}, \code{score}, \code{best_adj},
#'   \code{best_score}, \code{history}, \code{iterations}, and \code{distribution}.
#'
#' @export
tabu_search_bn_cpp <- function(X,
                               distribution = c("gaussian", "multinomial"),
                               max_iter = 200,
                               tabu_tenure = 10,
                               verbose = FALSE) {
  distribution <- match.arg(distribution)
  p <- ncol(X)

  if (distribution == "gaussian") {
    Xc <- as.matrix(X)
    if (!is.numeric(Xc)) {
      stop("tabu_search_bn_cpp: X must be numeric for distribution = 'gaussian'.")
    }
    X_disc <- matrix(0L, nrow = nrow(Xc), ncol = ncol(Xc))
    n_levels <- rep(1L, ncol(Xc))

    res <- tabu_search_bn_cpp_internal(
      X = Xc,
      X_disc = X_disc,
      n_levels = n_levels,
      max_iter = as.integer(max_iter),
      tabu_tenure = as.integer(tabu_tenure),
      verbose = verbose,
      score_type = 0L
    )
    res$distribution <- "gaussian"
    return(res)
  }

  # multinomial
  prep <- .prepare_discrete_cpp(X)
  X_disc <- prep$X_int
  n_levels <- prep$n_levels
  Xc <- matrix(0.0, nrow = nrow(X_disc), ncol = ncol(X_disc))

  res <- tabu_search_bn_cpp_internal(
    X = Xc,
    X_disc = X_disc,
    n_levels = n_levels,
    max_iter = as.integer(max_iter),
    tabu_tenure = as.integer(tabu_tenure),
    verbose = verbose,
    score_type = 1L
  )
  res$distribution <- "multinomial"
  res
}
