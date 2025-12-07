# Branch and Bound for Bayesian Network structure learning (discrete or Gaussian)
# Uses functions from CommonFunction.R (is_acyclic, local_bic_multinomial, local_bic_gaussian)
source("R/CommonFunction.R")

# Helper: generate all subsets of a set of indices represented as bitmasks 0..(2^m-1)
# but we'll generate only as needed.

# Precompute local scores for every node and every possible parent subset (excluding the node itself)
precompute_local_scores <- function(X, distribution = c("gaussian","multinomial")) {
  distribution <- match.arg(distribution)
  p <- ncol(X)
  nodes <- seq_len(p)

  # For each node j, we'll map a bitmask over p bits (0-based) to a score
  # but we will ignore the bit for the node itself (force it to 0).
  max_mask <- bitwShiftL(1L, p)
  scores_list <- vector("list", p)

  for (j in nodes) {
    scores <- numeric(max_mask)
    # iterate masks 0..2^p -1 but skip those with jth bit set
    for (mask in 0:(max_mask - 1L)) {
      if (bitwAnd(mask, bitwShiftL(1L, j - 1L)) != 0L) {
        scores[mask + 1L] <- -Inf
        next
      }
      # extract parent indices from mask
      if (mask == 0L) {
        parents <- integer(0)
      } else {
        parents_idx <- which(as.logical(intToBits(mask)[1:p]))
        parents <- parents_idx[parents_idx != j]
      }
      # call appropriate local score
      if (distribution == "gaussian") {
        scores[mask + 1L] <- local_bic_gaussian(j, parents, X)
      } else {
        scores[mask + 1L] <- local_bic_multinomial(j, parents, X)
      }
    }
    scores_list[[j]] <- scores
  }

  scores_list
}

# Convenience: convert set of indices to bitmask (1-based indices)
indices_to_mask <- function(indices, p) {
  if (length(indices) == 0L) return(0L)
  mask <- 0L
  for (i in indices) mask <- bitwOr(mask, bitwShiftL(1L, i - 1L))
  mask
}

# Given available predecessors (as a mask), compute the best local score for node j
best_score_given_predecessors <- function(scores_j, preds_mask) {
  # we need to consider only parent masks that are subsets of preds_mask and have 0 in node bit
  # iterate submasks of preds_mask
  best <- -Inf
  sub <- preds_mask
  repeat {
    s <- sub
    val <- scores_j[s + 1L]
    if (!is.infinite(val) && !is.na(val) && val > best) best <- val
    if (sub == 0L) break
    sub <- bitwAnd(sub - 1L, preds_mask)
  }
  best
}

# Build adjacency matrix from chosen parent masks per node
build_adj_from_parent_masks <- function(parent_masks, p) {
  adj <- matrix(0L, p, p)
  for (j in seq_len(p)) {
    mask <- parent_masks[j]
    if (mask == 0L) next
    parents <- which(as.logical(intToBits(mask)[1:p]))
    for (i in parents) adj[i, j] <- 1L
  }
  adj
}

# Branch and Bound over topological orders. We maintain a prefix (ordered nodes placed)
# For each next node appended, allowed parents are the nodes already placed (the prefix).
# For pruning we use an admissible upper bound: current_score + sum(best_possible[j] for remaining nodes),
# where best_possible[j] is the maximum local score for node j over all parent sets (precomputed).

branch_and_bound_bn <- function(X, distribution = c("gaussian","multinomial")) {
  distribution <- match.arg(distribution)
  p <- ncol(X)
  if (p <= 1L) {
    adj <- matrix(0L, p, p)
    return(list(best_adj = adj, best_score = bic_score_bn(adj, X, distribution), time = 0, nodes_explored = 1L, pruned = 0L))
  }

  cat("Precomputing local scores for every node and parent subset (2^p masks each)...\n")
  start_time <- Sys.time()
  scores_list <- precompute_local_scores(X, distribution)
  max_mask <- bitwShiftL(1L, p)

  # best possible per node over all masks (admissible optimistic bound)
  best_possible <- numeric(p)
  for (j in seq_len(p)) best_possible[j] <- max(scores_list[[j]], na.rm = TRUE)

  global_best_score <- -Inf
  global_best_parent_masks <- rep(0L, p)
  nodes_explored <- 0L
  pruned <- 0L

  # DFS recursion: prefix is vector of placed nodes in order,
  # parent_masks_so_far stores the chosen parent mask for each node once placed (0 for unplaced)
  dfs <- function(prefix, parent_masks_so_far, current_score, remaining_nodes_mask) {
    # remaining_nodes_mask is bitmask of nodes not yet placed
    nodes_explored <<- nodes_explored + 1L

    # compute optimistic upper bound
    if (remaining_nodes_mask == 0L) {
      # complete order: evaluate
      if (current_score > global_best_score) {
        global_best_score <<- current_score
        global_best_parent_masks <<- parent_masks_so_far
      }
      return()
    }

    # optimistic bound = current_score + sum(best_possible for remaining nodes)
    # We can tighten by using for the immediate children allowed preds = mask(prefix)
    # but this simple bound is admissible
    rem_indices <- which(as.logical(intToBits(remaining_nodes_mask)[1:p]))
    optimistic <- current_score + sum(best_possible[rem_indices])
    if (optimistic <= global_best_score) {
      pruned <<- pruned + 1L
      return()
    }

    # enumerate choices for next node: any remaining node
    # better heuristic: try nodes in order of their best_possible descending to reach good solutions earlier
    order_try <- rem_indices[order(-best_possible[rem_indices])]

    # compute mask of predecessors (nodes already placed)
    preds_mask <- 0L
    if (length(prefix) > 0L) preds_mask <- indices_to_mask(prefix, p)

    for (j in order_try) {
      # for chosen next j, allowed parent sets are subsets of preds_mask
      best_local <- best_score_given_predecessors(scores_list[[j]], preds_mask)
      # if there is no valid local score (shouldn't happen), skip
      if (is.infinite(best_local) || is.na(best_local)) next

      # to find which parent mask achieved that best_local, search subsets
      chosen_mask <- 0L
      sub <- preds_mask
      found <- FALSE
      repeat {
        s <- sub
        if (!is.infinite(scores_list[[j]][s + 1L]) && !is.na(scores_list[[j]][s + 1L]) && scores_list[[j]][s + 1L] == best_local) {
          chosen_mask <- s
          found <- TRUE
          break
        }
        if (sub == 0L) break
        sub <- bitwAnd(sub - 1L, preds_mask)
      }
      if (!found) {
        # fallback: iterate all masks (rare)
        for (m in 0:(max_mask - 1L)) {
          if (bitwAnd(m, bitwShiftL(1L, j - 1L)) != 0L) next
          if (bitwAnd(m, bitwNot(preds_mask)) != 0L) next # m must be subset of preds_mask
          if (!is.infinite(scores_list[[j]][m + 1L]) && !is.na(scores_list[[j]][m + 1L]) && scores_list[[j]][m + 1L] == best_local) {
            chosen_mask <- m
            break
          }
        }
      }

      # update and recurse
      parent_masks_so_far_new <- parent_masks_so_far
      parent_masks_so_far_new[j] <- chosen_mask
      prefix_new <- c(prefix, j)
      remaining_new <- bitwAnd(remaining_nodes_mask, bitwNot(bitwShiftL(1L, j - 1L)))

      dfs(prefix_new, parent_masks_so_far_new, current_score + best_local, remaining_new)
    }
  }

  # start recursion with empty prefix
  all_nodes_mask <- bitwShiftL(1L, p) - 1L
  dfs(integer(0), rep(0L, p), 0, all_nodes_mask)

  end_time <- Sys.time()
  exec_time <- end_time - start_time

  best_adj <- build_adj_from_parent_masks(global_best_parent_masks, p)

  list(best_adj = best_adj,
       best_score = global_best_score,
       time = exec_time,
       nodes_explored = nodes_explored,
       pruned = pruned)
}


