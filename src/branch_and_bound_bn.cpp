// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "CommonFunction.h"

#include <vector>
#include <map>
#include <string>
#include <functional>
#include <cmath>


using namespace Rcpp;
using namespace arma;

// ----------------- Precompute local scores (Gaussian & Multinom) -----------------
// Return a List where each element is a NumericVector of length max_mask (2^p)
// [[Rcpp::export]]
List precompute_local_scores_gaussian_cpp(const arma::mat &X) {
  int p = X.n_cols;
  int max_mask = 1 << p;
  List out(p);
  for (int j = 0; j < p; ++j) {
    NumericVector scores(max_mask, R_NegInf);
    for (int mask = 0; mask < max_mask; ++mask) {
      // skip masks with jth bit set
      if (mask & (1 << j)) {
        scores[mask] = R_NegInf;
        continue;
      }
      // extract parents indices (0-based)
      std::vector<unsigned int> parents_vec;
      for (int b = 0; b < p; ++b) if ((mask & (1 << b)) && b != j) parents_vec.push_back((unsigned int)b);
      arma::uvec parents;
      if (!parents_vec.empty()) parents = arma::uvec(parents_vec);
      else parents.reset();
      double sc = local_bic_gaussian_cpp(j, parents, X);
      scores[mask] = sc;
    }
    out[j] = scores;
  }
  return out;
}

// [[Rcpp::export]]
List precompute_local_scores_multinom_cpp(const arma::imat &X_disc, const arma::ivec &n_levels) {
  int p = X_disc.n_cols;
  int max_mask = 1 << p;
  List out(p);
  for (int j = 0; j < p; ++j) {
    NumericVector scores(max_mask, R_NegInf);
    for (int mask = 0; mask < max_mask; ++mask) {
      if (mask & (1 << j)) {
        scores[mask] = R_NegInf;
        continue;
      }
      std::vector<unsigned int> parents_vec;
      for (int b = 0; b < p; ++b) if ((mask & (1 << b)) && b != j) parents_vec.push_back((unsigned int)b);
      arma::uvec parents;
      if (!parents_vec.empty()) parents = arma::uvec(parents_vec);
      else parents.reset();
      double sc = local_bic_multinom_cpp(j, parents, X_disc, n_levels);
      scores[mask] = sc;
    }
    out[j] = scores;
  }
  return out;
}

// ----------------- best_score_given_predecessors_cpp (submask enumeration) -----------------
// [[Rcpp::export]]
double best_score_given_predecessors_cpp(const NumericVector &scores_j, int preds_mask) {
  double best = R_NegInf;
  int sub = preds_mask;
  while (true) {
    int s = sub;
    double val = scores_j[s];
    if (!R_IsNA(val) && std::isfinite(val) && val > best) best = val;
    if (sub == 0) break;
    sub = (sub - 1) & preds_mask;
  }
  return best;
}

// ----------------- build_adj_from_parent_masks -----------------
// parent_masks: IntegerVector length p
// [[Rcpp::export]]
arma::imat build_adj_from_parent_masks(const IntegerVector &parent_masks) {
  int p = parent_masks.size();
  arma::imat adj(p, p, arma::fill::zeros);
  for (int j = 0; j < p; ++j) {
    int mask = parent_masks[j];
    if (mask == 0) continue;
    for (int b = 0; b < p; ++b) {
      if (mask & (1 << b)) adj(b, j) = 1;
    }
  }
  return adj;
}

// ----------------- Branch and Bound (supports gaussian or multinom) -----------------
//
// For multinomial: provide X_disc (integer-coded 0..r-1) and n_levels
// For gaussian: provide X (numeric matrix)
//
// [[Rcpp::export]]
List branch_and_bound_bn_cpp(const Nullable<RObject> &X_nullable = R_NilValue,
                             const Nullable<arma::imat> &X_disc_nullable = R_NilValue,
                             const Nullable<arma::ivec> &n_levels_nullable = R_NilValue,
                             std::string distribution = "gaussian") {

  bool is_gauss = (distribution == "gaussian");
  bool is_multi = (distribution == "multinomial");
  if (!is_gauss && !is_multi) stop("distribution must be 'gaussian' or 'multinomial'");

  arma::mat X;
  arma::imat X_disc;
  arma::ivec n_levels;

  if (is_gauss) {
    if (X_nullable.isNull()) stop("For gaussian distribution, provide X (numeric matrix).");
    X = as<arma::mat>(X_nullable);
  } else {
    // multinomial: automatically convert X to X_disc and n_levels if not provided
    if (!X_disc_nullable.isNull() && !n_levels_nullable.isNull()) {
      X_disc = as<arma::imat>(X_disc_nullable);
      n_levels = as<arma::ivec>(n_levels_nullable);
    } else if (!X_nullable.isNull()) {
      // Convert R dataframe or matrix to integer matrix
      Rcpp::DataFrame df(X_nullable);
      int N = df.nrows();
      int p_col = df.size();
      X_disc.set_size(N, p_col);
      n_levels.set_size(p_col);

      for (int j = 0; j < p_col; ++j) {
        Rcpp::RObject col = df[j];
        Rcpp::IntegerVector col_int = Rcpp::as<Rcpp::IntegerVector>(col);
        int max_lvl = 0;
        for (int i = 0; i < N; ++i) {
          X_disc(i,j) = col_int[i] - 1; // assumes factors in R start at 1
          if (X_disc(i,j) > max_lvl) max_lvl = X_disc(i,j);
        }
        n_levels[j] = max_lvl + 1;
      }
    } else {
      stop("For multinomial provide X, or X_disc and n_levels.");
    }
  }

  int p = is_gauss ? X.n_cols : X_disc.n_cols;
  if (p <= 1) {
    arma::imat adj(p, p, arma::fill::zeros);
    double s = 0.0;
    if (p == 1) {
      arma::uvec empty; empty.reset();
      if (is_gauss) s = local_bic_gaussian_cpp(0, empty, X);
      else s = local_bic_multinom_cpp(0, empty, X_disc, n_levels);
    }
    return List::create(_["best_adj"] = adj,
                        _["best_score"] = s,
                        _["time"] = 0.0,
                        _["nodes_explored"] = 1,
                        _["pruned"] = 0);
  }

  Rcout << "Precomputing local scores for every node and parent subset (2^p masks each)...\n";
  auto t0 = std::chrono::high_resolution_clock::now();

  List scores_list;
  if (is_gauss) scores_list = precompute_local_scores_gaussian_cpp(X);
  else scores_list = precompute_local_scores_multinom_cpp(X_disc, n_levels);

  int max_mask = 1 << p;

  // best_possible per node
  NumericVector best_possible(p, R_NegInf);
  for (int j = 0; j < p; ++j) {
    NumericVector sj = scores_list[j];
    double m = R_NegInf;
    for (int i = 0; i < sj.size(); ++i) {
      double v = sj[i];
      if (!R_IsNA(v) && std::isfinite(v) && v > m) m = v;
    }
    best_possible[j] = m;
  }

  double global_best_score = R_NegInf;
  IntegerVector global_best_parent_masks(p);
  long nodes_explored = 0;
  long pruned = 0;
  int all_nodes_mask = (1 << p) - 1;

  // recursive DFS lambda
  std::function<void(std::vector<int>&, IntegerVector&, double, int)> dfs =
    [&](std::vector<int>& prefix, IntegerVector& parent_masks_so_far, double current_score, int remaining_nodes_mask) -> void {
      nodes_explored++;

      if (remaining_nodes_mask == 0) {
        if (current_score > global_best_score) {
          global_best_score = current_score;
          global_best_parent_masks = clone(parent_masks_so_far);
        }
        return;
      }

      // compute optimistic bound
      std::vector<int> rem_indices;
      for (int b = 0; b < p; ++b) if (remaining_nodes_mask & (1 << b)) rem_indices.push_back(b + 1);

      double optimistic = current_score;
      for (int idx : rem_indices) optimistic += best_possible[idx - 1];
      if (optimistic <= global_best_score) {
        pruned++;
        return;
      }

      // order_try: rem_indices sorted by descending best_possible
      std::sort(rem_indices.begin(), rem_indices.end(),
                [&](int a, int b){ return best_possible[a-1] > best_possible[b-1]; });

      // preds_mask from prefix
      int preds_mask = 0;
      for (int v : prefix) preds_mask |= (1 << (v - 1));

      for (int j : rem_indices) {
        NumericVector scores_j = scores_list[j - 1];
        double best_local = best_score_given_predecessors_cpp(scores_j, preds_mask);
        if (!std::isfinite(best_local)) continue;

        // find mask achieving best_local (search submasks)
        int chosen_mask = 0;
        int sub = preds_mask;
        bool found = false;
        while (true) {
          int s = sub;
          double val = scores_j[s];
          if (!R_IsNA(val) && std::isfinite(val) && val == best_local) {
            chosen_mask = s;
            found = true;
            break;
          }
          if (sub == 0) break;
          sub = (sub - 1) & preds_mask;
        }

        if (!found) {
          // fallback: iterate all masks
          for (int m = 0; m < max_mask; ++m) {
            if (m & (1 << (j - 1))) continue;
            if ((m & (~preds_mask)) != 0) continue;
            double val = scores_j[m];
            if (!R_IsNA(val) && std::isfinite(val) && val == best_local) {
              chosen_mask = m;
              break;
            }
          }
        }

        IntegerVector parent_masks_new = clone(parent_masks_so_far);
        parent_masks_new[j - 1] = chosen_mask;
        prefix.push_back(j);
        int remaining_new = remaining_nodes_mask & (~(1 << (j - 1)));
        dfs(prefix, parent_masks_new, current_score + best_local, remaining_new);
        prefix.pop_back();
      }
    };

    std::vector<int> start_prefix;
    IntegerVector start_parent_masks(p);
    dfs(start_prefix, start_parent_masks, 0.0, all_nodes_mask);

    auto t1 = std::chrono::high_resolution_clock::now();
    double exec_time = std::chrono::duration<double>(t1 - t0).count();

    arma::imat best_adj = build_adj_from_parent_masks(global_best_parent_masks);

    return List::create(_["best_adj"] = best_adj,
                        _["best_score"] = global_best_score,
                        _["time"] = exec_time,
                        _["nodes_explored"] = nodes_explored,
                        _["pruned"] = pruned);
}
