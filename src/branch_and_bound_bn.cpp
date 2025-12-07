// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "branch_and_bound_bn.h"  // your header file with declarations

#include <vector>
#include <map>
#include <string>
#include <functional>
#include <cmath>

// ----------------- cycle detection functions (0-based indices) -----------------
bool dfs_cycle(int u, const arma::imat &adj, std::vector<int> &state) {
  state[u] = 1;  // visiting

  int p = adj.n_cols;
  for (int v = 0; v < p; ++v) {
    if (adj(u, v) == 0) continue;
    if (state[v] == 1) return true;
    if (state[v] == 0 && dfs_cycle(v, adj, state)) return true;
  }

  state[u] = 2;
  return false;
}

bool is_acyclic_cpp(const arma::imat &adj) {
  int p = adj.n_cols;
  std::vector<int> state(p, 0);
  for (int u = 0; u < p; ++u) {
    if (state[u] == 0) {
      if (dfs_cycle(u, adj, state)) return false;
    }
  }
  return true;
}

// ----------------- Gaussian local BIC (your function, 0-based indices) -----------------
// [[Rcpp::export]]
double local_bic_gaussian_cpp(int j, const arma::uvec &parents, const arma::mat &X) {
  int N = X.n_rows;
  arma::vec y = X.col(j);

  if (parents.n_elem == 0) {
    double mu = arma::mean(y);
    arma::vec r = y - mu;
    double rss = arma::dot(r, r);
    double sigma2 = rss / N;
    if (sigma2 <= 0) sigma2 = 1e-12;
    double loglik = -0.5 * N * (std::log(2.0 * datum::pi * sigma2) + 1.0);
    int k = 2;
    return loglik - 0.5 * std::log((double)N) * k;
  }

  int q = parents.n_elem;
  arma::mat Xd(N, q + 1);
  Xd.col(0).ones();
  for (int k = 0; k < q; ++k) {
    Xd.col(k + 1) = X.col(parents(k));
  }

  arma::mat XtX = Xd.t() * Xd;
  arma::vec Xty = Xd.t() * y;

  arma::vec beta;
  bool solved = arma::solve(beta, XtX, Xty);
  if (!solved) beta = arma::pinv(XtX) * Xty;

  arma::vec r = y - Xd * beta;
  double rss = arma::dot(r, r);
  double sigma2 = rss / N;
  if (sigma2 <= 0) sigma2 = 1e-12;
  double loglik = -0.5 * N * (std::log(2.0 * datum::pi * sigma2) + 1.0);

  int k = q + 2;
  return loglik - 0.5 * std::log((double)N) * k;
}

// Graph-level Gaussian BIC using your local function
// [[Rcpp::export]]
double bic_score_gaussian_cpp(const arma::imat &adj, const arma::mat &X) {
  if (!is_acyclic_cpp(adj)) return R_NegInf;
  int p = X.n_cols;
  double total = 0.0;
  for (int j = 0; j < p; ++j) {
    arma::uvec parents = arma::find(adj.col(j) != 0);
    total += local_bic_gaussian_cpp(j, parents, X);
  }
  return total;
}

// ----------------- Multinomial local BIC (your function) -----------------
double local_bic_multinom_cpp(int j,
                              const arma::uvec &parents,
                              const arma::imat &X_disc,
                              const arma::ivec &n_levels) {
  int N = X_disc.n_rows;
  int r_j = n_levels(j);

  if (parents.n_elem == 0) {
    std::vector<int> counts(r_j, 0);
    for (int n = 0; n < N; ++n) {
      int s = X_disc(n, j);
      if (s >= 0 && s < r_j) counts[s]++;
    }
    double loglik = 0.0;
    int total = 0;
    for (int s = 0; s < r_j; ++s) total += counts[s];
    if (total == 0) return R_NegInf;

    for (int s = 0; s < r_j; ++s) {
      if (counts[s] > 0) {
        double p = (double)counts[s] / (double)total;
        loglik += counts[s] * std::log(p);
      }
    }
    int k = r_j - 1;
    return loglik - 0.5 * std::log((double)N) * k;
  }

  int q = parents.n_elem;
  std::map< std::vector<int>, std::vector<int> > table;

  for (int n = 0; n < N; ++n) {
    std::vector<int> key(q);
    for (int k = 0; k < q; ++k) key[k] = X_disc(n, parents(k));
    int child = X_disc(n, j);

    auto it = table.find(key);
    if (it == table.end()) {
      std::vector<int> vec(r_j, 0);
      if (child >= 0 && child < r_j) vec[child]++;
      table.insert(std::make_pair(key, vec));
    } else {
      if (child >= 0 && child < r_j) it->second[child]++;
    }
  }

  double loglik = 0.0;
  for (auto &kv : table) {
    std::vector<int> &vec = kv.second;
    int total = 0;
    for (int s = 0; s < r_j; ++s) total += vec[s];
    if (total == 0) continue;
    for (int s = 0; s < r_j; ++s) {
      if (vec[s] > 0) {
        double p = (double)vec[s] / (double)total;
        loglik += vec[s] * std::log(p);
      }
    }
  }

  int q_j = (int)table.size();
  int k = (r_j - 1) * q_j;
  return loglik - 0.5 * std::log((double)N) * k;
}

// Graph-level multinomial BIC
// [[Rcpp::export]]
double bic_score_multinom_cpp(const arma::imat &adj,
                              const arma::imat &X_disc,
                              const arma::ivec &n_levels) {
  if (!is_acyclic_cpp(adj)) return R_NegInf;
  int p = X_disc.n_cols;
  double total = 0.0;
  for (int j = 0; j < p; ++j) {
    arma::uvec parents = arma::find(adj.col(j) != 0);
    total += local_bic_multinom_cpp(j, parents, X_disc, n_levels);
  }
  return total;
}

// ----------------- Precompute local scores (Gaussian & Multinom) -----------------
// Return a List where each element is a NumericVector of length max_mask (2^p)
// [[Rcpp::export]]
List precompute_local_scores_gaussian(const arma::mat &X) {
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
List precompute_local_scores_multinom(const arma::imat &X_disc, const arma::ivec &n_levels) {
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

// ----------------- best_score_given_predecessors (submask enumeration) -----------------
// [[Rcpp::export]]
double best_score_given_predecessors(const NumericVector &scores_j, int preds_mask) {
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
List branch_and_bound_bn(const Nullable<arma::mat> &X_nullable = R_NilValue,
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
    if (X_disc_nullable.isNull() || n_levels_nullable.isNull()) stop("For multinomial provide X_disc and n_levels.");
    X_disc = as<arma::imat>(X_disc_nullable);
    n_levels = as<arma::ivec>(n_levels_nullable);
  }

  int p = is_gauss ? X.n_cols : X_disc.n_cols;
  if (p <= 1) {
    arma::imat adj(p, p, arma::fill::zeros);
    double s = 0.0;
    if (p == 1) {
      if (is_gauss) {
        arma::uvec empty; empty.reset();
        s = local_bic_gaussian_cpp(0, empty, X);
      } else {
        arma::uvec empty; empty.reset();
        s = local_bic_multinom_cpp(0, empty, X_disc, n_levels);
      }
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
  if (is_gauss) scores_list = precompute_local_scores_gaussian(X);
  else scores_list = precompute_local_scores_multinom(X_disc, n_levels);

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
        double best_local = best_score_given_predecessors(scores_j, preds_mask);
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
