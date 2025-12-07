// src/heuristic_cpp.cpp
//
// C++ implementations (Gaussian + Multinomial) for heuristic BN structure
// learning: hill-climbing & tabu search.
//
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <map>
#include <vector>
#include <algorithm>
#include <string>

using namespace Rcpp;
using namespace arma;

// ---------- Acyclicity ----------

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

// ---------- Gaussian BIC ----------

double local_bic_gaussian_cpp(int j, const arma::uvec &parents, const arma::mat &X) {
  int N = X.n_rows;
  arma::vec y = X.col(j);

  if (parents.n_elem == 0) {
    double mu = arma::mean(y);
    arma::vec r = y - mu;
    double rss = arma::dot(r, r);
    double sigma2 = rss / N;
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
  arma::vec beta = arma::solve(XtX, Xty);

  arma::vec r = y - Xd * beta;
  double rss = arma::dot(r, r);
  double sigma2 = rss / N;
  double loglik = -0.5 * N * (std::log(2.0 * datum::pi * sigma2) + 1.0);

  int k = q + 2;
  return loglik - 0.5 * std::log((double)N) * k;
}

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

// ---------- Multinomial BIC ----------

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
      vec[child]++;
      table.insert(std::make_pair(key, vec));
    } else {
      it->second[child]++;
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

  int q_j = table.size();
  int k = (r_j - 1) * q_j;
  return loglik - 0.5 * std::log((double)N) * k;
}

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

// ---------- Common helpers for HC / Tabu ----------

enum ScoreType { GAUSSIAN = 0, MULTINOM = 1 };

double score_dag(const arma::imat &adj,
                 const arma::mat &X,
                 const arma::imat &X_disc,
                 const arma::ivec &n_levels,
                 ScoreType type) {
  if (type == GAUSSIAN) {
    return bic_score_gaussian_cpp(adj, X);
  } else {
    return bic_score_multinom_cpp(adj, X_disc, n_levels);
  }
}

// [[Rcpp::export]]
Rcpp::List hill_climbing_bn_cpp_internal(const arma::mat &X,
                                         const arma::imat &X_disc,
                                         const arma::ivec &n_levels,
                                         int max_iter,
                                         bool verbose,
                                         int score_type) {

  int p = X.n_cols;
  arma::imat adj(p, p, fill::zeros);
  ScoreType type = (score_type == 0) ? GAUSSIAN : MULTINOM;

  double best_score = score_dag(adj, X, X_disc, n_levels, type);
  Rcpp::NumericVector history;
  history.push_back(best_score);

  if (verbose) {
    Rcpp::Rcout << "HC C++ (" << (type == GAUSSIAN ? "gaussian" : "multinom")
                << "): initial score = " << best_score << std::endl;
  }

  for (int iter = 0; iter < max_iter; ++iter) {
    double best_move_score = best_score;
    arma::imat best_adj = adj;
    int best_from = -1, best_to = -1, best_type = -1;

    for (int i = 0; i < p; ++i) {
      for (int j = 0; j < p; ++j) {
        if (i == j) continue;

        // ADD
        if (adj(i, j) == 0) {
          arma::imat new_adj = adj;
          new_adj(i, j) = 1;
          if (is_acyclic_cpp(new_adj)) {
            double sc = score_dag(new_adj, X, X_disc, n_levels, type);
            if (sc > best_move_score) {
              best_move_score = sc;
              best_adj = new_adj;
              best_from = i; best_to = j; best_type = 0;
            }
          }
        } else {
          // REMOVE
          {
            arma::imat new_adj = adj;
            new_adj(i, j) = 0;
            double sc = score_dag(new_adj, X, X_disc, n_levels, type);
            if (sc > best_move_score) {
              best_move_score = sc;
              best_adj = new_adj;
              best_from = i; best_to = j; best_type = 1;
            }
          }
          // REVERSE
          if (adj(j, i) == 0) {
            arma::imat new_adj = adj;
            new_adj(i, j) = 0;
            new_adj(j, i) = 1;
            if (is_acyclic_cpp(new_adj)) {
              double sc = score_dag(new_adj, X, X_disc, n_levels, type);
              if (sc > best_move_score) {
                best_move_score = sc;
                best_adj = new_adj;
                best_from = i; best_to = j; best_type = 2;
              }
            }
          }
        }
      }
    }

    if (best_type == -1 || best_move_score <= best_score) {
      if (verbose) {
        Rcpp::Rcout << "HC C++: no improvement at iteration "
                    << (iter + 1) << std::endl;
      }
      break;
    }

    adj = best_adj;
    best_score = best_move_score;
    history.push_back(best_score);

    if (verbose) {
      const char *typestr = (best_type == 0 ? "add"
                               : (best_type == 1 ? "remove" : "reverse"));
      Rcpp::Rcout << "HC C++: iter " << (iter + 1)
                  << " - move " << typestr << " " << (best_from + 1)
                  << " -> " << (best_to + 1)
                  << " - score = " << best_score << std::endl;
    }
  }

  return Rcpp::List::create(
    _["adj"] = adj,
    _["score"] = best_score,
    _["history"] = history,
    _["iterations"] = (int)history.size() - 1
  );
}

// [[Rcpp::export]]
Rcpp::List tabu_search_bn_cpp_internal(const arma::mat &X,
                                       const arma::imat &X_disc,
                                       const arma::ivec &n_levels,
                                       int max_iter,
                                       int tabu_tenure,
                                       bool verbose,
                                       int score_type) {

  int p = X.n_cols;
  arma::imat adj(p, p, fill::zeros);
  ScoreType type = (score_type == 0) ? GAUSSIAN : MULTINOM;

  double current_score = score_dag(adj, X, X_disc, n_levels, type);
  double best_score    = current_score;
  arma::imat best_adj  = adj;

  Rcpp::NumericVector history;
  history.push_back(current_score);

  if (verbose) {
    Rcpp::Rcout << "Tabu C++ (" << (type == GAUSSIAN ? "gaussian" : "multinom")
                << "): initial score = " << current_score << std::endl;
  }

  std::vector< std::string > tabu_list;

  for (int iter = 0; iter < max_iter; ++iter) {
    bool found = false;
    double best_candidate_score = R_NegInf;
    arma::imat best_candidate_adj = adj;
    std::string best_edge_key;
    int best_from = -1, best_to = -1, best_type = -1;

    for (int i = 0; i < p; ++i) {
      for (int j = 0; j < p; ++j) {
        if (i == j) continue;

        std::string edge_key = std::to_string(i) + "_" + std::to_string(j);

        // ADD
        if (adj(i, j) == 0) {
          arma::imat new_adj = adj;
          new_adj(i, j) = 1;
          if (is_acyclic_cpp(new_adj)) {
            double sc = score_dag(new_adj, X, X_disc, n_levels, type);
            bool is_tabu = (std::find(tabu_list.begin(), tabu_list.end(),
                                      edge_key) != tabu_list.end());
            bool aspir = (sc > best_score);

            if (!is_tabu || aspir) {
              if (sc > best_candidate_score) {
                best_candidate_score = sc;
                best_candidate_adj = new_adj;
                best_edge_key = edge_key;
                best_from = i; best_to = j; best_type = 0;
                found = true;
              }
            }
          }

        } else {
          // REMOVE
          {
            arma::imat new_adj = adj;
            new_adj(i, j) = 0;
            double sc = score_dag(new_adj, X, X_disc, n_levels, type);
            bool is_tabu = (std::find(tabu_list.begin(), tabu_list.end(),
                                      edge_key) != tabu_list.end());
            bool aspir = (sc > best_score);

            if (!is_tabu || aspir) {
              if (sc > best_candidate_score) {
                best_candidate_score = sc;
                best_candidate_adj = new_adj;
                best_edge_key = edge_key;
                best_from = i; best_to = j; best_type = 1;
                found = true;
              }
            }
          }

          // REVERSE
          if (adj(j, i) == 0) {
            arma::imat new_adj = adj;
            new_adj(i, j) = 0;
            new_adj(j, i) = 1;
            if (is_acyclic_cpp(new_adj)) {
              double sc = score_dag(new_adj, X, X_disc, n_levels, type);
              std::string edge_key_rev =
                std::to_string(j) + "_" + std::to_string(i);
              bool is_tabu = (std::find(tabu_list.begin(), tabu_list.end(),
                                        edge_key_rev) != tabu_list.end());
              bool aspir = (sc > best_score);

              if (!is_tabu || aspir) {
                if (sc > best_candidate_score) {
                  best_candidate_score = sc;
                  best_candidate_adj = new_adj;
                  best_edge_key = edge_key_rev;
                  best_from = i; best_to = j; best_type = 2;
                  found = true;
                }
              }
            }
          }
        }
      }
    }

    if (!found) {
      if (verbose) {
        Rcpp::Rcout << "Tabu C++: no admissible move at iteration "
                    << (iter + 1) << std::endl;
      }
      break;
    }

    adj = best_candidate_adj;
    current_score = best_candidate_score;
    history.push_back(current_score);

    if (current_score > best_score) {
      best_score = current_score;
      best_adj   = adj;
    }

    tabu_list.push_back(best_edge_key);
    if ((int)tabu_list.size() > tabu_tenure) {
      tabu_list.erase(tabu_list.begin());
    }

    if (verbose) {
      const char *typestr = (best_type == 0 ? "add"
                               : (best_type == 1 ? "remove" : "reverse"));
      Rcpp::Rcout << "Tabu C++: iter " << (iter + 1)
                  << " - move " << typestr << " " << (best_from + 1)
                  << " -> " << (best_to + 1)
                  << " - score = " << current_score
                  << " - best_score = " << best_score << std::endl;
    }
  }

  return Rcpp::List::create(
    _["adj"] = adj,
    _["score"] = current_score,
    _["best_adj"] = best_adj,
    _["best_score"] = best_score,
    _["history"] = history,
    _["iterations"] = (int)history.size() - 1
  );
}
