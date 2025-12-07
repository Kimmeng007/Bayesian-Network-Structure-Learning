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
#include "CommonFunction.h"

using namespace Rcpp;
using namespace arma;

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
