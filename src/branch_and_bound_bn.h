#ifndef BN_BNB_H
#define BN_BNB_H

#include <RcppArmadillo.h>
#include <vector>
#include <map>
#include <string>
#include <functional>

using namespace Rcpp;
using namespace arma;

// ----------------- Cycle detection -----------------
bool dfs_cycle(int u, const arma::imat &adj, std::vector<int> &state);
bool is_acyclic_cpp(const arma::imat &adj);

// ----------------- Gaussian BIC -----------------
double local_bic_gaussian_cpp(int j, const arma::uvec &parents, const arma::mat &X);
double bic_score_gaussian_cpp(const arma::imat &adj, const arma::mat &X);

// ----------------- Multinomial BIC -----------------
double local_bic_multinom_cpp(int j,
                              const arma::uvec &parents,
                              const arma::imat &X_disc,
                              const arma::ivec &n_levels);
double bic_score_multinom_cpp(const arma::imat &adj,
                              const arma::imat &X_disc,
                              const arma::ivec &n_levels);

// ----------------- Precompute local scores -----------------
List precompute_local_scores_gaussian(const arma::mat &X);
List precompute_local_scores_multinom(const arma::imat &X_disc,
                                      const arma::ivec &n_levels);

// ----------------- Helper functions -----------------
double best_score_given_predecessors(const NumericVector &scores_j, int preds_mask);
arma::imat build_adj_from_parent_masks(const IntegerVector &parent_masks);

// ----------------- Branch and Bound DAG search -----------------
// For Gaussian: provide X
// For Multinomial: provide X_disc and n_levels
List branch_and_bound_bn(const Nullable<arma::mat> &X_nullable = R_NilValue,
                         const Nullable<arma::imat> &X_disc_nullable = R_NilValue,
                         const Nullable<arma::ivec> &n_levels_nullable = R_NilValue,
                         std::string distribution = "gaussian");

// ----------------- Score type enum -----------------
enum ScoreType { GAUSSIAN = 0, MULTINOM = 1 };

// General DAG scoring wrapper
double score_dag(const arma::imat &adj,
                 const arma::mat &X,
                 const arma::imat &X_disc,
                 const arma::ivec &n_levels,
                 ScoreType type);

#endif // BN_BNB_H
