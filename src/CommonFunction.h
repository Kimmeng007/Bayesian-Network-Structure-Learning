#ifndef COMMONFUNCTION_H
#define COMMONFUNCTION_H

#include <RcppArmadillo.h>
#include <map>
#include <vector>
#include <algorithm>
#include <string>

using namespace Rcpp;
using namespace arma;

// ---------- Acyclicity ----------
bool dfs_cycle(int u, const arma::imat &adj, std::vector<int> &state);

bool is_acyclic_cpp(const arma::imat &adj);

// ---------- Gaussian BIC ----------
double local_bic_gaussian_cpp(int j,
                              const arma::uvec &parents,
                              const arma::mat &X);

double bic_score_gaussian_cpp(const arma::imat &adj,
                              const arma::mat &X);

// ---------- Multinomial BIC ----------
double local_bic_multinom_cpp(int j,
                              const arma::uvec &parents,
                              const arma::imat &X_disc,
                              const arma::ivec &n_levels);

double bic_score_multinom_cpp(const arma::imat &adj,
                              const arma::imat &X_disc,
                              const arma::ivec &n_levels);

#endif // COMMONFUNCTION_H
