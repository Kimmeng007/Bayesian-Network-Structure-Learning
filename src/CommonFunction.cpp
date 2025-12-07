// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <map>
#include <vector>
#include <algorithm>
#include <string>
#include "CommonFunction.h"

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
