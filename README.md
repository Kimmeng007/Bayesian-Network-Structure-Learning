# Bayesian Network Structure Learning Package in R 

**M2BayesNet** is an R package for **Bayesian Network structure learning** using **score-based methods**. It provides several algorithms to learn the structure of a Bayesian Network and returns both the **adjacency matrix** and the **BIC score** of the learned network. The package supports **Gaussian** and **Multinomial** data only.

## Features

- Score-based Bayesian Network structure learning
- Uses **BIC (Bayesian Information Criterion)** as the scoring function
- Multiple learning strategies:
  - Naive method
  - Branch and Bound
  - Heuristic methods:
    - Hill Climbing
    - Tabu Search
- Outputs:
  - Adjacency matrix of the learned DAG
  - BIC score of the learned structure
- Supported distributions:
  - Gaussian
  - Multinomial

## Installation
The package can be used in any R environment, following the steps below.
1. Install required dependencies
```r
install.packages(c("Rcpp","RcppArmadillo","devtools","roxygen2","testthat"))
```

2. Install the package from GitHub
```r
devtools::install_github("Kimmeng007/Bayesian-Network-Structure-Learning")
```
3. Load the package
```r
library(M2BayesNet)
```
4. (For Windows user) If you get a compiler error, install Rtools:
https://cran.r-project.org/bin/windows/Rtools/

## Supported Data Types

M2BayesNet currently supports:

- **Gaussian data** (continuous variables)
- **Multinomial data** (categorical variables)

Other distributions are not supported.

## Structure Learning Methods

### 1. Naive Method

A simple structure learning approach that evaluates candidate network structures directly using the BIC score.

Returns:
- Adjacency matrix
- BIC score

### 2. Branch and Bound

An exact search method that reduces the search space using bounding techniques to efficiently find high-scoring structures.

Returns:
- Adjacency matrix
- BIC score

### 3. Heuristic Methods

Designed for larger networks where exact search is computationally expensive.

#### Hill Climbing
- Greedy local search
- Iteratively improves the network by local structure modifications

#### Tabu Search
- Local search with memory
- Uses a tabu list to avoid cycling and local optima

Returns (for all heuristic methods):
- Adjacency matrix
- BIC score

## Example Usage

```r
# Example dataset (Gaussian or Multinomial)
data <- your_data_matrix

# Naive method
result_naive <- naive_bn_structure_learning(data, distribution = "gaussian")

# Branch and Bound
result_bb <-  branch_and_bound_bn(data, distribution = "multinomial")

# Hill Climbing
result_hc <- hill_climbing_bn(X_boston, distribution = "gaussian", max_iter = 100, verbose = FALSE)

# Tabu Search
result_tabu <- tabu_search_bn(X_boston, distribution = "gaussian", max_iter = 50, tabu_tenure = 10, verbose = FALSE)

# Extract results
adjacency_matrix <- result_bb$best_adj
bic_score <- result_bb$best_score
```

## Output Format

Each structure learning function returns a list containing:

- `adjacency_matrix`  
  A square matrix representing the learned Bayesian Network structure

- `bic_score`  
  The BIC score of the learned network

---

## Limitations

- Supports only Gaussian and Multinomial distributions
- Structure learning only (no parameter learning)
- No built-in network visualization

## Contributors
HONG Kimmeng, KOH Tito, NOUV Ratanakmuny
