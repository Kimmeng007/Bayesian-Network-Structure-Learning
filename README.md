# Structure Learning for Bayesian Networks

This repository contains an implementation and experimental study of **structure learning algorithms for Bayesian Networks**, focusing on learning Directed Acyclic Graphs (DAGs) from data using score-based and constraint-based approaches. The project was developed as an academic submission and explores both theoretical foundations and practical performance of Bayesian network structure learning methods.

## Project Overview

Bayesian Networks (BNs) are probabilistic graphical models that represent conditional dependencies between random variables using a Directed Acyclic Graph (DAG). Learning the structure of a BN from data is a challenging combinatorial optimization problem.

This project investigates:
- How Bayesian Network structures can be learned from data
- The trade-offs between different structure learning strategies
- Empirical evaluation of learned structures on benchmark datasets

## Key Concepts

- Bayesian Networks (BNs)
- Directed Acyclic Graphs (DAGs)
- Structure learning
- Score-based methods
- Constraint-based methods
- Conditional independence tests

## Methods Implemented

The project focuses on **structure learning**, including:

### 1. Score-Based Learning
- Searches for the best DAG by optimizing a scoring function
- Typical scores include likelihood-based or information-theoretic criteria
- Uses heuristic search strategies due to the exponential search space

### 2. Constraint-Based Learning
- Relies on statistical conditional independence tests
- Builds the graph structure by identifying dependencies and independencies
- Often faster but sensitive to statistical errors

## Experiments & Evaluation

The implementation is evaluated through:
- Synthetic or real-world datasets
- Comparison of learned graph structures
- Analysis of accuracy, complexity, and robustness

Evaluation criteria may include:
- Structural Hamming Distance (SHD)
- Computational efficiency
- Stability across runs

## How to use the package?
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

## Contributors
HONG Kimmeng, KOH Tito, NOUV Ratanakmuny
