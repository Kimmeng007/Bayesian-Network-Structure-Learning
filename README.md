# Structure Learning for Bayesian Networks

This repository contains an implementation and experimental study of **structure learning algorithms for Bayesian Networks**, focusing on learning Directed Acyclic Graphs (DAGs) from data using score-based and constraint-based approaches. The project was developed as an academic submission and explores both theoretical foundations and practical performance of Bayesian network structure learning methods.

## Project Overview

Bayesian Networks (BNs) are probabilistic graphical models that represent conditional dependencies between random variables using a Directed Acyclic Graph (DAG). Formally, given a set of random variables  
\[
\mathbf{X} = \{X_1, X_2, \dots, X_n\},
\]
a Bayesian Network defines the joint distribution as:
\[
P(\mathbf{X}) = \prod_{i=1}^{n} P\left(X_i \mid \mathrm{Pa}_G(X_i)\right),
\]
where \( \mathrm{Pa}_G(X_i) \) denotes the parents of \( X_i \) in the DAG \( G \).

Learning the structure of a BN from data is a challenging combinatorial optimization problem.

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

Score-based methods search for the optimal DAG \( G \) by maximizing a scoring function over the space of all possible DAGs:
\[
G^* = \arg\max_{G \in \mathcal{G}} \; \text{Score}(G \mid \mathcal{D}),
\]
where \( \mathcal{D} \) is the observed dataset.

A commonly used score is the Bayesian Information Criterion (BIC):
\[
\text{BIC}(G \mid \mathcal{D}) = \log P(\mathcal{D} \mid G, \hat{\theta})
- \frac{k}{2} \log N,
\]
with \( \hat{\theta} \) the maximum likelihood parameters,  
\( k \) the number of free parameters, and  
\( N \) the sample size.

Due to the super-exponential size of the DAG space, heuristic search strategies are typically employed.

### 2. Constraint-Based Learning

Constraint-based methods infer the graph structure by testing conditional independence relationships of the form:
\[
X_i \perp\!\!\!\perp X_j \mid \mathbf{Z},
\]
where \( \mathbf{Z} \subset \mathbf{X} \setminus \{X_i, X_j\} \).

Edges are removed or oriented based on these tests under the **faithfulness assumption**, yielding a graph consistent with the observed independencies.

These methods are often faster but can be sensitive to statistical errors in independence testing.

## Experiments & Evaluation

The implementation is evaluated through:
- Synthetic or real-world datasets
- Comparison of learned graph structures
- Analysis of accuracy, complexity, and robustness

Structural accuracy is commonly measured using the **Structural Hamming Distance (SHD)**:
\[
\text{SHD}(G, \hat{G}) = \#(\text{edge additions, deletions, or reversals})
\]
required to transform the learned graph \( \hat{G} \) into the true graph \( G \).

Additional evaluation criteria include:
- Computational efficiency
- Stability across runs

## How to use the package?
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
