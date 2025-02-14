
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CausalBootstrap

## Overview

**CausalBootstrap** is an R package designed for causal inference in
stratified randomized experiments. It implements a sharp variance
estimator and a bootstrap method to estimate treatment effects from data
organized by grouping variables ($M$) and treatment indicators ($Z$).

## Installation

You can install the development version of **CausalBootstrap** from
GitHub using the `devtools` package:

``` r
# Install devtools if you haven't already
install.packages("devtools")

# Install CausalBootstrap
devtools::install_github("yu-hao-yang/CausalBootstrap")
```

## Example

Here’s a simple example of how to use the `CausalBootstrap` function:

``` r
library(CausalBootstrap)

# pair = FALSE
CausalBootstrap(Y = test_stratified$Y, Z = test_stratified$Z, M = test_stratified$M)

# pair = TRUE
CausalBootstrap(Y = test_pair$Y, Z = test_pair$Z, M = test_pair$M, pair = TRUE)
```

## Reference

Yu, H., Zhu, K., & Liu, H. (2025+). Sharp variance estimator and causal
bootstrap in stratified randomized experiments. arXiv preprint
[arXiv:2401.16667](https://arxiv.org/abs/2401.16667).
