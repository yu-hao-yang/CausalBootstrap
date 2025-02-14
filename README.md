
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

The first dataset `test_stratified` comes from a stratified randomized
experiment (McClure et al., 2014) and the second dataset `test_pair`
comes from a paired randomized experiment (King et al., 2000).

## Reference

Yu, H., Zhu, K., & Liu, H. (2025+). Sharp variance estimator and causal
bootstrap in stratified randomized experiments. arXiv preprint
[arXiv:2401.16667](https://arxiv.org/abs/2401.16667).

McClure, E. A., Sonne, S. C., Winhusen, T., Carroll, K. M., Ghitza, U.
E., McRae-Clark, A. L., … & Gray, K. M. (2014). Achieving Cannabis
Cessation—Evaluating N-acetylcysteine Treatment (ACCENT): Design and
implementation of a multi-site, randomized controlled study in the
National Institute on Drug Abuse Clinical Trials Network. Contemporary
clinical trials, 39(2), 211-223.

King, G., Gakidou, E., Ravishankar, N., Moore, R. T., Lakin, J., Vargas,
M., … & Llamas, H. H. (2007). A “politically robust” experimental design
for public policy evaluation, with application to the Mexican universal
health insurance program. Journal of Policy Analysis and Management,
26(3), 479-506.
