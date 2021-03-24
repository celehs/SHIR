
# SHIR

<!-- badges: start -->

<!-- badges: end -->

## Overview

The algorithm performs SHIR, a novel federated learning approach used
for aggregating high dimensional and heterogeneous data from local
sites, in order to obtain efficient estimators of the logistic model.
SHIR protects individual data through a summary-statistics-based
integrating procedure. At each local site, it fits LASSO to derive
summary data that is free of the individual level information.
Susequently, at the central node, it aggregates the derived statistics
at the local sites, and produces the integrative estimator.

## Background

Meta-analyzing multiple studies allows for more precise estimates and
enables investigation of generalizability. However, in the presence of
heterogeneity across studies and high dimensional predictors, such
integrative analysis is highly challenging. An major application of such
integrative analysis is to develop generalizable predictive models using
electronic health records (EHR) data from different hospitals. EHR data
is subject to high dimensional features and important privacy
constraints. The procedure SHIR protects individual data through
summary-statistics-based integrating procedure, accommodates between
study heterogeneity in both the covariate distribution and model
parameters, and attains consistent variable selection.

## Algorithm description

![Algorithm flowchart](./figures/Flowchart_SHIR.png)

## Installation

You can install the stable version of SHIR from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("SHIR")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("celehs/SHIR")
```

## Getting started

Follow the main steps displayed in the example, in which we apply SHIR
to a simulated dataset.

## Citation

Cai, T., Liu, M., & Xia, Y. (2019). Individual Data Protected
Integrative Regression Analysis of High-dimensional Heterogeneous Data.
arXiv: Methodology.
