
<!-- README.md is generated from README.Rmd. Please edit that file -->
aggrmodel
=========

<!-- badges: start -->
<!-- badges: end -->
The goal of aggrmodel is to ...

Installation
------------

You can install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("gabrielfranco89/aggrmodel")
```

Example
-------

This is a basic example which shows you how to solve a common problem:

``` r
library(aggrmodel)

df = subset(simuData, Cluster==1)
mkt = subset(market, Cluster==1)
mkt = subset(mkt, select=-Cluster)

aggrFit = aggrmodel(data = df, market = mkt, Y = Load, timeVar = Time, groupVar = Group, repVar = Rep, n_basis = 7)
```
