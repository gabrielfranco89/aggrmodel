
<!-- README.md is generated from README.Rmd. Please edit that file -->
aggrmodel <img src='man/figures/logo.png' align="right" height="139" />
=======================================================================

<!-- badges: start -->
<!-- badges: end -->
A R package to fit aggregated data model.

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

aggrFit = aggrmodel(data = df,        ## Your data
                    market = mkt,     ## Data market
                    Y = Load,         ## Dep. variable
                    timeVar = Time,   ## Time variable
                    groupVar = Group, ## Group variable
                    repVar = Rep,     ## Replicate variable
                    n_basis = 7)      ## Number of basis functions
```
