# aggrmodel

A R Package for aggregated data

## Installation

You should install it directly from github

``` r
install.packages('devtools') ## if not installed
library(devtools)
install_github("gabrielfranco89/aggrmodel")
```

## Examples

``` r
library(aggrmodel)

## Load datasets
df = data(simuData, package='aggrmodel')
mkt = data(market, package='aggrmodel')

## Simple model
aggrFit = aggrmodel(data = df, market = mkt, Y = Load, timeVar = Time, groupVar = Group, repVar = Rep, n_basis = 7)
```
