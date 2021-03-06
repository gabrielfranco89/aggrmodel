% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/likelihoodFunctions.R
\name{loglikWrapper}
\alias{loglikWrapper}
\title{Wrapper for log-likelihood evaluation}
\usage{
loglikWrapper(
  pars,
  dataWrap,
  mktWrap,
  sCovWrap,
  covWrap,
  corWrap,
  betaWrap,
  designWrap,
  optWrap = TRUE,
  nCons,
  parallel = FALSE,
  nCores = parallel::detectCores() - 1,
  nBasisCov,
  nOrderCov,
  verbWrap,
  positive = FALSE,
  cicle = FALSE,
  truncateDec = NULL
)
}
\arguments{
\item{pars}{parameters to be evaluated in optim}

\item{dataWrap}{dataset}

\item{mktWrap}{market dataset}

\item{sCovWrap}{Sample covariance matrix}

\item{covWrap}{covariance structure}

\item{corWrap}{correlation structure}

\item{betaWrap}{beta parameter}

\item{designWrap}{Design matrix}

\item{optWrap}{Optmization criterion should be via sample cov. matrix (TRUE) or via likelihood (more sensitive)}

\item{nCons}{number of types of consumers}

\item{parallel}{TRUE/FALSE if computation should be parallel}

\item{nCores}{Number of clusters. Default: parallel::detectCores()}

\item{nBasisCov}{number of basis functions for functional variance}

\item{nOrderCov}{order of basis functions for functional variance}

\item{verbWrap}{TRUE/FALSE prints likelihood values during optimization}

\item{positive}{Use positive restriction?}

\item{truncateDec}{Integer: Decimal to be truncated in exponential functional}
}
\description{
Wrapper for log-likelihood evaluation
}
