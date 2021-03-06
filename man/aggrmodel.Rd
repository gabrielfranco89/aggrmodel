% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggrmodel.R
\name{aggrmodel}
\alias{aggrmodel}
\title{Fit Aggregated Model}
\usage{
aggrmodel(
  formula = NULL,
  data,
  market,
  Y = NULL,
  timeVar,
  timeVar2 = NULL,
  groupVar,
  repVar,
  n_basis,
  n_basis_cov = NULL,
  basisFunction = "B-Splines",
  n_basis2 = NULL,
  n_order2 = NULL,
  basisFunction2 = NULL,
  cicleRep = FALSE,
  n_order = 4,
  covType = "Homog_Uniform",
  corType = "exponential",
  optSampleCovMatrix = FALSE,
  sigPar_init = NULL,
  corPar_init = NULL,
  tauPar_init = NULL,
  betaCov_init = NULL,
  returnFitted = TRUE,
  positive_restriction = FALSE,
  optimMethod = "L-BFGS-B",
  use_parallel = FALSE,
  n_cores = parallel::detectCores() - 1,
  truncateDec = NULL,
  verbose = FALSE,
  optVerbose = FALSE,
  useGrad = FALSE,
  diffTol = 1e-05,
  itMax = 100
)
}
\arguments{
\item{formula}{building...}

\item{data}{Dataset containing group, replicates (if any), time and aggregated signal}

\item{market}{Market data frame. MUST be a 3 column dataframe with the following order: Group, Type and Number of subjects}

\item{Y}{Dependent variable: aggregated signal}

\item{timeVar}{Name of time variable}

\item{timeVar2}{Name of second functional}

\item{groupVar}{Name of grouping variable}

\item{repVar}{Name of replicates variable}

\item{n_basis}{Number of basis functions for basis expansion}

\item{n_basis_cov}{Number of basis functions for variance functional expansion}

\item{basisFunction}{Character indicating which basis: 'B-Splines' or 'Fourier'}

\item{n_basis2}{Number of basis for second functional}

\item{n_order2}{Order for second functional expansion}

\item{basisFunction2}{Character indicating which basis: 'B-Splines' or 'Fourier'}

\item{cicleRep}{Indicator TRUE/FALSE if replicates are cyclical}

\item{n_order}{Order of basis Splines (Default: 4)}

\item{covType}{Covariance functional type. One of "Homog_Uniform", "Homog" or "Heterog"}

\item{corType}{Correlation structure type. One of "periodic" or "exponential"}

\item{optSampleCovMatrix}{Optmization criterion via sample covariance matrix convergence (TRUE and default) or via likelihood (more sensitive)}

\item{sigPar_init}{Inital values for sigma}

\item{corPar_init}{Numeric: Initial value for correlation parameters (default:20)}

\item{tauPar_init}{Numeric: Initial value for expoent parameters of complete covariance (default:0.5)}

\item{betaCov_init}{Inital values for variance functional expansion}

\item{returnFitted}{Should the fitted values be returned in output?}

\item{positive_restriction}{TRUE/FALSE if mean curves are strictly positive}

\item{optimMethod}{Choose optim method (Default: L-BFGS-B)}

\item{use_parallel}{TRUE/FALSE if computation should be parallel}

\item{n_cores}{Number of clusters. Default: parallel::detectCores()}

\item{truncateDec}{Decimal to be truncated at covariance matrix}

\item{verbose}{TRUE/FALSE prints likelihood values during optimization}

\item{optVerbose}{Print parameters while in optmization}

\item{useGrad}{Use gradient function approximation for optimization? (Default: FALSE)}

\item{diffTol}{Tolerance of model covergence}

\item{itMax}{Maximum number of iterations (Default: 100)}
}
\value{
An aggrmodel object
}
\description{
Fit Aggregated Model
}
\examples{

set.seed(81453)
df <- createSimuData(B1 = 8,
                    nRep=10)

mkt <- attr(df,"market")
df <- subset(df, group <= 8) # get only cluster 1
mkt <- subset(mkt,group<=8)

fit <-
 aggrmodel(
   data = df,
   market = mkt,
   Y = obs,
   timeVar = time,
   groupVar = group,
   repVar = rep,
   n_basis = 8,
   covType = "Homog_Uniform",
   corType = "exponential",
   returnFitted = TRUE,
   use_parallel = TRUE
 )

plot(fit)

}
