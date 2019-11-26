#' Simulated clustered aggregated signal
#'
#' A dataset containing 8 groups divided into 2 clusters with 4 replications each.
#'
#' @format A data frame with 1536 obs. of 8 variables:
#' \describe{
#'   \item{Group}{Grouping variable with 8 different groups}
#'   \item{C1}{Market information: customer of type C1}
#'   \item{C2}{Market information: customer of type C2}
#'   \item{C3}{Market information: customer of type C3}
#'   \item{Rep}{Replication variable with 4 different replicates}
#'   \item{Cluster}{Clustering variable with 2 clusters}
#'   \item{Time}{Time variable with 48 samples times}
#'   \item{Load}{Simulated dependent aggregated variable}
#' }
#' @source Simulated data
"simuData"

#' Market for simulated data
#'
#' A dataset containing market information for simulated aggregated data.
#'
#' @format A data frame with 24 obs. of 4 variables:
#' \describe{
#' \item{Group}{Grouping variable with 8 different groups}
#' \item{Cluster}{Clustering variable with 2 clusters}
#' \item{Type}{Subject type}
#' \item{Num}{Number of subjects of this type in this}
#' }
#' @source Simulated data
"market"
