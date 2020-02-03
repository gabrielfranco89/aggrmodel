#' Simulated clustered aggregated signal
#'
#' A dataset containing 10 groups divided into 2 clusters (B1=4,B2=6) with 20 replications each.
#'
#' @format A data frame with 9600 obs. of 7 variables:
#' \describe{
#'   \item{group}{Grouping variable with 8 different groups}
#'   \item{rep}{Replication variable with 4 different replicates}
#'   \item{temperature}{Temperature information}
#'   \item{cluster}{Clustering variable with 2 clusters}
#'   \item{time}{Time variable with 48 samples times}
#'   \item{y}{Simulated dependent aggregated variable}
#'   \item{y_real}{Non noisy dependent aggregated variable}
#' }
#'
#' The parameters and market used for simulation are in attributes.
#' @source Simulated data by createSimuData function
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


#' Simulated Mean Curves for simulation analysis
#'
#' A dataset containing two different sets of three mean curves for simulation. Can be used for cluster analysis or standard aggregated fit
#'
#' @format A data frame with 144 obs. of 4 variables:
#' \describe{
#' \item{Time}{Observed times in 30 minutes interval}
#' \item{Type}{Type flag: A, B or C}
#' \item{Cluster1}{Simple mean curve with featured peaks}
#' \item{Cluster2}{Mean curves used in previous trafo's research.}
#' }
#' @source Simulated data in Cluster 1 and estimated mean curves in Cluster 2 obtained as in  Lenzi et al. (2017) and Dias et al. (2019,2013).
"simulatedMeanCurves"

#' Simulated Temperature for simulation analysis
#'
#' A dataset containing temperature for 20 days observed in 30 minutes intervals
#'
#' @format A data frame with 960 obs. of 3 variables:
#' \describe{
#' \item{Time}{Observed times in 30 minutes interval}
#' \item{Day}{Day flag: from 1 to 20}
#' \item{Temp}{Temperature: day 1 to 10 is a Brazilian summer scenario and day 11 to 20 a Brazilian winter scenario}
#' }
"simuTemperature"
