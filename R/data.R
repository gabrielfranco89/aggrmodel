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
