#' Plot fit's mean curve
#'
#' @name plot.aggrmodel
#'
#' @param object aggrmodel object fit
#' @param scales Character value indicating if scales should be 'fixed' or 'free' in the ggplot
#'
#' @return ggplot plot object
#' @import ggplot2
#' @export
plot.aggrmodel <- function(object, scales = 'fixed'){
    require(ggplot2)
    mcMtx <- object$mc
    ## Plot
    p <- ggplot(aes(x=time, y=mc), data=mcMtx) +
        geom_line() +
        facet_wrap(.~type, scales = scales)
    p
}

#' Plot fitted values over observed values
#' @name plotFitted
#' @param object aggrmodel object fit
#' @param obsAlpha Numeric between 0 and 1 indicating observed curves transparency
#' @param predColor Character or numeric indicating fitted curves color
#' @import ggplot2
#' @export
plotFitted <- function(object,
                       obsAlpha = .5,
                       predColor = 'magenta'){
    require(ggplot2)
    dd <- object#$fitted
    p <- dd %>% ggplot(aes(x=time,y=y, group=rep)) +
        geom_line(alpha=obsAlpha) +
        geom_line(aes(y=pred), col=predColor) +
        facet_wrap(.~group)
    p
}

#' Plot mean curves of Aggregated Model with clusters
#'
#' @param object aggrmodel_cluster object
#' @param scales Character value indicating if scales should be 'fixed' or 'free' in the ggplot
#'
#' @return A ggplot2 object
#' @export
#' @import ggplot2
#'
#' @examples
#' library(dplyr)
#' data = simuData %>%
#'   select(group=Group, rep=Rep, time=Time, y=Load)
#' df=data
#' mkt = market
#' mkt$Cluster=NULL
#' colnames(mkt) = c('group','type','num')
#'
#' ## Go for a walk after run the code below
#' fitCluster = aggrmodel_cluster(data = df, market=mkt, Y = 'y', timeVar =
#' 'time',groupVar = 'group', repVar = 'rep', n_basis = 7,n_basis_cov = NULL,
#' n_cluster = 2,n_trials = 1000, n_order = 4, corType = 'periodic', verbose=TRUE)
#'
#' plot(fitCluster)
plot.aggrmodel_cluster <- function(object,scales='fixed'){
    require(ggplot2)
    p <- object$mc %>%
        ggplot(aes(x=time,y=mc)) +
        geom_line() +
        facet_grid(cluster~type, scales=scales)
    p
}

#' Plot fitted values over observed values for an aggrmodel_cluster object
#' @name plotFitted_cluster
#' @param object aggrmodel object fit
#' @param obsAlpha Numeric between 0 and 1 indicating observed curves transparency
#' @param predColor Character or numeric indicating fitted curves color
#' @export
#' @import ggplot2
plotFitted_cluster <- function(object,
                       obsAlpha = .5,
                       predColor = 'magenta'){
    require(ggplot2)
    dd <- object$predData
    p <- dd %>% ggplot(aes(x=time,y=y, group=rep)) +
        geom_line(alpha=obsAlpha) +
        geom_line(aes(y=pred), col=predColor) +
        facet_wrap(.~group)
    p
}
