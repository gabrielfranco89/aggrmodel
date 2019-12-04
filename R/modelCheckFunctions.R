#' Plot fit's mean curve
#'
#' @name plotMC
#' @param object aggrmodel object fit
#'
#' @return ggplot plot object
plotMC <- function(object){
    require(ggplot2)
    mcMtx <- object$mc
    ## Plot
    p <- ggplot(aes(x=time, y=mc), data=mcMtx) +
        geom_line() +
        facet_wrap(.~type)
    p
}


