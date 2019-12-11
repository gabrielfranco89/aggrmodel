#' Plot fit's mean curve
#'
#' @name plotMC
#' @param object aggrmodel object fit
#' @param scales Character value indicating if scales should be 'fixed' or 'free' in the ggplot
#'
#' @return ggplot plot object
plotMC <- function(object, scales = 'fixed'){
    require(ggplot2)
    mcMtx <- object$mc
    ## Plot
    p <- ggplot(aes(x=time, y=mc), data=mcMtx) +
        geom_line() +
        facet_wrap(.~type, scales = scales)
    p
}


