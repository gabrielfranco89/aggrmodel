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

#' Plot fitted values over observed values
#' @name plotFitted
#' @param object aggrmodel object fit
#' @param obsAlpha Numeric between 0 and 1 indicating observed curves transparency
#' @param predColor Character or numeric indicating fitted curves color
plotFitted <- function(object,
                       obsAlpha = .5,
                       predColor = 'magenta'){
    require(ggplot2)
    dd <- object$fitted
    p <- dd %>% ggplot(aes(x=time,y=y, group=rep)) +
        geom_line(alpha=obsAlpha) +
        geom_line(aes(y=pred), col=predColor) +
        facet_wrap(.~group)
    p
}
