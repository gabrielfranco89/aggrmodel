#' Plot fit's mean curve
#'
#' @name plot.aggrmodel
#'
#' @param object aggrmodel object fit
#' @param scales Character value indicating if scales should be 'fixed' or 'free' in the ggplot
#' @param CI Logical. Plot mean curves confidence interval (default:TRUE)
#' @param temps Points of second functional component to be evaluated. If empty, it will use the 3 first quartiles
#' @param sub_temp When with double functional component, plot only the 3 first quartiles (default:TRUE)
#'
#' @return ggplot plot object
#' @import ggplot2
#' @import fda
#' @export
plot.aggrmodel <- function(object, scales = 'fixed', CI = TRUE, sub_temp=TRUE, temps=NULL){
    mc <- object$mc
    if(is.null(mc$time2)){
      p <- ggplot(aes(x=time, y=mc), data=mc) +
        geom_line()
      if(CI) p <- p +
          geom_line(aes(x=time, y=mc_lwr), linetype=2, alpha = .4)+
          geom_line(aes(x=time, y=mc_upr), linetype=2, alpha = .4)
      p + facet_wrap(.~type, scales = scales)
    }else{
      if(sub_temp){
        C = length(unique(mc$type))
        tLen = length(unique(mc$time))
        if(is.null(temps))
          temps = as.numeric(summary(mc$time2)[2:4])
        basisObj = fda::create.bspline.basis(range(mc$time),
                                             nbasis = object$n_basis,
                                             norder = object$n_order)
        B = predict(basisObj, rep(unique(mc$time),times=length(temps)))
        t2 <- rep(temps, each=tLen)
        basisObj = fda::create.bspline.basis(range(mc$time2),
                                        nbasis = object$n_basis2,
                                        norder = object$n_order2)
        B2 = predict(basisObj, t2)
        nc <- ncol(B)
        BB2 <- cbind(B,B2)
        B <- apply(BB2,1,function(x) x[c(1:nc)] %x% x[-c(1:nc)])
        B <- t(B)
        rm(BB2)
        beta_surf = object$beta[1:c(object$n_basis*object$n_basis2*C)]
        beta_upr <- beta_surf + object$betaSE[1:c(object$n_basis*object$n_basis2*C)]
        beta_lwr <- beta_surf - object$betaSE[1:c(object$n_basis*object$n_basis2*C)]
        my_mc = B %*% matrix(beta_surf, ncol=C)
        mc_upr <- B %*% matrix(beta_upr,ncol=C)
        mc_lwr <- B %*% matrix(beta_lwr,ncol=C)
        ### Mean curves data frame
        mc = data.frame(mc =c(my_mc),
                        mc_upr = c(mc_upr),
                        mc_lwr = c(mc_lwr),
                        type = rep(unique(mc$type), each=tLen*length(temps)),
                        time = rep(unique(mc$time), times = length(temps)*C),
                        time2 = rep(rep(temps,each=tLen),times=C))
      }
      p <- ggplot(aes(x=time, y=mc,group=time2), data=mc) +
        geom_line()
      if(CI)
        p <- p +
        geom_line(aes(x=time, y=mc_lwr), linetype=2, alpha = .4)+
        geom_line(aes(x=time, y=mc_upr), linetype=2, alpha = .4)
      p + facet_wrap(type~time2, scales = scales, labeller = label_both)
    }
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
                       predColor = 'magenta',
                       scales = "fixed"){
    dd <- object$fitted
    if(is.null(object$mc$time2)){
    p <- dd %>% ggplot(aes(x=time,y=y, group=rep)) +
        geom_line(alpha=obsAlpha) +
        geom_line(aes(y=pred), col=predColor) +
        facet_wrap(.~group, scales=scales)
    }else{
    p <- dd %>% ggplot(aes(x=time,y=y, group=rep)) +
        geom_line(alpha=obsAlpha) +
        geom_line(aes(y=pred,col=time2)) +
        facet_wrap(.~group,scales=scales)
    }
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
plot.aggrmodel_cluster <- function(object,scales='fixed',IC=TRUE){
  require(ggplot2)
  p <- object$mc %>%
    ggplot(aes(x=time,y=mc)) +
    geom_line()
  if(IC) p <- p +
    geom_line(aes(x=time, y=mc_lwr), linetype=2, alpha = .4)+
    geom_line(aes(x=time, y=mc_upr), linetype=2, alpha = .4)
  p + facet_grid(cluster~type, scales=scales)

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
