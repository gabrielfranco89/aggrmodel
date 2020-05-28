##' Aux function to get initial values in aggrmodel function
##'
##' @title get_inits
##'
##' @param X Design matrix
##' @param I Number of replicates
##' @param y Dependent variable
##' @param C Number of subject types
##' @param covType Covariance structure type
##' @param market Model market
##' @param corPar_init Initial value for correlation parameters
##' @param tauPar_init Initial value for covariance parameter
##' @param truncateDec Decimal to be truncated
##' @param n_basis Number of basis function in model fit
##' @param n_basis_cov Number of basis function in heterogeneous case
##' @param t time vector
##' @param n_order Order of basis
##' @param sigPar_init inital values for sigma parameter
##' @param betaCov_init inital values for variance functional
##' @param basisFunc Basis function type
##'
##' @return An object with initial values
##' @author Gabriel Franco
get_inits <- function(X, I, data, C,
                      covType,
                      market,
                      sigPar_init,
                      corPar_init,
                      # tauPar_init,
                      betaCov_init,
                      truncateDec,
                      n_basis, n_basis_cov,
                      t, n_order, basisFunc
){
    ## Get beta init -----------------------------
    fit_init <- lm(data$y~X-1)
    beta_init <- coef(fit_init)
    sigma_fit <- summary(fit_init)$sigma/sqrt(I)

    if(covType == 'Homog_Uniform'){
        sigPar <- ifelse(is.null(sigPar_init), sigma_fit, sigPar_init)
        if(is.null(corPar_init)) corPar <- 1
        else corPar <- corPar_init
        parIn <- c(sigPar, corPar)
        lowBoundVec <- c(1e-4, 1e-4)
        ubCor <- ifelse(is.null(truncateDec), Inf, log(10^truncateDec))
        upperBoundVec <- c(Inf, ubCor)
    }

    if(covType == 'Homog'){
        if(is.null(sigPar_init))
            sigPar <- rep(sigma_fit,C)
        else{
          if(length(sigPar_init)!=C) stop("sigPar_init must have length equal",C)
            sigPar <- sigPar_init
        }
        if(is.null(corPar_init))
          corPar <- rep(1,C)
        else{
          if(length(corPar_init)!=C) stop("corPar_init must gave length equal",C)
          corPar <- corPar_init
        }
        parIn <- c(sigPar, corPar)
        lowBoundVec <- rep(1e-4,2*C)
        ubCor <- ifelse(is.null(truncateDec), Inf, log(10^truncateDec))
        upperBoundVec <- c(rep(Inf,C), rep(ubCor,C))
    }
    if(covType == 'Heterog'){
        if(is.null(n_basis_cov)){
            warning("Using the same number of basis of model fit:",n_basis)
            n_basis_cov <- n_basis
        }
        if(is.null(betaCov_init)){
             betaCov_init <- rep(log(sigma_fit), C*n_basis_cov)
        }
        else{
            if(length(betaCov_init)!=C*n_basis_cov) stop("betaCov_init must have the same length as number of basis for covariance")
        } # end if/else
        if(is.null(corPar_init))
            corPar <- rep(1, C)
        else
            corPar <- corPar_init
        parIn <- c(betaCov_init,
                   corPar)
        lowBoundVec <- c(rep(-Inf,times=(C*n_basis_cov)), ##beta
                         rep(1e-4,C)) ## corPar
        ubCor <- ifelse(is.null(truncateDec), Inf, log(10^truncateDec))
        upperBoundVec <- c(rep(Inf,times=(C*n_basis_cov)),
                           rep(ubCor,C))
     }# end if heterog

    output <- list()
    output$X <- X
    output$beta <- beta_init
    output$lb <- lowBoundVec
    output$ub <- upperBoundVec
    output$par <- parIn
    output$n_basis_cov <- n_basis_cov
    output
}


#' Fit a simple aggregated model
#' @name simpleAggrmodel
#'
#' @param market Market data frame. MUST be a 3 column dataframe with the following order: Group, Type and Number of subjects
#' @param n_basis  Number of basis functions for basis expansion
#' @param n_order Order of basis Splines (Default: 4)
#' @param basisFunction Character indicating which basis: 'B-Splines' (default) or 'Fourier'
#' @param data Dataset containing group, replicates (if any), time and aggregated signal
#'
#' @importFrom tidyr spread
#' @return lm object
#'
simpleAggrmodel <- function(data,
                            market,
                            n_type,
                            n_basis,
                            n_order,
                            basisFunction){
  mktLong <- tidyr::spread(market,type,num)
  groupList <- subset(data, select=c(group))
  mktLong <- merge(groupList,mktLong)
  X <- buildX(market=mktLong,
                  timeVec = data$time,
                  nType = n_type,
                  n_basis = n_basis,
                  n_order = n_order,
                  basis = basisFunction)
  # I <- length(unique(data$rep))
  # X <- lapply(XList,
  #             function(x)
  #               do.call(rbind,replicate(n=I,
  #                                       expr=x,
  #                                       simplify=FALSE)
  #               )
  # )
  # X <- do.call(rbind, X)
  X <- data.frame(y=data$y, X)
  fitAggr <- lm(y~.-1, data = X)
  fitAggr
}


#' Fit multiple simple models to obtain the clustering with minor SRE
#' @name getClusterInitials
#'
#' @param data Must be an ordered data frame by rep > group > time with theses names
#' @param market A data frame with columns group, type and num
#' @param n_knots Number of knots to be fitted on simple model
#' @param n_trials Number of random clustering configurations to be fitted
#' @param bFunc Basis function to be used: "B-Splines" or "Fourier"
#' @param n_order Order of basis Splines (Default: 4)
#' @param n_cluster Number of grouping clusters
#'
#' @return A list containing: data set with "group" and "cluster" indicating the fit with minor squared root error and an lm object of the chosen clustering configuration
#'
#' @importFrom purrr map2
#'
getClusterInitials <- function(data,
                               market,
                               n_knots,
                               n_order,
                               n_cluster,
                               n_trials,
                               bFunc
                               ){
  ## Separate dataset in clusters
  C <- length(unlist(unique(market[,2])))
  J <- length(unique(data$group))
  smpObl <- rep(1:n_cluster, each=C)
  if((J-length(smpObl)) <0)
    stop('Cannot fit clusters if number of groups at each cluster are smaller than number of clusters')
  smpVec <- matrix(nrow=J, ncol = n_trials)
  itmax=0
  for(i in 1:n_trials){
    if(i==1)
      smpVec[,i] <-  base::sample(c(smpObl,
                                 sample(x=1:n_cluster,
                                        size = c(J-length(smpObl)),
                                        replace = TRUE) )
      )
    else{
      flag = FALSE
      while(!flag & itmax < c(n_trials + 1e3)){
        cl <- base::sample(c(smpObl,
                                 sample(x=1:n_cluster,
                                        size = c(J-length(smpObl)),
                                        replace = TRUE) )
        )
        cl <- as.numeric(factor(cl, levels = unique(cl)))
        if(i==2){
          flag <- identical(cl, smpVec[,1])
        } else {
          flag <- apply(smpVec[,1:c(i-1)],2, function(x) identical(x,cl))
        }
        flag <- all(flag == FALSE)
        if(flag) smpVec[,i] <- cl
        itmax=itmax+1
      }
    }
  }
  if(any(is.na(smpVec))){
    smpVec <- matrix(smpVec[!is.na(smpVec)],nrow=J)
    n_trials <- ncol(smpVec)
  }
  grps <- unique(data$group)
  trial_dd <- data.frame(
    group = rep(grps, times = n_trials),
    cluster= c(smpVec),
    trial=rep(1:n_trials, each = length(grps))
  )
  mktLong <- tidyr::spread(market,type,num)
  grp_dd <- data.frame(group = data$group)
  mktLong <- merge(grp_dd,mktLong)
  X <- buildX(market=mktLong,
              timeVec = data$time,
              nType = C,
              n_basis = n_knots,
              n_order = n_order,
              basis = bFunc)
  fitList <- by(trial_dd,INDICES = trial_dd$trial,
                FUN = function(x){
                  mktc <- merge(mktLong, x, by="group")
                  x_spliter <- mktc$cluster
                  fit_dd <- data.frame(y=data$y, X)
                  fit_dd_split <- split(fit_dd, x_spliter)
                  fit_list <- lapply(fit_dd_split,
                                     function(ff){
                                       lm(y~.-1,data=ff)
                                     })
                  SR_list <- lapply(fit_list, function(x) sum(abs(resid(x))))
                  Reduce('+',SR_list)
                })

  ## Choose smaller error and output the selected cluster config
  minErr <- as.integer(which.min(fitList))
  clusterOut <- subset(trial_dd, trial == minErr)
  x <- clusterOut
  clusterOut$trial <- NULL
  rownames(clusterOut) <- NULL
  mktc <- merge(mktLong, x, by="group")
  x_spliter <- mktc$cluster
  fit_dd <- data.frame(y=data$y, X)
  fit_dd_split <- split(fit_dd, x_spliter)
  fitListOut <- lapply(fit_dd_split,
                     function(ff){
                       lm(y~.-1,data=ff)
                     })
   return(list('clusterConfig'=clusterOut,
              'fitList' = fitListOut))
}
