## Correlation Matrices ============================================

#' Create auxiliary correlation base matrix
#'
#' @name createCorBase
#' @param timeVec Vector of time
#' @return Square correlation matrix of size \code{length(unique(timeVec))}
#' @examples
#' createCorBase(seq(0,1, length.out=12))
#' @export
createCorBase <- function(timeVec){

    time <- unique(timeVec)
    T <- length(time)
    tmp <- data.frame(i = rep(time, each=T),
                 j = rep(time, times = T)
                 )
    tmp$diff = tmp$i - tmp$j


    mtxOut  <- matrix(tmp[,3],
                      nrow = T,
                      ncol = T)


    return(mtxOut)
}

## Periodic ----------------------------------------

#' Periodic correlation matrix
#'
#' @name periodicCorMtx
#' @param timeVec Vector of time
#' @param corPar Covariance parameter
#' @param truncateDec Decimal to be truncated. Default is none.
#'
#' @details
#' Create a covariance matrix of type \emph{periodic} with elements
#'
#' \deqn{r(s,t | \theta) = exp \left\{-2\theta\, sen^2\left(\pi\frac{(t-s)}{T}\right) \right\}.}
#'
#' @return Square correlation matrix of size \code{length(unique(timeVec))} of type periodic.
#'
#' @examples
#' myTime <- seq(0,1, length.out = 12)
#' periodicCorMtx(myTime, corPar = 8, truncateDec = 4)
#' @export
periodicCorMtx <- function(timeVec,
                           corPar,
                           truncateDec = NULL){
    ## require(Matrix)

    baseMtx <- createCorBase(timeVec)

    mtx <- apply(baseMtx, 2, function(x) exp(-2* (sin(pi*x)^2)*abs(corPar^2 )))

    if(!is.null(truncateDec)) mtx <- round(mtx,truncateDec)

    ## mtx <- as(mtx, 'symmetricMatrix')

    return(mtx)

}

## Ornstein–Uhlenbeck (exponential) ------------------

#' Exponential correlation matrix
#'
#' @name expCorMtx
#' @param timeVec Vector of time
#' @param corPar Covariance parameter
#' @param truncateDec Decimal to be truncated. Default is none.
#'
#' @details
#' Create a covariance matrix of type \emph{exponential} with elements
#'
#' \deqn{r(s,t | \theta) = exp \left\{-2\theta\, \left(\frac{|t-s|}{T}\right) \right\}.}
#'
#' @return Square correlation matrix of size \code{length(unique(timeVec))} of type exponential.
#'
#' @examples
#' myTime <- seq(0,1, length.out = 12)
#' expCorMtx(myTime, corPar = 8, truncateDec = 4)
#' @export
expCorMtx <- function(timeVec, corPar,
                      truncateDec = NULL){
    ## require(Matrix)

    baseMtx <- createCorBase(timeVec)

    mtx <- apply(baseMtx, 2, function(x) exp(-abs(x)*abs(corPar)))

    if(!is.null(truncateDec)) mtx <- round(mtx,truncateDec)

    ## mtx <- as(mtx, 'symmetricMatrix')

    return(mtx)

}


## VARIANCE MATRIX ========================================

#' Functional variance matrix
#'
#' @name createVarMtx
#' @param functionalVec Vector containing values for the matrix diagonal
#' @param sigPar Scale parameter (see details)
#' @param tauPar Non-negative power parameter (see details)
#'
#' @details
#' The functional variance matrix is a matrix with diagonal of elements:
#'
#' \deqn{v(t) = \sigma \, (\eta(t))^{-\tau}}
#'
#' @return A square functional variace matrix with size = \code{length(functionalVec)}
#' @export
createVarMtx <- function(functionalVec,
                            sigPar,
                            tauPar){
    diagVec = sigPar*(functionalVec^(-tauPar))
    mtx <- diag(x=diagVec) ##, names=FALSE)
    return(mtx)
}


######################################################################
#   ____             ____  _                   _
#  / ___|_____   __ / ___|| |_ _ __ _   _  ___| |_ _   _ _ __ ___  ___
# | |   / _ \ \ / / \___ \| __| '__| | | |/ __| __| | | | '__/ _ \/ __|
# | |__| (_) \ V /   ___) | |_| |  | |_| | (__| |_| |_| | | |  __/\__ \
#  \____\___/ \_/   |____/ \__|_|   \__,_|\___|\__|\__,_|_|  \___||___/
######################################################################

#' Covariance matrix for aggregated model
#'
#' @name covMatrix
#' @param market Market information (improve)
#' @param group.name empty
#' @param type.name empty
#' @param mkt.name empty
#' @param timeVec empty
#' @param sigPar improve
#' @param tauPar empty
#' @param corPar empty
#' @param funcMtx empty
#' @param covType Any of 'Homog_uniform', 'Homog' or 'Heterog'.
#' @param corType Any of 'periodic' (default) or 'exponential'.
#' @param nKnots empty
#' @param truncateDec empty
#'
#' @details building...
#'
#' @return A list of C covariance matrices of size \code{ncol(funcMtx)}
#'
#' @examples
#' set.seed(2019)
#'
#' ## Create market
#' mkt = data.frame(group = rep(1:3, each=2),
#'                  type = rep(1:2, times = 3),
#'                  value = sample(1:20, 6))
#'
#' myTimevec = seq(0,1, length.out = 4)
#' mySigPar = matrix(c(2,3), ncol=2)
#' myTauPar = matrix(c(.2, .2), ncol=2)
#' myCorPar = matrix(c(8, 12), ncol=2)
#' myFuncMtx = matrix(runif(8), nrow = 4, ncol=2)
#'
#' ## Homogeneous example
#' homogMtx = covMatrix(market = mkt, group.name = 'group', type.name = 'type', mkt.name = 'value', timeVec = myTimevec, sigPar = mySigPar, tauPar = myTauPar, corPar = myCorPar, covType = 'Homog', corType = 'periodic')
#' @export
#'
#' @import Matrix
#' @import dplyr
#' @import purrr
#' @import tidyr
covMatrix <- function(market,
                      group.name,
                      type.name,
                      mkt.name,
                      timeVec,
                      sigPar,
                      tauPar,
                      corPar,
                      funcMtx = NULL,
                      covType,
                      corType = 'periodic',
                      nKnots = NULL,
                      truncateDec = NULL
                      ){
    require(Matrix,quietly=TRUE)
    require(dplyr,quietly=TRUE)
    require(purrr,quietly=TRUE)
    require(tidyr,quietly=TRUE)
    select <- dplyr::select

    ## Preamble
    myMkt <- market
    colnames(myMkt) = c('group', 'type', 'num')
    t = unique(timeVec)
    T = length(t)
    C = length(unique(myMkt$type))
    J = length(unique(myMkt$group))

    mktComp <- myMkt  %>%
        spread(key = type,
               value = num) %>%
    ##    select(-group) %>%
    ##    as.matrix() %>%
    ##    split(1:J)
        split(.$group) %>%
        map(function(x) x %>% select(-group) %>% as.matrix())

    ## Homog Unif ::::::::::::::::::::::::::::::
    if(covType == 'Homog_Uniform'){
        if(any(length(sigPar) != 1,
               length(corPar) != 1))
            stop('Please, check number of parameters for homogeneous uniform model!')
        vc <- createVarMtx(functionalVec = rep(1,T),
                           sigPar = sigPar,
                           tauPar = 0)
        if(corType == 'periodic')
            cc <- periodicCorMtx(timeVec = t,
                                 corPar = corPar,
                                 truncateDec = truncateDec)
        if(corType == 'exponential')
            cc <- expCorMtx(timeVec = t,
                            corPar = corPar,
                            truncateDec = Par)
        covMtx  <- vc %*% cc %*% vc
        covMtxList <- tapply(market[,3], market[,1],
                       function(m){
                           tmp <- lapply(m, function(mjc) mjc*covMtx)
                           Reduce('+', tmp)
                       })
    } # end if homog unif

    ## Homog ::::::::::::::::::::::::::::::::::::::::
    if(covType == 'Homog'){
     if(any(length(sigPar) != C,
            length(corPar) != C))
            stop('Please, check number of parameters for homogeneous model!')
        covMtxListC <- lapply(1:C,
                              function(c, sigParIn, corParIn, trc, t){

                                  vc <- createVarMtx(functionalVec = rep(1,T),
                                                     sigPar = sigParIn[c],
                                                     tauPar = rep(0,C))

                                  if(corType == 'periodic')
                                      cc <- periodicCorMtx(timeVec = t,
                                                           corPar = corParIn[c],
                                                           truncateDec = trc)

                                  if(corType == 'exponential')
                                      cc <- expCorMtx(timeVec = t,
                                                      corPar = corParIn[c],
                                                      truncateDec = trc)

                                  return( vc %*% cc %*% vc)

                              },
                              corParIn = corPar,
                              sigParIn = sigPar,
                              trc = truncateDec,
                              t=t)
        covMtxList = lapply(mktComp,
                     function(mj, cMtx, C){
                         mm = lapply(1:C,
                                    function(c) mj[c] * cMtx[[c]])
                         return(Reduce('+', mm))
                     },
                     cMtx = covMtxListC,
                     C=C
                     )
    } # end if homog

    ## Heterogeneous ::::::::::::::::::::::::::::::::::::
    if(covType == 'Heterog'){

        ## NOTE ----------------------------------------
        ##   The object funcVec must be a matrix TxC of
        ##        fitted tipologies.
        ##   Can be any functional of representation Tx1. This allows
        ##     heterogeneous variance structures as in
        ##     Dias, Garcia, Schmidt (2011)
        if(any(ncol(funcMtx) != C,
               nrow(funcMtx) != T,
               length(corPar) != C))
            stop('Please, check number of parameters for proportional model!')
        covMtxListC <-
            lapply(1:C,
                   function(c, sigParIn, corParIn,
                            tauParIn, funcVecIn,
                            trc, t){
                       vc <- createVarMtx(functionalVec = funcVecIn[,c],
                                          sigPar = sigParIn[c],
                                          tauPar = tauParIn[c])
                       if(corType == 'periodic')
                           cc <- periodicCorMtx(timeVec = t,
                                                corPar = corParIn[c],
                                                truncateDec = trc)
                       if(corType == 'exponential')
                           cc <- expCorMtx(timeVec = t,
                                           corPar = corParIn[c],
                                           truncateDec = trc)
                       return( vc %*% cc %*% vc)
                   },
                   corParIn = corPar,
                   sigParIn = sigPar,
                   tauParIn = tauPar,
                   funcVecIn = funcMtx,
                   trc = truncateDec,
                   t=t) # end lapply
        covMtxList = lapply(mktComp,
                     function(mj, cMtx, C){
                         mm = lapply(1:C,
                                    function(c) mj[c] * cMtx[[c]])
                         return(Reduce('+', mm))
                     },
                     cMtx = covMtxListC,
                     C=C
                     )
    } # end if hetero
    return(covMtxList)
}








