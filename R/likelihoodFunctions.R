#' Build design matrix to fit aggregated model
#'
#' @name buildX
#' @param market Market data frame. MUST be a 3 column dataframe with the following order: Group, Type and Number of subjects
#' @param timeVec Vector of sampled times
#' @param n_basis Number of basis functions for basis expansion
#' @param n_order Order of basis splines (Default: 4)
#'
#' @return List of basis matrices for each group in market[,1]
#' @examples
#' ## Create market example and time vector
#' myMkt = data.frame(group = rep(1:4, each=3),
#'                    type = rep(1:3, times=4),
#'                    mkt = sample(1:20, size=12, replace=TRUE))
#' myTimeVec = seq(0,1, length.out=12)
#'
#' designList = buildX(myMkt, myTimeVec, n_basis = 5)
#' str(designList)

buildX <- function(market,
                   timeVec,
                   n_basis,
                   n_order = 4 ## Cubic Splines (default)
                   ){

    ## Preamble
    require(fda, quietly=TRUE)
    require(tidyr, quietly=TRUE)
    C <- length(unique(market[,2]))
    J <- length(unique(market[,1]))

    ## Create basis matrix B
    t <- unique(timeVec)
    basisObj = create.bspline.basis(range(t),
                                    nbasis = n_basis,
                                    norder = n_order)
    B = predict(basisObj, t)
    
    ## Kronecker product with market    
    X <- tapply(market[,3],
                market[,1],
                FUN=function(m) t(m) %x% B,
                simplify=FALSE)
    return(X)        
}

#' Compute log-likelihood for aggregated model
#'
#' @name logLikelihood
#'
#'
#' @param data Data Frame with 4 columns in the following order: Group, Replicates, Time, Signal
#' @param muVecList List of \eqn{X\beta} for each group
#' @param covMatrixList List of covariance matrices for each group
#'
#' @return Log-likehood value

logLikelihood <- function(data,
                          muVecList,
                          covMtxList){

    ## Preamble
    require(mvtnorm, quietly=TRUE)

    grps <- unique(data[,1])

    ## Sort data to match time, group and replicates
    ##    data <- data[order(data[,1], data[,2], data[,3]),]
    
    sumLogLik <- 0
    for(j in grps){

        subData <- subset(data, data[,1]==j)
        
        actualMu <- muVecList[[j]]
        actualSigma <- as.matrix(covMtxList[[j]])


        logLik <- tapply(subData[,4], subData[,2],
                         FUN=dmvnorm,
                         mean = actualMu,
                         sigma = actualSigma,
                         log=TRUE)

        sumLogLik <- sumLogLik + sum(logLik)
        
    }

    return(-sumLogLik)
    

}



loglikWrapper <- function(pars,
                          dataWrap,
                          mktWrap,
                          covWrap,
                          corWrap,
                          betaWrap,
                          designListWrap,
                          nCons){
    muList <- lapply(designListWrap,
                     function(x) x %*% betaWrap)
    if(covWrap == 'Homog_Uniform'){
        sigmaList <- covMatrix(market = mktWrap,
                               group.name = 'Group',
                               type.name = 'type',
                               mkt.name = 'mkt',
                               timeVec = dataWrap$time,
                               sigPar = pars[1],
                               corPar = pars[2],
                               covType = 'Homog_Uniform',
                               corType = corWrap )
    }
    if(covWrap == 'Homog'){
        C <- length(unique(mktWrap[,2]))
        sigmaList <- covMatrix(market = mktWrap,
                                  group.name = 'Group',
                                  type.name = 'type',
                                  mkt.name = 'mkt',
                                  timeVec = dataWrap$time,
                                  sigPar = pars[1:C],
                                  corPar = pars[(C+1):(2*C)],
                                  covType = 'Homog',
                                  corType = corWrap )
    }
    lk <- logLikelihood(data = dataWrap,
                        muVecList = muList,
                        covMtxList = sigmaList)
    return(lk)
}
