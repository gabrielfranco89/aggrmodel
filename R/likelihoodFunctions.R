#' Build design matrix to fit aggregated model
#'
#' @name buildX
#' @param market Market data frame. MUST be a 3 column dataframe with the following order: Group, Type and Number of subjects
#' @param timeVec Vector of sampled times
#' @param basis Character indicating which basis: 'B-Splines' or 'Fourier'
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
#' @import fda
#' @import tidyr
#' @import mvtnorm
#' @export

buildX <- function(market,
                   timeVec,
                   n_basis,
                   basis = 'B-Splines',
                   n_order = 4 ## Cubic Splines (default)
                   ){

    ## Preamble
    require(fda, quietly=TRUE)
    require(tidyr, quietly=TRUE)
    C <- length(unique(market[,2]))
    J <- length(unique(market[,1]))

    ## Create basis matrix B
    t <- unique(timeVec)
    if(basis=='B-Splines')
        basisObj = create.bspline.basis(range(t),
                                        nbasis = n_basis,
                                        norder = n_order)
    if(basis=='Fourier'){
        basisObj = create.fourier.basis(range(t),
                                        nbasis = n_basis)
    }
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
#' @param covMtxList List of covariance matrices for each group
#'
#' @return Log-likehood value
#' @import mvtnorm
#' @export

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


#' Wrapper for log-likelihood evaluation
#' @name loglikWrapper
#'
#' @param pars parameters to be evaluated in optim
#' @param dataWrap dataset
#' @param mktWrap market dataset
#' @param covWrap covariance structure
#' @param corWrap correlation structure
#' @param betaWrap beta parameter
#' @param designListWrap list of design matrices for each group
#' @param nCons number of types of consumers
#' @param nBasisCov number of basis functions for functional variance
#' @export
loglikWrapper <- function(pars,
                          dataWrap,
                          mktWrap,
                          covWrap,
                          corWrap,
                          betaWrap,
                          designListWrap,
                          nCons,
                          nBasisCov ## for heterog model
                          ){
    muList <- lapply(designListWrap,
                     function(x) as.matrix(x) %*% matrix(betaWrap,
                                                         ncol=1))
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
    if(covWrap == 'Heterog'){
        ## <-------------------- UNDER CONSTRUCTION
        C <- length(unique(mktWrap[,2]))
        tvec <- unique(dataWrap$time)
        basisObj = create.bspline.basis(range(tvec),
                                        nbasis = n_basis,
                                        norder = n_order)
        B <- predict(basisObj, tvec)
        betaMC <- pars[1:(C*nBasisCov)]
        betaMtx <- cbind(beta=as.matrix(betaMC),
                         type=rep(1:C, each=n_basis))
        mcMtx <- tapply(betaMtx[,1],
                        betaMtx[,2],
                        function(x) B %*% x)
        funcVarIn <- matrix(unlist(mcMtx), ncol = C)
        sigmaList <- covMatrix(market = mktWrap,
                               group.name = 'Group',
                               type.name = 'type',
                               mkt.name = 'mkt',
                               timeVec = dataWrap$time,
                               funcMtx = funcVarIn,
                               corPar = pars[(C*nBasisCov+1):(length(pars)-C)],
                               tauPar = pars[((length(pars)-C+1):length(pars))],
                               covType = 'Heterog',
                               corType = corWrap )
        ## UNDER CONSTRUCTION -------------------->
    }
    lk <- logLikelihood(data = dataWrap,
                        muVecList = muList,
                        covMtxList = sigmaList)
    return(lk)
}


## ================================================
##   ╔═╗╔╦╗  ╔═╗┬ ┬┌┐┌┌─┐┌┬┐┬┌─┐┌┐┌┌─┐
##   ║╣ ║║║  ╠╣ │ │││││   │ ││ ││││└─┐
##   ╚═╝╩ ╩  ╚  └─┘┘└┘└─┘ ┴ ┴└─┘┘└┘└─┘
## ================================================

#' Title
#'
#' @param data
#' @param sigmaList
#' @param xbetaList
#' @param probTab
#' @param B
#'
#' @return
#' @export
#' @importFrom mvtnorm dmvnorm
#'
#' @examples
Q_function <- function(data,sigmaList,xbetaList,probTab,B){
    logLikOut=0
    for(j in unique(data$group)){
        for(i in unique(data$rep)){
            yij <- subset(data, group==j&rep==i)$y
            for(b in 1:B){
                xbeta_jb <- xbetaList[[j]][,b]
                prob_ijb <- subset(probTab, grps==j&reps==i)[,b+2]
                prob_ijb <- as.numeric(prob_ijb)
                sigma_jb <- as.matrix(sigmaList[[b]][[j]])
                logLikOut <- logLikOut +
                    prob_ijb* mvtnorm::dmvnorm(x=yij,
                                                    mean = xbeta_jb,
                                                    sigma= sigma_jb,
                                                    log = TRUE)
            } # end for b
        } # end for i
    } # end for j
    # message(paste('',round(-logLikOut,6)))
    -logLikOut ## minus for maximization
}


#' Title
#'
#' @param covPar
#' @param data
#' @param market
#' @param betaPar
#' @param piPar
#' @param B
#' @param t
#' @param K
#' @param basisFunction
#' @param n_order
#' @param pTab
#' @param I
#' @param J
#'
#' @return
#' @export
#'
#' @examples
Q_wrapper <- function(covPar,data,market,betaPar,piPar,pTab,B,t,K,I,J,
                      basisFunction,n_order){
    covPar <- matrix(covPar, nrow=B)
    ## Get probTabs
    sigMtxList <- lapply(1:B, function(b){
        sig <- covMatrix(market = market,group.name = 'group',
                         type.name = 'type',mkt.name = 'num',
                         timeVec = t,sigPar = covPar[b,1],
                         tauPar = NULL,corPar = covPar[b,2],
                         funcMtx = NULL,covType = 'Homog_Uniform',
                         corType = 'periodic',nKnots = NULL,
                         truncateDec = 8)
        sig
    })
    XList <- buildX(market = market,timeVec = t,n_basis = K,
                    basis = basisFunction,n_order = n_order)
    # X <- lapply(XList,
    #             function(x)
    #               do.call(rbind,replicate(n=I,
    #                                       expr=x,
    #                                       simplify=FALSE)
    #               )
    # )
    X <- XList
    xbeta <- lapply(X, function(x) apply(betaPar, 2, function(bt) x %*% bt))
    # yj <- split(y,grps)
    # probTab_init <- matrix(nrow=I*J, ncol=(B))
    # probTab_names <- data.frame(reps = vector(mode=class(data$rep),length = I*J),
    #                             grps = vector(mode=class(data$group),length = I*J))
    # k=1
    # for(i in unique(data$rep)){
    #     for(j in unique(data$group)){
    #         for(b in 1:B){
    #             probTab_names$reps[k] = i
    #             probTab_names$grps[k] = j
    #             yij <- subset(data, rep==i & group==j)$y
    #             probTab_init[k,b] <- mvtnorm::dmvnorm(x = as.numeric(yij),
    #                                                   mean = as.numeric(c(xbeta[[j]][,b])),
    #                                                   sigma = as.matrix(sigMtxList[[b]][[j]]),
    #                                                   log=TRUE)
    #
    #         }
    #         k=k+1
    #     }
    # }
    # for(k in 1:(I*J)){
    #     for(b in 1:B){
    #         probTab_init[k,b] <- probTab_init[k,b] + log(piPar[1,b])
    #     }
    # }
    # denomProb <- log(rowSums(exp(probTab_init)))
    # for(k in 1:(I*J)){
    #     for(b in 1:B){
    #         probTab_init[k,b] <- probTab_init[k,b] - denomProb[k]
    #     }
    #     if(all(exp(probTab_init[k,])==0)){ ## if all probs are zero
    #         b1 <- which.min(probTab_init[k,])
    #         probTab_init[k,] <- rep(0,times=B)
    #         probTab_init[k,b1] <- 1
    #     }else{
    #         probTab_init[k,] <- exp(probTab_init[k,])
    #     }
    # }
    # pTab_init <- cbind(probTab_names, probTab_init)
    ## SEND TO Q_Function
    Q_function(data=data,
               sigmaList = sigMtxList,
               xbetaList = xbeta,
               probTab = pTab,
               B=B)

    }
