#' Build design matrix to fit aggregated model
#'
#' @name buildX
#'
#' @param timeVec Vector of sampled times
#' @param basis Character indicating which basis: 'B-Splines' or 'Fourier'
#' @param n_basis Number of basis functions for basis expansion
#' @param n_order Order of basis splines (Default: 4)
#' @param marketLong Market indexed by group and rep
#' @param nType Number of subjects type
#' @param timeVec2 Vector of second functional component
#' @param n_basis2 Number of basis function for timeVec2 expansion
#' @param basis2 Character indicating which basis: 'B-Splines' or 'Fourier'
#' @param n_order2 Order of B-Splines (Default:4)
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
#' @export

buildX <- function(marketLong,
                   nType,
                   timeVec,
                   n_basis,
                   basis = 'B-Splines',
                   n_order = 4,
                   timeVec2=NULL,
                   n_basis2=NULL,
                   basis2=NULL,
                   n_order2=NULL
                   ){
    C <- nType
    t <- timeVec
    if(basis=='B-Splines')
        basisObj = create.bspline.basis(range(t),
                                        nbasis = n_basis,
                                        norder = n_order)
    if(basis=='Fourier'){
        basisObj = create.fourier.basis(range(t),
                                        nbasis = n_basis)
    }
    B = predict(basisObj, t)
    ## Second functional component, if available
    if(!is.null(timeVec2)){
        n_basis2 <- ifelse(is.null(n_basis2),n_basis,n_basis2)
        basis2 <- ifelse(is.null(basis2),basis,basis2)
        n_order2 <- ifelse(is.null(n_order2),n_order,n_order2)
        t2 <- timeVec2
        if(basis2=='B-Splines')
            basisObj = create.bspline.basis(range(t2),
                                            nbasis = n_basis2,
                                            norder = n_order2)
        if(basis=='Fourier'){
            basisObj = create.fourier.basis(range(t2),
                                            nbasis = n_basis2)
        }
        B2 = predict(basisObj, t2)
        nc <- ncol(B)
        BB2 <- cbind(B,B2)
        B <- apply(BB2,1,function(x) x[c(1:nc)] %x% x[-c(1:nc)])
        B <- t(B)
        rm(BB2)
    }
    tmp <- list()
    for(c in 1:C){
        tmp[[c]] <- marketLong[,1+c]*B
    }
    do.call(cbind,tmp)
    # X <- cbind(marketLong[,-c(1:2)],B)
    #X <- apply(X,1,function(x) x[1:C]%x%x[-c(1:C)])
    #t(X)
}

#' Compute log-likelihood for aggregated model
#'
#' @name logLikelihood
#'
#'
#' @param data Data Frame with 4 columns in the following order: Group, Replicates, Time, Signal
#' @param muVec Vector \eqn{X\beta} of fitted values
#' @param covMtxList List of covariance matrices for each group
#'
#' @return Log-likehood value
#' @importFrom mvtnorm dmvnorm
#' @export
logLikelihood <- function(data,
                          muVec,
                          covMtxList){
    data$flag <- as.factor(paste(data$rep,data$group))
    data$flag <- as.integer(data$flag)
    data$mu <- muVec
    lk <- by(data = data,
                INDICES = data$flag,
                FUN = function(dt){
                    jj <- dt$group[1]
                    actualMu <- dt$mu
                    actualSigma <- covMtxList[[jj]]
                    mvtnorm::dmvnorm(x = dt$y,
                                     mean = actualMu,
                                     sigma = actualSigma,
                                     log = TRUE)
                })
    -sum(unlist(lk))
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
#' @param designWrap Design matrix
#' @param nCons number of types of consumers
#' @param nBasisCov number of basis functions for functional variance
#' @param nOrderCov order of basis functions for functional variance
#' @param verbWrap TRUE/FALSE prints likelihood values during optimization
#' @param truncateDec Integer: Decimal to be truncated in exponential functional
#' @export
loglikWrapper <- function(pars,
                          dataWrap,
                          mktWrap,
                          covWrap,
                          corWrap,
                          betaWrap,
                          designWrap,
                          nCons,
                          nBasisCov,
                          nOrderCov, ## for heterog model
                          verbWrap,
                          positive = FALSE,
                          truncateDec = NULL
                          ){
    if(positive){
        nCoef <- ncol(designWrap)
        betaWrap <- pars[c(1:nCoef)]
        pars <- pars[-c(1:nCoef)]
    }
    mu <- designWrap %*% betaWrap
    if(positive) mu <- exp(mu)
    if(covWrap == 'Homog_Uniform'){
        sigmaList <- covMatrix(market = mktWrap,
                               group.name = 'Group',
                               type.name = 'type',
                               mkt.name = 'mkt',
                               timeVec = dataWrap$time,
                               sigPar = pars[1],
                               corPar = ifelse(positive,exp(pars[2]),pars[2]),
                               covType = 'Homog_Uniform',
                               corType = corWrap,
                               truncateDec = truncateDec)
    }
    if(covWrap == 'Homog'){
        C <- length(unique(mktWrap[,2]))
        cp <- pars[(C+1):(2*C)]
        if(positive) cp <- exp(cp)
        sigmaList <- covMatrix(market = mktWrap,
                               group.name = 'Group',
                               type.name = 'type',
                               mkt.name = 'mkt',
                               timeVec = dataWrap$time,
                               sigPar = pars[1:C],
                               corPar = cp,
                               covType = 'Homog',
                               corType = corWrap,
                               truncateDec = truncateDec)
    }
    if(covWrap == 'Heterog'){
        ## <-------------------- UNDER CONSTRUCTION
        C <- length(unique(mktWrap[,2]))
        tvec <- unique(dataWrap$time)
        basisObj = create.bspline.basis(range(tvec),
                                        nbasis = nBasisCov,
                                        norder = nOrderCov)
        B <- predict(basisObj, tvec)
        betaMC <- pars[1:(C*nBasisCov)]
        funcVarIn <- B %*% matrix(betaMC, ncol=3)
        # betaMtx <- cbind(beta=as.matrix(betaMC),
        #                  type=rep(1:C, each=nBasisCov))
        # mcMtx <- tapply(betaMtx[,1],
        #                 betaMtx[,2],
        #                 function(x) B %*% x)
        # funcVarIn <- matrix(unlist(mcMtx), ncol = C)
        sigParIn <- pars[(C*nBasisCov+1):(length(pars)-(2*C))]
        corParIn <- pars[(C*nBasisCov+C+1):(length(pars)-C)]
        tauParIn  <- pars[((length(pars)-C+1):length(pars))]
        if(positive){
        corParIn <- exp(corParIn)
        tauParIn <- exp(tauParIn)
        }
        sigmaList <- covMatrix(market = mktWrap,
                               group.name = 'Group',
                               type.name = 'type',
                               mkt.name = 'mkt',
                               timeVec = dataWrap$time,
                               funcMtx = funcVarIn,
                               sigPar = sigParIn,
                               corPar = corParIn,
                               tauPar = tauParIn,
                               covType = 'Heterog',
                               corType = corWrap,
                               truncateDec = truncateDec)
    }
    lk <- logLikelihood(data = dataWrap,
                        muVec = mu,
                        covMtxList = sigmaList
                        )
    if(is.infinite(lk))
        lk <- .Machine$double.xmax

    if(all(verbWrap & covWrap=='Heterog'))
        message("\n lk =", round(lk,6),
                "\n mean(betaPar) =", paste(round(apply(funcVarIn,2,mean),
                                               4),
                                         collapse=','),
                "\n mean(sigPar) =", paste(round(sigParIn,4), collapse=','),
                "\n corPar =", paste(round(corParIn,4), collapse=','),
                "\n tauPar =", paste(round(tauParIn,4), collapse=',')
                )
        if(all(verbWrap & !covWrap=='Heterog'))
        message("\n lk =", round(lk,6),
                "\n pars =", paste(round(pars,2),
                                           collapse=',')
                )
    return(lk)
}


## ================================================
##   ╔═╗╔╦╗  ╔═╗┬ ┬┌┐┌┌─┐┌┬┐┬┌─┐┌┐┌┌─┐
##   ║╣ ║║║  ╠╣ │ │││││   │ ││ ││││└─┐
##   ╚═╝╩ ╩  ╚  └─┘┘└┘└─┘ ┴ ┴└─┘┘└┘└─┘
## ================================================

#' Function to be maximized in aggrmodel_cluster function
#'
#' @param data Data Frame with 4 columns in the following order: Group, Replicates, Time, Signal
#' @param sigmaList List of covariance matrices for each group
#' @param xbetaList List of estimated signal for each group
#' @param probTab Data frame of expected probabilities from E-Step containing a group variable, replicate variable and cluster probabilities
#' @param B Number of grouping clusters
#'
#' @return A scalar likelihood value
#' @export
#' @importFrom mvtnorm dmvnorm
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


#' Wrapper for M-Step function evaluation
#'
#' @param covPar Covariance parameters
#' @param data Data Frame with 4 columns in the following order: Group, Replicates, Time, Signal
#' @param market Data Frame with 4 columns in the following order: Group, Replicates, Time, Signal
#' @param betaPar Mean curves parameters
#' @param piPar Probabilities parameters
#' @param B Number of grouping clusters
#' @param t Time vector
#' @param K Number of basis functions
#' @param basisFunction Character indicating which basis: 'B-Splines' (default) or 'Fourier'
#' @param n_order Order of basis Splines (Default: 4)
#' @param pTab Data frame of expected probabilities from E-Step containing a group variable, replicate variable and cluster probabilities
#' @param I Number of replicates
#' @param J Number of groups
#' @param C Number of consumer types
#' @param covWrap covariance structure
#' @param corWrap correlation structure
#'
#' @return
#' @export
Q_wrapper <- function(covPar,data,market,betaPar,piPar,pTab,B,t,K,C,I,J,
                      basisFunction,n_order,
                      covWrap, corWrap){
    if(covWrap == 'Homog_Uniform'){
    covPar <- matrix(covPar, nrow=B)
    sigMtxList <- lapply(1:B, function(b){
        sig <- covMatrix(market = market,group.name = 'group',
                         type.name = 'type',mkt.name = 'num',
                         timeVec = t,sigPar = covPar[b,1],
                         tauPar = NULL,corPar = covPar[b,2],
                         funcMtx = NULL,covType = 'Homog_Uniform',
                         corType = corWrap, nKnots = NULL,
                         truncateDec = 8)
        sig
    })
    } ## end if cov homog_uniform
    if(covWrap == 'Homog'){
        covPar <- matrix(covPar, ncol=B,byrow=TRUE)
        sigMtxList <- lapply(1:B, function(b){
            sig <- covMatrix(market = market,group.name = 'group',
                             type.name = 'type',mkt.name = 'num',
                             timeVec = t,sigPar = covPar[1:C,b],
                             tauPar = NULL,corPar = covPar[(C+1):(2*C),b],
                             funcMtx = NULL,covType = 'Homog',
                             corType = corWrap, nKnots = NULL,
                             truncateDec = 8)
            sig
        })
    }
    XList <- buildX(market = market,timeVec = t,n_basis = K,
                    basis = basisFunction,n_order = n_order)
    X <- XList
    xbeta <- lapply(X, function(x) apply(betaPar, 2, function(bt) x %*% bt))
    ## SEND TO Q_Function
    Q_function(data=data,
               sigmaList = sigMtxList,
               xbetaList = xbeta,
               probTab = pTab,
               B=B)
}
