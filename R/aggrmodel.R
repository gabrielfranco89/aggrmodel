#' Fit Aggregated Model
#'
#' @name aggrmodel
#' @param formula building...
#' @param data Dataset containing group, replicates (if any), time and aggregated signal
#' @param market Market data frame. MUST be a 3 column dataframe with the following order: Group, Type and Number of subjects
#' @param Y Dependent variable: aggregated signal
#' @param timeVar Name of time variable
#' @param groupVar Name of grouping variable
#' @param repVar Name of replicates variable
#' @param basisFunction Character indicating which basis: 'B-Splines' or 'Fourier'
#' @param n_basis Number of basis functions for basis expansion
#' @param n_basis_cov Number of basis functions for variance functional expansion
#' @param n_order Order of basis Splines (Default: 4)
#' @param cicleRep Indicator TRUE/FALSE if replicates are cyclical
#' @param covType Covariance functional type. One of "Homog_Uniform", "Homog" or "Heterog"
#' @param corType Correlation structure type. One of "periodic" or "exponential"
#' @param diffTol Tolerance of model covergence
#'
#' @return An aggrmodel object
#' @examples
#' df = subset(simuData, Cluster==1)
#' mkt = subset(market, Cluster==1)
#' mkt = subset(mkt, select=-Cluster)
#'
#' aggrFit = aggrmodel(data = df, market = mkt, Y = Load, timeVar = Time, groupVar = Group, repVar = Rep, n_basis = 7)
#' @import Matrix
#' @export

aggrmodel <- function(formula=NULL,
                      data,
                      market,
                      Y = NULL,
                      timeVar,
                      groupVar,
                      repVar,
                      n_basis,
                      n_basis_cov = NULL,
                      basisFunction = 'B-Splines',
                      cicleRep = FALSE,
                      n_order = 4,
                      covType = 'Homog_Uniform',
                      corType = 'periodic',
                      diffTol = 1e-6){


    ## Checklist
    covFlag = covType %in% c("Homog_Uniform", "Homog", "Heterog")
    if(!covFlag)
        stop("Wrong covType option! Argument must be one of Homog_Uniform, Homog or Heterog.")
    corFlag = corType %in% c('periodic', 'exponential')
    if(!corFlag)
        stop("Wrong corType option! Argument must be periodic or exponential.")
    if(basisFunction == 'Fourier' && n_basis%%2!=1)
        stop('For Fourier basis, n_basis must be an odd number!')

    ## Preamble
    require(Matrix,quietly=TRUE)

    y = data[[substitute(Y)]]
    t = data[[substitute(timeVar)]]
    t = t/max(t)
    grps = data[[substitute(groupVar)]]
    reps = data[[substitute(repVar)]]

    J = length(unique(grps))
    I = length(unique(reps))
    T = length(unique(t))
    C = length(unique(market[,2]))

    dd <- data.frame(group=grps,
                     rep=reps,
                     time=t,
                     y = y)
    dd <- dd[order(dd[,2], dd[,1], dd[,3]),]
    y <- dd[,4]

    ## Check market
    if(ncol(market) != 3)
        stop('Market must have 3 columns in the following order: Group, Type and Number of subjects. Please check your market input.')

    ## Build disgragation basis expansion design matrix
    XList <- buildX(market=market,
                    timeVec = t,
                    n_basis = n_basis,
                    n_order = n_order,
                    basis = basisFunction)
    if(!is.null(formula)){
        cvrtMtx <- as.data.frame(model.matrix(formula, data = data))
        cvrtMtx <- cbind(cvrtMtx,
                         group=grps,
                         rep=reps,
                         time=t
                         )
        cvrtMtx <- cvrtMtx[order(cvrtMtx[['rep']],
                                 cvrtMtx[['group']],
                                 cvrtMtx[['time']]),]
        ## same for all replicates
        cvrtMtx <- subset(cvrtMtx,
                          rep == unique(reps)[1])
        cvrtMtx <- split(cvrtMtx, f=cvrtMtx$group)
        for(j in names(XList)){
            XList[[j]] <- cbind(as.data.frame(XList[[j]]),
                                cvrtMtx[[j]]
                                )
            XList[[j]] <-as.matrix(
                subset(XList[[j]],
                       select=-c(rep, group, time))
                )
            }
        rm(cvrtMtx)
    } ## end if is.null(formula)
    X <- lapply(XList,
                 function(x)
                     do.call(rbind,replicate(n=I,
                                             expr=x,
                                             simplify=FALSE)
                             )
                 )
    X <- do.call(rbind, X)
    X <- cbind(dd, X)
    X <- X[,-c(1:3)]
    fit_init <- lm(y~.-1, data = X)
    beta_init <- coef(fit_init)
    sigma_init <- sqrt(summary(fit_init)$sigma)

    ## Step 1: obtain Sigma estimates
    ## 1.1 Create parVec according to covType and corType
    if(covType == 'Homog_Uniform'){
        sigPar <- sigma_init
        corPar <- 180
        parIn <- c(sigPar, corPar)
        lowBoundVec <- c(-Inf, 0)
    }

    if(covType == 'Homog'){
        sigPar <- rep(sigma_init,C)
        corPar <- rep(180, C)
        parIn <- c(sigPar, corPar)
        lowBoundVec <- c(rep(-Inf,C), rep(0,C))
    }
    if(covType == 'Heterog'){
        ## <----------------------- UNDER CONSTRUCTION
        if(is.null(n_basis_cov)){
            warning(paste("Using the same number of basis of model fit:",n_basis))
            n_basis_cov <- n_basis
            sigPar <- beta_init
        } else {
            ## fit new object for heterog.
            XListCov <- buildX(market=market,
                            timeVec = t,
                            n_basis = n_basis_cov,
                            n_order = n_order,
                            basis = basisFunction)
            Xcov <- lapply(XListCov,
                        function(x)
                            do.call(rbind,replicate(n=I,
                                                    expr=x,
                                                    simplify=FALSE)
                                    )
                        )
            Xcov <- do.call(rbind, Xcov)
            Xcov <- cbind(dd, Xcov)
            Xcov <- Xcov[,-c(1:3)]
            fit_init_cov <- lm(y~.-1, data = Xcov)
            beta_init_cov <- coef(fit_init)
            sigPar <- beta_init_cov
        } # end if/else
        sigPar <- beta_init
        corPar <- rep(180, C)
        tauPar <- rep(.5, C)
        parIn <- c(sigPar, corPar, tauPar)
        lowBoundVec <- c(rep(-Inf,C*n_basis_cov),
                         rep(0,2*C)) ## tau's > 0?
        ## UNDER CONSTRUCTION ---------------------->
    }# end if heterog

    ## While preamble
    betaIn = beta_init
    lkDiff = diffTol + 1
    lkIn = Inf

    ## While loop
    while(lkDiff > diffTol){

        ## W.1 Obtain covariance estimates via optim
        opt <- optim(par = parIn,
                     fn = loglikWrapper,
                     covWrap = covType,
                     corWrap = corType,
                     dataWrap = dd,
                     mktWrap = market,
                     betaWrap = betaIn,
                     designListWrap = XList,
                     nCons = C,
                     lower = lowBoundVec,
                     method = "L-BFGS-B",
                     nBasisCov = n_basis_cov)
        parOut <- opt$par
        lkOut <- opt$value
        ## message('\nLikehood =', lkOut)

        ## W.2 Update Sigma estimates
        if(covType == 'Homog_Uniform'){
            sigmaOutList <- covMatrix(market = market,
                                      group.name = 'Group',
                                      type.name = 'type',
                                      mkt.name = 'mkt',
                                      timeVec = dd$time,
                                      sigPar = parOut[1],
                                      corPar = parOut[2],
                                      covType = 'Homog_Uniform',
                                      corType = corType )
        }
        if(covType == 'Homog'){
            sigmaOutList <- covMatrix(market = market,
                                      group.name = 'Group',
                                      type.name = 'type',
                                      mkt.name = 'mkt',
                                      timeVec = dd$time,
                                      sigPar = parOut[1:C],
                                      corPar = parOut[(C+1):(2*C)],
                                      covType = 'Homog',
                                      corType = corType )
        }

        sigmaInvList <- lapply(sigmaOutList, qr.solve)
        ## W.3 Obtain beta estimates with previous steps

        sigmaInv <- bdiag(sigmaInvList)
        sigmaInv <- bdiag(replicate(n=I,
                                    expr=sigmaInv,
                                    simplify=FALSE))
        X <- do.call(rbind, XList)
        X <- do.call(rbind, replicate(n=I,
                                      expr=X,
                                      simplify=FALSE))
        betaLeft <- t(X)%*%sigmaInv%*%X
        betaRight <- t(X)%*%sigmaInv%*%y
        if(cicleRep){
            ## Note: Here we use a(0) = a(T) restriction
            if(basisFunction=='B-Splines')
                basisObj = create.bspline.basis(range(t),
                                                nbasis = n_basis,
                                                norder = n_order)
            if(basisFunction=='Fourier'){
                basisObj = create.fourier.basis(range(t),
                                                nbasis = n_basis)
            }
            B = predict(basisObj, t)
            deltaVec <- matrix(B[1,] - B[nrow(B),],nrow=1)
            ## Get lambdas
            ## I can vectorize later (maybe)
            m1List <- split(x=qr.solve(betaLeft),
                            f=rep(1:C, each=n_basis))
            m2List <- split(x=betaRight,
                            f=rep(1:C, each=n_basis))
            lambda <- matrix(numeric(C),nrow=1)
            for(c in 1:C){
                lUp <- deltaVec%*%
                    matrix(m1List[[c]],
                           nrow=n_basis)[,(n_basis*(c-1)+1):(c*n_basis)]%*%
                    matrix(m2List[[c]],ncol=1)
                lDown <- deltaVec%*%
                    matrix(m1List[[c]],
                           nrow=n_basis)[,(n_basis*(c-1)+1):(c*n_basis)]%*%
                    t(deltaVec)
                lambda[1,c] <- lUp/lDown
            }
            Rvec <- lambda %x% deltaVec
            betaOut <- Matrix::solve(betaLeft, (betaRight - t(Rvec)))
        } else {
            betaOut <- Matrix::solve(betaLeft, betaRight)
        } ## end if/else cicle

        ## W.4 UPDATES
        lkDiff <- abs(lkOut - lkIn)
        betaIn <- as.numeric(betaOut)
        parIn  <- parOut
        lkIn <- lkOut
    } ## End while(lkDiff > diffTol)

    ## Output Mean curves
    if(basisFunction=='B-Splines')
        basisObj = create.bspline.basis(range(t),
                                        nbasis = n_basis,
                                        norder = n_order)
    if(basisFunction=='Fourier'){
        basisObj = create.fourier.basis(range(t),
                                        nbasis = n_basis)
    }
    tuni <- unique(t)
    B = predict(basisObj, tuni)
    ## Separate betas
    betaMC <- betaOut[1:(C*n_basis)]
    betaMtx <- cbind(beta=as.matrix(betaMC),
                     type=rep(1:C, each=n_basis))
    mcMtx <- tapply(betaMtx[,1],
                    betaMtx[,2],
                    function(x) B %*% x)
    mcMtx <- data.frame(mc=unlist(mcMtx),
                        time=rep(tuni,times=C),
                        type=rep(unique(market[,2]), each=length(tuni)))
    ## Output data with predicted values
    dd$pred <- as.numeric(X %*% betaOut)
    ## Return
    outList <- list('beta' = as.matrix(betaOut),
                'pars' = parOut,
                'mc' = mcMtx,
                'n_basis' = n_basis,
                'n_order' = n_order,
                'fitted' = dd)
    if(!is.null(formula))
        outList[['formula']] <- formula
    class(outList)='aggrmodel'
    return(outList)
}


