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
#' @param n_basis Number of basis functions for basis expansion
#' @param n_order Order of basis Splines (Default: 4)
#' @param cicleRep Indicator TRUE/FALSE if replicates are cyclical
#' @param covType Covariance functional type. One of "Homog_Uniform", "Homog" or "Heterog"
#' @param corType Correlation structure type. One of "periodic" or "exponential"
#' @param diffTol Tolerance of model covergence
#'
#' @return beta and sigma pars
#' @examples
#' df = subset(simuData, Cluster==1)
#' mkt = subset(market, Cluster==1)
#' mkt = subset(mkt, select=-Cluster)
#'
#' aggrFit = aggrmodel(data = df, market = mkt, Y = Load, timeVar = Time, groupVar = Group, repVar = Rep, n_basis = 7)

aggrmodel <- function(formula=NULL,
                      data,
                      market,
                      Y = NULL,
                      timeVar,
                      groupVar,
                      repVar,
                      n_basis,
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
                         n_order = n_order)
    bb <- lapply(XList,
                 function(x)
                     do.call(rbind,replicate(n=I,
                                             expr=x,
                                             simplify=FALSE)
                             )
                 )
    bb <- do.call(rbind, bb)
    bb <- cbind(dd, bb)
    bb <- bb[,-c(1:3)]
    fit_init <- lm(y~.-1, data = bb)
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
        ## Get mean curves
        require(fda,quietly=TRUE)
        basisObj = create.bspline.basis(range(t),
                                        nbasis = n_basis,
                                        norder = n_order)
        B = predict(basisObj, t)
        mc <- tapply(beta_init,
                     rep(1:C,
                         each=n_basis),
                     function(b) B %*% b)

        ## <-- Under construction


        sigPar <- rep(sigma_init,C)
        corPar <- rep(-180, C)
        tauPar <- rep(.5, C)
        funcVecIn <- fit

        ## -->

    }

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
                     method = "L-BFGS-B")
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
        if(cicleRep){
        ## Note: Here we use a(0) = a(T) restriction
        } else {
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
            betaOut <- Matrix::solve(betaLeft, betaRight)
        } ## end if/else cicle

        ## W.4 UPDATES
        lkDiff <- abs(lkOut - lkIn)
        betaIn <- as.numeric(betaOut)
        parIn  <- parOut
        lkIn <- lkOut
    } ## End while(lkDiff > diffTol)

    ## Output Mean curves
    basisObj = create.bspline.basis(range(t),
                                    nbasis = n_basis,
                                    norder = n_order)
    B = predict(basisObj, t)
    ## Separate betas
    betaMtx <- cbind(beta=as.matrix(betaOut),
                     type=rep(1:C, each=n_basis))
    mcMtx <- tapply(betaMtx[,1],
                    betaMtx[,2],
                    function(x) B %*% x)
    mcMtx <- data.frame(mc=unlist(mcMtx),
                        time=rep(t,times=C),
                        type=rep(unique(market[,2]), each=length(t)))

    ## Return
    return(list('beta' = as.matrix(betaOut),
                'pars' = parOut,
                'mc' = mcMtx,
                'n_basis' = n_basis,
                'n_order' = n_order))
}


