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
#' @param corPar_init Numeric: Initial value for correlation parameters (default:20)
#' @param tauPar_init Numeric: Initial value for expoent parameters of complete covariance (default:0.5)
#' @param diffTol Tolerance of model covergence
#' @param verbose TRUE/FALSE prints likelihood values during optimization
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
                      timeVar2=NULL,
                      groupVar,
                      repVar,
                      n_basis,
                      n_basis_cov = NULL,
                      basisFunction = 'B-Splines',
                      timeVec2 = NULL,
                      n_basis2 = NULL,
                      n_order2 = NULL,
                      basisFunction2 = NULL,
                      cicleRep = FALSE,
                      n_order = 4,
                      covType = 'Homog_Uniform',
                      corType = 'periodic',
                      corPar_init = 20,
                      tauPar_init = .5,
                      optimMethod = "L-BFGS-B",
                      truncateDec = NULL,
                      verbose = FALSE,
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
    if(ncol(market) != 3)
        stop('Market must have 3 columns in the following order: Group, Type and Number of subjects. Please check your market input.')

    ## Preamble
    y = data[[substitute(Y)]]
    t = data[[substitute(timeVar)]]
    t = t/max(t)
    if(is.null(timeVar2)) t2 <- NULL
    else t2 <- data[[substitute(timeVar2)]]
    grps = as.factor(data[[substitute(groupVar)]])
    reps = as.factor(data[[substitute(repVar)]])
    J = length(unique(grps))
    I = length(unique(reps))
    T = length(unique(t))
    C = length(unique(market[,2]))
    dd <- data.frame(group=as.integer(grps),
                     rep=as.integer(reps),
                     time=t,
                     y = y)
    if(!is.null(timeVar2)) dd$time2 <- t2
    dd <- dd[order(dd$group, dd$rep, dd$time),]
    y <- dd$y
    colnames(market) <- c("group","type","num")
    market$group <- as.integer(factor(market$group,
                                      levels=levels(grps)))
    ## Build disagragation basis expansion design matrix
    mktLong <- spread(market,type,num)
    groupnrep <- subset(dd, select=c(group,rep))
    mktLong <- merge(groupnrep,mktLong)
    X <- buildX(market=mktLong,
                    nType = C,
                    timeVec = t,
                    n_basis = n_basis,
                    n_order = n_order,
                    basis = basisFunction,
                    timeVec2 = t2,
                    n_basis2 = n_basis2,
                    n_order2 = n_order2,
                    basis2 = basisFunction2)
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
    ## Get initial values -------------------------------------
    init <- get_inits(X=X, I=I, y=y, C=C,
                      covType=covType,
                      market=market,
                      corPar_init=corPar_init,
                      tauPar_init=tauPar_init,
                      truncateDec=truncateDec,
                      n_basis=n_basis, n_basis_cov=n_basis_cov,
                      t=t, n_order=n_order, basisFunc=basisFunc
                      )
    # X <- init$X
    beta_init <- init$beta
    parIn <- init$par
    lowerBoundVec <- init$lb
    upperBoundVec <- init$ub
    n_basis_cov <- init$n_basis_cov

    ## While ----------------------------------------
    betaIn = beta_init
    lkDiff = diffTol + 1
    lkIn = Inf
    while(lkDiff > diffTol){
        ## W.1 Obtain covariance estimates via optim
        opt <- optim(par = parIn,
                     fn = loglikWrapper,
                     covWrap = covType,
                     corWrap = corType,
                     dataWrap = dd,
                     mktWrap = market,
                     betaWrap = betaIn,
                     designWrap = X,
                     nCons = C,
                     lower = lowerBoundVec,
                     upper = upperBoundVec,
                     method = optimMethod,
                     nBasisCov = n_basis_cov,
                     nOrderCov = n_order,
                     truncateDec = truncateDec,
                     verbWrap = verbose)
        parOut <- opt$par
        lkOut <- opt$value

        ## W.2 Update Sigma estimates
        if(covType == 'Homog_Uniform'){
            sigmaOutList <- covMatrix(market = market,
                                      group.name = 'group',
                                      type.name = 'type',
                                      mkt.name = 'mkt',
                                      timeVec = dd$time,
                                      sigPar = parOut[1],
                                      corPar = parOut[2],
                                      covType = 'Homog_Uniform',
                                      corType = corType,
                                      truncateDec = truncateDec)
        }
        if(covType == 'Homog'){
            sigmaOutList <- covMatrix(market = market,
                                      group.name = 'group',
                                      type.name = 'type',
                                      mkt.name = 'mkt',
                                      timeVec = dd$time,
                                      sigPar = parOut[1:C],
                                      corPar = parOut[(C+1):(2*C)],
                                      covType = 'Homog',
                                      corType = corType,
                                      truncateDec = truncateDec)
        }
        if(covType == 'Heterog'){
            tvec <- unique(t)
            basisObj = create.bspline.basis(range(tvec),
                                            nbasis = n_basis_cov,
                                            norder = n_order)
            B <- predict(basisObj, tvec)
            betaMC <- parOut[1:(C*n_basis_cov)]
            betaMtx <- cbind(beta=as.matrix(betaMC),
                             type=rep(1:C, each=n_basis_cov))
            mcMtx <- tapply(betaMtx[,1],
                            betaMtx[,2],
                            function(x) B %*% x)
            funcVarIn <- matrix(unlist(mcMtx), ncol = C)
            sigParIn <- parOut[(C*n_basis_cov+1):(length(parOut)-(2*C))]
            corParIn <- parOut[(C*n_basis_cov+C+1):(length(parOut)-C)]
            tauParIn  <- parOut[((length(parOut)-C+1):length(parOut))]
            sigmaOutList <- covMatrix(market = market,
                                   group.name = 'group',
                                   type.name = 'type',
                                   mkt.name = 'mkt',
                                   timeVec = dd$time,
                                   funcMtx = funcVarIn,
                                   sigPar = sigParIn,
                                   corPar = corParIn,
                                   tauPar = tauParIn,
                                   covType = 'Heterog',
                                   corType = corType,
                                   truncateDec = truncateDec)
        }


        sigmaInvList <- lapply(sigmaOutList, qr.solve)

        ## W.3 Obtain beta estimates with previous steps
        sigmaInv <- lapply(sigmaInvList, function(x)
            bdiag(replicate(n=I,
                      expr=x,
                      simplify=FALSE))
            )
        sigmaInv <- bdiag(sigmaInv)
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

    ## Output ----------------------
    if(is.null(timeVar2)){
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
    }
    if(!is.null(timeVar2)){
        if(basisFunction=='B-Splines')
            basisObj = create.bspline.basis(range(t),
                                            nbasis = n_basis,
                                            norder = n_order)
        if(basisFunction=='Fourier'){
            basisObj = create.fourier.basis(range(t),
                                            nbasis = n_basis)
        }
        tuni <- unique(t)
        ## Second component
        n_basis2 <- ifelse(is.null(n_basis2),n_basis,n_basis2)
        basis2 <- ifelse(is.null(basisFunction2),basisFunction,basisFunction2)
        n_order2 <- ifelse(is.null(n_order2),n_order,n_order2)
        tuni2 <- unique(t2)
        if(basis2=='B-Splines')
            basisObj2 = create.bspline.basis(range(tuni2),
                                            nbasis = n_basis2,
                                            norder = n_order2)
        if(basis2=='Fourier'){
            basisObj2 = create.fourier.basis(range(tuni2),
                                            nbasis = n_basis2)
        }
        t1len <- length(tuni)
        t2len <- length(tuni2)
        t1 <- rep(tuni, times = t2len)
        t2 <- rep(tuni2, each = t1len)
        B1 = predict(basisObj, t1)
        B2 = predict(basisObj2, t2)
        nc <- ncol(B1)
        BB2 <- cbind(B1,B2)
        B <- apply(BB2,1,function(x) x[c(1:nc)] %x% x[-c(1:nc)])
        B <- t(B)
        rm(BB2)
        ## Separate betas
        L <- ncol(B2)
        betaMC <- betaOut[1:(C*L*n_basis)]
        betaMtx <- cbind(beta=as.matrix(betaMC),
                         type=rep(1:C, each=L*n_basis))
        mcMtx <- tapply(betaMtx[,1],
                        betaMtx[,2],
                        function(x) B %*% x)
        mcMtx <- data.frame(time=rep(t1,times=C),
                            time2=rep(t2,times=C),
                            type=rep(unique(market[,2]), each=length(t1)),
                            mc=unlist(mcMtx))
        # colnames(mcMtx) <- c(timeVar,timeVar2,type,mc)
    }
    ## Output data with predicted values
    dd$pred <- as.numeric(X %*% betaOut)
    dd$group <- as.factor(dd$group)
    levels(dd$group) <- levels(grps)
    dd$rep <- as.factor(dd$rep)
    levels(dd$rep) <- levels(reps)
    ## Return
    outList <- list('beta' = as.matrix(betaOut),
                'pars' = parOut,
                'mc' = mcMtx,
                'n_basis' = n_basis,
                'n_order' = n_order,
                'fitted' = dd)
    if(!is.null(formula))
        outList[['formula']] <- formula
    if(covType=='Heterog')
        outList["covCurves"] <- funcVarIn
    class(outList)='aggrmodel'
    return(outList)
}

##' Aux function to get initial values in aggrmodel function
##'
##' @title get_inits
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
##' @param basisFunc Basis function type
##' @return An object with initial values
##' @author Gabriel Franco
get_inits <- function(X, I, y, C,
                      covType,
                       market,
                       corPar_init, tauPar_init,
                       truncateDec,
                       n_basis, n_basis_cov,
                       t, n_order, basisFunc
                       ){

    # X <- lapply(XList,
    #             function(x)
    #                 do.call(rbind,replicate(n=I,
    #                                         expr=x,
    #                                         simplify=FALSE)
    #                         )
    #             )
    # X <- do.call(rbind, X)
    ## Get beta init -----------------------------
    ddfit_init <- data.frame(y=y,X)
    fit_init <- lm(y~.-1, data = ddfit_init)
    beta_init <- coef(fit_init)
    sigma_init <- sqrt(summary(fit_init)$sigma)
    rm(ddfit_init)

    if(covType == 'Homog_Uniform'){
        sigPar <- sigma_init
        corPar <- corPar_init
        parIn <- c(sigPar, corPar)
        lowBoundVec <- c(-Inf, 1e-8)
        ubCor <- ifelse(is.null(truncateDec), Inf, log(10^truncateDec))
        upperBoundVec <- c(Inf, ubCor)
    }

    if(covType == 'Homog'){
        sigPar <- rep(sigma_init,C)
        corPar <- rep(corPar_init, C)
        parIn <- c(sigPar, corPar)
        lowBoundVec <- c(rep(-Inf,C), rep(1e-20,C))
        ubCor <- ifelse(is.null(truncateDec), Inf, log(10^truncateDec))
        upperBoundVec <- c(rep(Inf,C), rep(ubCor,C))
    }
    if(covType == 'Heterog'){
        ## <----------------------- UNDER CONSTRUCTION
        if(is.null(n_basis_cov)){
            warning("Using the same number of basis of model fit:",n_basis)
            n_basis_cov <- n_basis
            betaCov_init <- beta_init
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
            betaCov_init <- coef(fit_init_cov)
        } # end if/else
        sigPar <- rep(sigma_init, C)/unlist(tapply(betaCov_init,
                                                   rep(1:C, each=n_basis_cov),FUN=mean))
        corPar <- rep(corPar_init, C)
        tauPar <- rep(tauPar_init, C)
        parIn <- c(betaCov_init,sigPar, corPar, tauPar)
        lowBoundVec <- c(rep(-Inf,times=(C*n_basis_cov+C)),
                         rep(0,2*C)) ## tau's > 0?
        ubCor <- ifelse(is.null(truncateDec), Inf, log(10^truncateDec))
        upperBoundVec <- c(rep(Inf,times=(C*n_basis_cov+C)),
                         rep(ubCor,C), rep(Inf,C)) ## tau's > 0?
    }# end if heterog

    output <- list()
    ## output$sigPar <- sigPar
    ## output$corPar <- corPar
    ## output$tauPar <- tauPar
    output$X <- X
    output$beta <- beta_init
    output$lb <- lowBoundVec
    output$ub <- upperBoundVec
    output$par <- parIn
    output$n_basis_cov <- n_basis_cov
    output
}
