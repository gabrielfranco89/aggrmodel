#' Fit Aggregated Model
#'
#' @name aggrmodel
#'
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
#' @param returnFitted Should the fitted values be returned in output?
#' @param diffTol Tolerance of model covergence
#' @param itMax Maximum number of iterations (Default: 100)
#' @param verbose TRUE/FALSE prints likelihood values during optimization
#' @param timeVar2 Name of second functional
#' @param n_basis2 Number of basis for second functional
#' @param n_order2 Order for second functional expansion
#' @param basisFunction2 Character indicating which basis: 'B-Splines' or 'Fourier'
#' @param sigPar_init Inital values for sigma
#' @param betaCov_init Inital values for variance functional expansion
#' @param positive_restriction TRUE/FALSE if mean curves are strictly positive
#' @param optimMethod Choose optim method (Default: L-BFGS-B)
#' @param truncateDec Decimal to be truncated at covariance matrix
#' @param optSampleCovMatrix Optmization criterion via sample covariance matrix convergence (TRUE and default) or via likelihood (more sensitive)
#' @param optVerbose Print parameters while in optmization
#' @param useGrad Use gradient function approximation for optimization? (Default: FALSE)
#'
#' @return An aggrmodel object
#' @examples
#' df = simuData
#' mkt = attr(df, "market")
#' df = subset(df, cluster==2)
#' mkt = subset(mkt, group %in% unique(df$group))
#'
#' aggrFit = aggrmodel(data = df, market = mkt, Y = y, timeVar = time, groupVar = group, corType = 'exponential',repVar = rep, n_basis = 7)
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
                      n_basis2 = NULL,
                      n_order2 = NULL,
                      basisFunction2 = NULL,
                      cicleRep = FALSE,
                      n_order = 4,
                      covType = 'Homog_Uniform',
                      corType = 'exponential',
                      optSampleCovMatrix = FALSE,
                      sigPar_init = NULL,
                      corPar_init = NULL,
                      tauPar_init = NULL,
                      betaCov_init = NULL,
                      returnFitted = TRUE,
                      positive_restriction = FALSE,
                      optimMethod = "L-BFGS-B",
                      truncateDec = NULL,
                      verbose = FALSE,
                      optVerbose = FALSE,
                      useGrad = FALSE,
                      diffTol = 1e-5,
                      itMax = 100){


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
    C = length(unlist(unique(market[,2])))
    dd <- data.frame(group=as.integer(grps),
                     rep=as.integer(reps),
                     time=t,
                     y = y)
    if(!is.null(timeVar2)) dd$time2 <- t2
    dd <- dd[order(dd$group, dd$rep, dd$time),]
    y <- dd$y
    t <- dd$time
    colnames(market) <- c("group","type","num")
    market$group <- as.integer(factor(market$group,
                                      levels=levels(grps)))
    ## Build disagragation basis expansion design matrix
    mktLong <- tidyr::spread(market,type,num)
    grouprep <- subset(dd, select=c(group,rep))
    grouprep$flag <- rep(1:c(I*J), each=T)
    groupList <- subset(grouprep, select=group)
    mktLong <- merge(groupList,mktLong)
    X <- buildX(market=mktLong,
                    nType = C,
                    timeVec = dd$time,
                    n_basis = n_basis,
                    n_order = n_order,
                    basis = basisFunction,
                    timeVec2 = t2,
                    n_basis2 = n_basis2,
                    n_order2 = n_order2,
                    basis2 = basisFunction2)
    if(!is.null(formula)){
        cvtMtx <- as.data.frame(model.matrix(formula, data = data))
        cvtMtx <- cbind(cvtMtx,
                         group=grps,
                         rep=reps,
                         time=dd$time
                         )
        cvtMtx <- cvtMtx[order(cvtMtx[['group']],
                                 cvtMtx[['rep']],
                                 cvtMtx[['time']]),]
        X <- cbind(X,subset(cvtMtx,
                            select=-c(rep,group,time)))
        X <- as.matrix(X)
    } ## end if is.null(formula)
    XList <- split(X, grouprep$flag)
    XList <- lapply(XList, matrix, nrow=T)
    ## Get initial values -------------------------------------
    init <- get_inits(X=X, I=I,
                      data = dd,
                      #y=y,
                      C=C,
                      covType=covType,
                      market=mktLong,
                      sigPar_init=sigPar_init,
                      corPar_init=corPar_init,
                      #                    tauPar_init=tauPar_init,
                      betaCov_init=betaCov_init,
                      truncateDec=truncateDec,
                      n_basis=n_basis,
                      n_basis_cov=n_basis_cov,
                      t=t, n_order=n_order, basisFunc=basisFunction
                      )
    beta_init <- init$beta
    if(positive_restriction)
        beta_init <- solve(t(X)%*%X) %*% (t(X)%*%log(y))
    parIn <- init$par
    lowerBoundVec <- init$lb
    upperBoundVec <- init$ub
    n_basis_cov <- init$n_basis_cov

    ## While ----------------------------------------
    betaIn = beta_init
    lkDiff = diffTol + 1
    lkIn = Inf
    itCount = 1
    lkVec = normVec = numeric(itMax)
    if(!useGrad) gradLK <- NULL
    while(lkDiff > diffTol & itCount < itMax){
        ## W.1 Obtain covariance estimates via optim
        if(positive_restriction){
            parIn <- c(betaIn,parIn)
            pLen <- length(parIn)
            lowerBoundVec <- rep(-Inf, pLen)
            upperBoundVec <- rep(Inf, pLen)
        }
        if(!optSampleCovMatrix){
            residListGroup <- NULL
        } else {
            dd$resid <- dd$y - as.numeric(X %*% betaIn)
            residList <- by(dd, INDICES = dd$group,FUN = function(x){
                by(x, x$rep, function(d){
                    rr <- matrix(d$resid, ncol=1)
                    rr %*% t(rr)
                })
            })
            residListGroup <- lapply(residList, function(x){
                nRep <- length(x)
                Reduce("+",x) / (nRep - 1)
            }
            )
            rm(residList)
        } # end if optSCovMatrix

        opt <- optim(par = parIn,
                     fn = loglikWrapper,
                     gr = gradLK,
                     covWrap = covType,
                     corWrap = corType,
                     dataWrap = dd,
                     mktWrap = market,
                     betaWrap = betaIn,
                     sCovWrap = residListGroup,
                     designWrap = X,
                     optWrap = optSampleCovMatrix,
                     nCons = C,
                     lower = lowerBoundVec,
                     upper = upperBoundVec,
                     method = optimMethod,
                     nBasisCov = n_basis_cov,
                     nOrderCov = n_order,
                     truncateDec = truncateDec,
                     positive = positive_restriction,
                     verbWrap = optVerbose,
                     hessian=TRUE)
        parOut <- opt$par

        ## W.2 Update Sigma estimates ----
        if(positive_restriction){
            nCoef <- ncol(X)
            betaOut <- parOut[c(1:nCoef)]
            parOut <- parOut[-c(1:nCoef)]
            if(covType == 'Homog_Uniform')
                parOut[2] <- log(parOut[2])
            if(covType == 'Homog')
                parOut[(C+1):(2*C)] <- log(parOut[(C+1):(2*C)])
            if(covType == 'Heterog'){
                 parOut[(C*n_basis_cov+1):(length(parOut))] <- log( parOut[(C*n_basis_cov+1):(length(parOut))])
                 # parOut[((length(parOut)-C+1):length(parOut))] <- log( parOut[((length(parOut)-C+1):length(parOut))])
            }
            parsSE <-  tryCatch(sqrt(abs(diag(solve(opt$hessian)))),error=function(e) e )
            betaSE <- parsSE[1:nCoef]
        }
        if(!positive_restriction){
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
                funcVarIn <- B%*% matrix(betaMC,ncol=C)
                # betaMtx <- cbind(beta=as.matrix(betaMC),
                #                  type=rep(1:C, each=n_basis_cov))
                # mcMtx <- tapply(betaMtx[,1],
                #                 betaMtx[,2],
                #                 function(x) exp(B %*% x))
                # funcVarIn <- matrix(unlist(mcMtx), ncol = C)
                # sigParIn <- parOut[(C*n_basis_cov+1):(length(parOut)-(2*C))]
                corParIn <- parOut[(C*n_basis_cov+1):(length(parOut))]
                # tauParIn  <- parOut[((length(parOut)-C+1):length(parOut))]
                sigmaOutList <- covMatrix(market = market,
                                          group.name = 'group',
                                          type.name = 'type',
                                          mkt.name = 'mkt',
                                          timeVec = dd$time,
                                          funcMtx = funcVarIn,
                                          sigPar = rep(1,C),#sigParIn,
                                          corPar = corParIn,
                                          tauPar =rep(1,C), # tauParIn,
                                          covType = 'Heterog',
                                          corType = corType,
                                          truncateDec = truncateDec)
            }
            ## W.3 Obtain beta estimates with previous steps
            gr1 <- unique(grouprep)
            sigmaInvList <- lapply(sigmaOutList, solve)
            betaLeft <- matrix(0,nrow=ncol(X),ncol=ncol(X))
            betaRight <- matrix(0,nrow=ncol(X),ncol=1)
            se <- betaLeft
            for(k in unique(gr1$flag)){
                i <- gr1[gr1$flag==k,"rep"]
                j <- gr1[gr1$flag==k,"group"]
                x_sinv <- t(XList[[k]]) %*% sigmaInvList[[j]]
                mainLeft <- x_sinv %*%  XList[[k]]
                yij <- matrix( dd[dd$group==j&dd$rep==i,"y"], ncol=1)
                betaLeft <- betaLeft + mainLeft
                betaRight <- betaRight + x_sinv %*% yij
                se <- se + x_sinv %*% sigmaOutList[[j]]%*% t(x_sinv)
            }
            inv_tl <- solve(betaLeft)
            se <- inv_tl %*% se %*% t(inv_tl)
            betaSE <- sqrt(diag(se))
            # sigmaInv <- lapply(sigmaInvList, function(x)
            #     bdiag(replicate(n=I,
            #                     expr=x,
            #                     simplify=FALSE))
            # )
            # sigmaInv <- bdiag(sigmaInv)
            # betaLeft <- t(X)%*%sigmaInv%*%X
            # betaRight <- t(X)%*%sigmaInv%*%y
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
                Rvec <- matrix(B[1,] - B[nrow(B),],nrow=1)
                Rvec <- matrix(1,ncol=C) %x% Rvec
                if(!is.null(formula)) Rvec <- cbind(Rvec, matrix(0,ncol=c(ncol(X) - C*n_basis)))
                lambda <- solve(betaLeft, betaRight)
                lambda <- (Rvec %*% lambda) / (Rvec%*%t(Rvec))
                lambda <- as.numeric(lambda)
                ## Get lambdas
                # ## I can vectorize later (maybe)
                # m1List <- split(x=solve(betaLeft),
                #                 f=rep(1:C, each=n_basis))
                # m2List <- split(x=betaRight,
                #                 f=rep(1:C, each=n_basis))
                # lambda <- matrix(numeric(C),nrow=1)
                # for(c in 1:C){
                #     lUp <- deltaVec%*%
                #         matrix(m1List[[c]],
                #                nrow=n_basis)[,(n_basis*(c-1)+1):(c*n_basis)]%*%
                #         matrix(m2List[[c]],ncol=1)
                #     lDown <- deltaVec%*%
                #         matrix(m1List[[c]],
                #                nrow=n_basis)[,(n_basis*(c-1)+1):(c*n_basis)]%*%
                #         t(deltaVec)
                #     lambda[1,c] <- lUp/lDown
                # }
                # Rvec <- lambda %x% deltaVec
                betaOut <- solve(betaLeft, (betaRight - lambda*t(Rvec)))
            } else {
                betaOut <- solve(betaLeft, betaRight)
            } ## end if/else cicle
        } ## end if !positive_restriction
        ## W.4 UPDATES ----
        betaIn <- as.numeric(betaOut)
        parIn  <- parOut
        if(optSampleCovMatrix){
            lkValue <- logLikelihood(data = dd,
                                   muVec = X %*% betaIn,
                                   covMtxList = sigmaOutList
            )
            normVec[itCount] <- opt$value
        }
        else{
            lkValue <- opt$value
        }
        lkOut <- opt$value
        if(verbose){
            if(optSampleCovMatrix) message("norm value: ",lkOut)
            else message("lk value: ",lkOut)
        }
        lkVec[itCount] <- lkValue
        lkDiff <- abs(lkOut - lkIn)
        #if(positive_restriction) lkDiff = 0
        lkIn <- lkOut
        itCount <- itCount + 1
        if(itCount == itMax)
            warning(paste("Loop stopped at max iteration count:",itMax))
    } ## End while(lkDiff > diffTol)

    ## Output ----------------------
    # ## Beta Std Error ----
    # sigmaum <- lapply(sigmaOutList, function(x)
    #     bdiag(replicate(n=I,
    #                     expr=x,
    #                     simplify=FALSE))
    # )
    # sigmaum <- bdiag(sigmaum)
    # betaMult <- solve(t(X)%*%sigmaInv%*%X) %*% t(X)%*%sigmaInv
    # betaSE <- sqrt(diag(betaMult%*%sigmaum%*%t(betaMult)))
    ## Mean curves ----
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
        betaMC_SE <- betaSE[1:(C*n_basis)]
        betaLwr <- betaMC - qnorm(.975)*betaMC_SE
        betaUpr <- betaMC + qnorm(.975)*betaMC_SE
        mc <- B %*% matrix(betaMC, ncol=C)
        mc_lwr <- B %*% matrix(betaLwr, ncol=C)
        mc_upr <- B %*% matrix(betaUpr, ncol=C)
        mcMtx <- data.frame(mc= as.numeric(mc),# unlist(mcMtx),
                            mc_lwr = as.numeric(mc_lwr),
                            mc_upr = as.numeric(mc_upr),
                            time=rep(tuni,times=C),
                            type=rep(unlist(unique(market[,2])), each=length(tuni)))
        if(positive_restriction){
          mcMtx[,1:3] <- exp(mcMtx[,1:3])
        }
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
        # betaMtx <- cbind(beta=as.matrix(betaMC),
        #                  type=rep(1:C, each=L*n_basis))
        # mcMtx <- tapply(betaMtx[,1],
        #                 betaMtx[,2],
        #                 function(x) B %*% x)
        mcMtx <- B%*%matrix(betaMC,ncol=C)
        mcMtx <- data.frame(time=rep(t1,times=C),
                            time2=rep(t2,times=C),
                            type=rep(unique(market[,2]), each=length(t1)),
                            mc=c(mcMtx))
        # colnames(mcMtx) <- c(timeVar,timeVar2,type,mc)
    }
    ## Output data with predicted values
    dd$pred <- as.numeric(X %*% betaOut)
    if(positive_restriction) dd$pred <- exp(dd$pred)
    dd$group <- as.factor(dd$group)
    levels(dd$group) <- levels(grps)
    dd$rep <- as.factor(dd$rep)
    levels(dd$rep) <- levels(reps)
    ## Return ----
    parsSE <-  tryCatch(sqrt(abs(diag(solve(opt$hessian)))),error=function(e) e )
    outList <- list('beta' = as.matrix(betaOut),
                'pars' = parOut,
                'parsSE' = ifelse(positive_restriction, parsSE[-c(1:nCoef)], parsSE),
                'mc' = mcMtx,
                'n_basis' = n_basis,
                'n_order' = n_order,
                'lkVec' = lkVec,
                'betaSE' = ifelse(positive_restriction, parsSE[1:nCoef], betaSE))
    if(returnFitted) outList[["fitted"]] <- dd
    else outList[["residuals"]] <- dd$resid
    if(optSampleCovMatrix) outList[["normValues"]] <- normVec
    if(!is.null(formula))
        outList[['formula']] <- formula
    if(covType=='Heterog'){
        tvec <- unique(dd$time)
        basisObj = create.bspline.basis(range(tvec),
                                        nbasis = n_basis_cov,
                                        norder = n_order)
        B <- predict(basisObj, tvec)
        betaCov <- parOut[1:(C*n_basis_cov)]
        betaCov_lwr <- betaCov - qnorm(.975)*outList[["parsSE"]][1:(C*n_basis_cov)]
        betaCov_upr <- betaCov + qnorm(.975)*outList[["parsSE"]][1:(C*n_basis_cov)]
        ddCov <- data.frame(funcVar = c( B %*% matrix(betaCov, ncol=C)),
                            funcVar_lwr = c( B %*% matrix(betaCov_lwr, ncol=C)),
                            funcVar_upr = c( B %*% matrix(betaCov_upr, ncol=C)),
                            type = rep(1:C, each=length(tvec)),
                            time = rep(tvec,times=C)
        )
        outList[["covCurves"]] <- ddCov
    }
    class(outList)='aggrmodel'
    return(outList)
}

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
    # y <- data$y
    # ddfit_init <- data.frame(y=y,X)
    fit_init <- lm(data$y~X-1)
    beta_init <- coef(fit_init)
    sigma_fit <- sqrt(summary(fit_init)$sigma/(I-1))
    # rm(ddfit_init)

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
        ## fit new object for heterog ----
            # data$resid <- resid(fit_init)
            # # sCovDiag <- by(data, INDICES = data$group,FUN = function(x){
            # #     by(x, x$rep, function(y){
            # #         rr <- matrix(y$resid, ncol=1)
            # #         rr %*% t(rr)
            # #     })
            # # })
            # # sCovDiag <- lapply(sCovDiag, function(x){
            # #     nRep <- length(x)
            # #     out <- Reduce("+",x) / (nRep -1)
            # #     diag(out)
            # # }
            # # )
            # # sCovDiag <- lapply(sCovDiag, function(x) do.call(c, replicate(n=I, x, simplify = FALSE)))
            # # sCovDiag <- unlist(sCovDiag)
            # Xcov <- buildX(market=market,
            #                nType = C,
            #                timeVec = t,
            #                n_basis = n_basis_cov,
            #                n_order = n_order,
            #                basis = basisFunc)
            # # Xcov <- lapply(XListCov,
            # #             function(x)
            # #                 do.call(rbind,replicate(n=I,
            # #                                         expr=x,
            # #                                         simplify=FALSE)
            # #                         )
            # #             )
            # # Xcov <- do.call(rbind, Xcov)
            # # Xcov <- data.frame(sVar = log(sCovDiag), Xcov)
            # fit_init_cov <- lm(I(data$resid^2)~Xcov-1)#, data = Xcov)
            # beta1 <- coef(fit_init_cov)
            # # bObj <- fda::create.bspline.basis(0:1,nbasis = n_basis_cov,
            # #                                   norder = n_order,
            # #                                   basis = basisFunc)
            # # bMtx <- predict(bObj, unique(data$time))
            # # nu_hat <- Matrix::bdiag(replicate(n=C,bMtx,simplify = FALSE)) %*% beta1
            # # nu_hat <- matrix(sqrt(nu_hat), ncol = C)
            # # betaCov_init <- apply(nu_hat,2,function(nn){
            # #     dd_tmp <- data.frame(y=nn, bMtx)
            # #     coef(lm(y~.-1, data =dd_tmp))
            # # })
            # # betaCov_init <- c(betaCov_init)
            # betaCov_init <- beta1
            betaCov_init <- rep(sigma_fit, C*n_basis_cov)
        }
        else{
            if(length(betaCov_init)!=C*n_basis_cov) stop("betaCov_init must have the same length as number of basis for covariance")
        } # end if/else
        # if(is.null(sigPar_init))
        #     sigPar <- rep(1,C)
        #     # sigPar <- rep(sigma_fit, C)/unlist(tapply(betaCov_init,
        #     #                                           rep(1:C, each=n_basis_cov),FUN=mean))
        # else
        #     sigPar <- sigPar_init
        if(is.null(corPar_init))
            corPar <- rep(1, C)
        else
            corPar <- corPar_init
        # if(is.null(tauPar_init))
        #     tauPar <- rep(1, C)
        # else
        #     tauPar<- tauPar_init
        parIn <- c(betaCov_init,
                   # sigPar,
                   # tauPar,
                   corPar)
        lowBoundVec <- c(rep(-Inf,times=(C*n_basis_cov)), ##beta
                         rep(1e-4,C), ## sigPar & corPar
                         rep(0,C)) ## tauPar
        ubCor <- ifelse(is.null(truncateDec), Inf, log(10^truncateDec))
        upperBoundVec <- c(rep(Inf,times=(C*n_basis_cov)),
                           rep(ubCor,C))
                           # rep(Inf,C)) ## tau's > 0?
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
