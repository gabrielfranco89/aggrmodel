#' Fit aggregated model with clusters
#'
#' @param formula building...
#' @param data  Dataset containing group, replicates (if any), time and aggregated signal
#' @param market Market data frame. MUST be a 3 column dataframe with the following order: Group, Type and Number of subjects
#' @param Y Dependent variable: aggregated signal
#' @param timeVar Name of time variable
#' @param groupVar Name of grouping variable
#' @param repVar Name of replicates variable
#' @param n_basis Number of basis functions for basis expansion
#' @param n_basis_cov Number of basis functions for variance functional expansion
#' @param basisFunction Character indicating which basis: 'B-Splines' (default) or 'Fourier'
#' @param n_order Order of basis Splines (Default: 4)
#' @param covType Covariance functional type. One of "Homog_Uniform" (default), "Homog" or "Heterog"
#' @param corType Correlation structure type. One of "periodic" (default) or "exponential"
#' @param diffTol Tolerance of model covergence (Default: 1e-06)
#' @param n_cluster Number of grouping clusters
#' @param n_trials Number of random grouping trials for cluster initial values (Default: 42 and don't forget your towel!)
#' @param itMax Number of maximum iterations of EM algorithm (Default: 100)
#' @param verbose TRUE/FALSE indicating if steps of optimization should be printed as messages (Default: FALSE)
#' @param cicleRep Indicator TRUE/FALSE if replicates are cyclical
#' @param optMethod
#' @param returnFitted
#' @param optLk
#' @param corPar_init
#'
#' @return An aggrmodel_cluster object
#'
#' @name aggrmodel_cluster
#'
#' @import optimx
#' @importFrom mvnfast dmvn
#' @export
#' @examples
#'
#' set.seed(81453)
#' df <- createSimuData(nRep=5)
#' mkt <- attr(df,"market")
#'
#' fit_cl = aggrmodel_cluster(data = df, market=mkt, Y = obs, timeVar =
#' time,groupVar = group, repVar = rep, n_basis = 9, n_basis_cov = NULL,
#' n_cluster = 2,n_trials = 500, n_order = 4, corType = 'exponential',
#' returnFitted = TRUE, verbose=TRUE)
#'
#' plot(fitCluster, scales="free")
aggrmodel_cluster <- function(formula=NULL,
                              data,
                              market,
                              Y = NULL,
                              timeVar,
                              groupVar,
                              repVar,
                              n_basis,
                              n_basis_cov = NULL,
                              n_cluster = NULL,
                              n_trials = 42,
                              basisFunction = 'B-Splines',
                              cicleRep = FALSE,
                              n_order = 4,
                              covType = 'Homog_Uniform',
                              corType = 'exponential',
                              optMethod = "L-BFGS-B",
                              returnFitted = FALSE,
                              optLk = TRUE,
                              sigPar_init = NULL,
                              corPar_init = NULL,
                              diffTol = 1e-6,
                              itMax = 100,
                              verbose=FALSE){
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
  grps = as.factor(data[[substitute(groupVar)]])
  reps = as.factor(data[[substitute(repVar)]])
  J = length(unique(grps))
  I = length(unique(reps))
  T = length(unique(t))
  C = length(unlist(unique(market[,2])))
  B = n_cluster
  K = n_basis
  dd <- data.frame(group=as.integer(grps),
                   rep=as.integer(reps),
                   time=t,
                   y = y)
  dd <- dd[order(dd$group, dd$rep, dd$time),]
  y <- dd$y
  t <- dd$time
  colnames(market) <- c("group","type","num")
  market$group <- as.integer(factor(market$group,
                                    levels=levels(grps)))
  grp_rep <- rep(1:(I*J),each=T)
  dd_order <- data.frame(rep=rep(1:I,times=J),
                         group=rep(1:J,each=I),
                         id=1:(I*J)
                         )
  ## Set all initials parameters
  cl_init <- getClusterInitials(data=dd,market=market,
                                n_knots=n_basis,n_order=n_order,
                                n_trials=n_trials,n_cluster=B,
                                bFunc = basisFunction)
  cluster_init <- cl_init[[1]]
  beta_init <- unlist(lapply(cl_init[[2]],coef))
  beta_init <- matrix(beta_init, ncol=B)
  if(is.null(sigPar_init)){
    sigPar_init <- unlist(lapply(cl_init[[2]],function(x) summary(x)$sigma))
    sigPar_init <- rep(mean(sqrt(c(sigPar_init))), times=B)
    if(covType=="Homog") sigPar_init <- rep(sigPar_init,C)
  } else
    if(!all(length(sigPar_init)==B,length(sigPar_init==C*B)))
      stop(paste0("sigPar_init must have length ",B," or ",C*B))
  if(is.null(corPar_init)){
    corPar_init <- rep(1,times=B)
    if(covType=="Homog")
      covPar_init <- rep(corPar_init,times=C)
  }
  else
    if(!all(length(corPar_init)==B,length(corPar_init==C*B)))
      stop(paste0("corPar_init must have length ",B," or ",C*B))
  cp_in <- c(sigPar_init,corPar_init)
  # pi_init <- prop.table(table(cluster_init[,2]))
  # pi_init <- matrix(pi_init, ncol=B)
  pi_init <- matrix(nrow=J,ncol=B)
  for(b in 1:B){
    for(j in 1:J){
      pi_init[j,b]<-ifelse(cluster_init[j,2]==b,.9,(.1/(B-1)))
    }
  }
  ## While loop ------
  lkIn <- Inf
  lkDiff = diffTol+1
  n_it = 1
  mktLong <- tidyr::spread(market,type,num)
  groupList <- subset(dd, select=c(group,rep))
  mktLong <- merge(groupList,mktLong)
  mktList <- subset(mktLong, rep == as.integer(reps)[1])
  mktList <- subset(mktList, select = -rep)
  mktList <- split(mktList, mktList$group)
  X <- lapply(mktList, function(m){
    buildX(market = m, nType= C,
           timeVec = unique(t),n_basis = K,
           basis = basisFunction,n_order = n_order)
  }
  )
  mktLong$rep <- NULL
  XLong <-  buildX(market = mktLong, nType= C,
                   timeVec = t,n_basis = K,
                   basis = basisFunction,n_order = n_order)

  lkVec <- numeric(itMax)
  pi_out <- pi_init
  while(lkDiff > diffTol & n_it < itMax){
    if(verbose)
      message(paste("\nIteration num ",n_it))
    # E-Step --------------------------------------------------------------
    ## Compute prob for E-Step ----
    if(covType == 'Homog_Uniform'){
      covPar <- matrix(cp_in, nrow=B)
      sigMtxList <- lapply(1:B, function(b){
        sig <- covMatrix(market = market,group.name = 'group',
                         type.name = 'type',mkt.name = 'num',
                         timeVec = t,sigPar = covPar[b,1],
                         tauPar = NULL,corPar = covPar[b,2],
                         funcMtx = NULL,covType = 'Homog_Uniform',
                         corType = corType,nKnots = NULL,
                         truncateDec = 8)
        sig
      })
    }
    if(covType == 'Homog'){
      covPar <- matrix(cp_in, ncol=B,byrow=TRUE)
      sigMtxList <- lapply(1:B, function(b){
        sig <- covMatrix(market = market,group.name = 'group',
                         type.name = 'type',mkt.name = 'num',
                         timeVec = t,sigPar = covPar[1:C,b],
                         tauPar = NULL,corPar = covPar[(C+1):(2*C),b],
                         funcMtx = NULL,covType = 'Homog',
                         corType = corType,nKnots = NULL,
                         truncateDec = 8)
        sig
      })
    }
    xbeta <- lapply(X, function(x)  x %*% beta_init)
    ## Posteriori probability --------------------------------------------------
    probTab <- matrix(nrow=J,ncol=c(B+1))
    probTab[,1] <- 1:J
    dens_b <- matrix(nrow=I*J,ncol=B)
    for(b in 1:B){
      dens_b[,b] <- by(data = dd,INDICES = rep(1:c(I*J),each=T),FUN = function(dt){
        j <- dt$group[1]
        mvnfast::dmvn(
          X = dt$y,
          mu = xbeta[[j]][,b],
          sigma = sigMtxList[[b]][[j]],
          log = FALSE
        )
      }) # end by
    } # end for b
    my_constant <- 1/apply(dens_b,1,median)
    for(b in 1:B){
      dens_b[,b] <- dens_b[,b]*my_constant
      p_jb <- tapply(X = dens_b[,b],INDEX = rep(1:J,each=I),FUN = prod)
      numerator <- p_jb * pi_init[,b]
      probTab[,c(b+1)] <- numerator
    }
    denominator_j <- apply(X = probTab[,-1],MARGIN = 1,FUN = sum)
    for(b in 1:B) probTab[,c(b+1)] <- probTab[,c(b+1)]/denominator_j
    ## M-Step ==================================================================
    if(optLk){
      sigma_hat = NULL
    } else{
      pred <- XLong %*% beta_init
      sigma_hat <- list()
      for(b in 1:B){
        sigma_hat[[b]] <- list()
        dd$resid <- dd$y - pred[,b]
        tmp <- by(data = dd,INDICES = grp_rep,FUN = function(df){
          by(df, df$rep, function(x){
            out <- list()
            i <- x$rep[1]
            j <- x$group[1]
            rr <- matrix(x$resid,ncol=1)
#            prb <- subset(probTab,subset = reps==i & grps==j)[,c(2+b)]
            prb <- probTab[j,c(b+1)]
            out[[1]] <- prb * (rr %*%t (rr))
            out[[2]] <- prb
            out
          })
        })
        upr = lwr = 0
        for(j in 1:J){
          for(i in 1:I){
            id <- subset(dd_order,rep==i&group==j)$id
            upr <- upr + tmp[[id]][[1]]
            lwr <- lwr + tmp[[id]][[2]]
          } # end for i
          sigma_hat[[b]][[j]] <- upr/lwr
        } # end for j
      } # end for b
    } # end if optLk
    if(verbose) message("Optimizing...")
    if(optMethod=="bobyqa"){
      opt <- optimx::optimx(par = cp_in,
                            fn = Q_wrapper,
                            lower= rep(1e-4, length(cp_in)),
                            method = "bobyqa",
                            data = dd,
                            market = market,
                            optWrap = optLk,
                            piPar = pi_init,
                            pTab = probTab,
                            sCovList = sigma_hat,
                            B=B,
                            t=t,
                            K=n_basis,
                            xbeta =xbeta,
                            I=I,
                            J=J,
                            C=C,
                            basisFunction=basisFunction,
                            n_order=n_order,
                            covWrap=covType,
                            corWrap=corType,
                            hessian = TRUE
      )
    }
    else{
      opt <-optim(par = cp_in,
                  fn = Q_wrapper,
                  lower= rep(1e-4, length(cp_in)),
                  method = "L-BFGS-B",
                  data = dd,
                  market = market,
                  optWrap = optLk,
                  piPar = pi_init,
                  pTab = probTab,
                  sCovList = sigma_hat,
                  B=B,
                  t=t,
                  K=n_basis,
                  xbeta =xbeta,
                  I=I,
                  J=J,
                  C=C,
                  basisFunction=basisFunction,
                  n_order=n_order,
                  covWrap=covType,
                  corWrap=corType,
                  hessian = TRUE
      )
    }
    lkOut <- opt$value
    if(optMethod=="bobyqa")
      cp_out <- as.numeric(opt[1:length(cp_in)])
    else
      cp_out <-opt$par
    if(verbose) {
      message(paste("\nlk value:", round(lkOut,6)))
      message("parameters: ", paste(round(cp_out,4),collapse = ","))
    }
    ## Covariance matrices out ----
    if(covType=='Homog_Uniform'){
      covPar <- matrix(cp_out, nrow=B)
      sigMtxList_out <- lapply(1:B, function(b){
        sig <- covMatrix(market = market,group.name = 'group',
                         type.name = 'type',mkt.name = 'num',
                         timeVec = t,sigPar = covPar[b,1],
                         tauPar = NULL,corPar = covPar[b,2],
                         funcMtx = NULL,covType = 'Homog_Uniform',
                         corType = corType,nKnots = NULL,
                         truncateDec = 8)
        sig
      })
    }
    if(covType=='Homog'){
      covPar <- matrix(cp_out, ncol=B,byrow=TRUE)
      sigMtxList_out <- lapply(1:B, function(b){
        sig <- covMatrix(market = market,group.name = 'group',
                         type.name = 'type',mkt.name = 'num',
                         timeVec = t,sigPar = covPar[1:C,b],
                         tauPar = NULL,corPar = covPar[(C+1):(2*C),b],
                         funcMtx = NULL,covType = 'Homog',
                         corType = corType,nKnots = NULL,
                         truncateDec = 8)
        sig
      })
    }
    ## COMPUTE BETAs ----
    betaSE <- beta_out <- beta_init ## initiate
    for(b in 1:B){
      sigInvList <- lapply(sigMtxList_out[[b]], solve)
      tmpLeft <- matrix(0,nrow=nrow(beta_out),ncol=nrow(beta_out))
      tmpRight <- matrix(0,nrow=nrow(beta_out),ncol=1)
      se <- tmpLeft
      for(j in 1:J){
        x_sinv <- t(X[[j]]) %*% sigInvList[[j]]
        mainLeft <- x_sinv %*%  X[[j]]
        for(i in 1:I){
          pijb <- probTab[j,c(b+1)] #subset(probTab,subset = reps==i & grps==j)[,c(2+b)]
          yij <- matrix( dd[dd$group==j&dd$rep==i,"y"], ncol=1)
          tmpLeft <- tmpLeft + pijb*mainLeft
          tmpRight <- tmpRight + pijb*(x_sinv %*% yij)

          bij <- pijb*x_sinv
          se <- se + bij %*% sigMtxList_out[[b]][[j]]%*% t(bij)
        } # end for i
      } # end for j
      inv_tl <- solve(tmpLeft)
      se <- inv_tl %*% se %*% t(inv_tl)
      betaSE[,b] <- sqrt(diag(se))
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
        Phi = predict(basisObj, t)
        deltaVec <- matrix(Phi[1,] - Phi[nrow(Phi),],nrow=1)
        ## Get lambdas
        ## I can vectorize later (maybe)
        m1List <- split(x=qr.solve(tmpLeft),
                        f=rep(1:C, each=n_basis))
        m2List <- split(x=tmpRight,
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
        beta_out[,b] <- as.numeric(solve(tmpLeft, (tmpRight - t(Rvec))))
      } else {
        beta_out[,b] <- as.numeric(solve(tmpLeft,tmpRight))
      } ## end if/else cicle
    }

    pi_out <- probTab[,-1]
    # pi_out <- apply(probTab[,-c(1:2)],2,mean)
    xbeta_out <- lapply(X, function(x)  x %*% beta_out)
    ## Check convergence & updates
     lkDiff <- abs(lkIn - lkOut)
    beta_init <- beta_out
    pi_init <- pi_out
    cp_in <- cp_out
    lkVec[n_it] <- lkOut
    n_it <- n_it+1
    lkIn <- lkOut
  } # end while loop
  ## OUTPUTS ----
  predList <- XLong %*% beta_out
  pizao <- apply(pi_out,2,rep,each=T*I)
  pred <- numeric(nrow(dd))
  for(b in 1:B)
    pred <- pred + predList[,b]*pizao[,b]
  dd$pred <- pred
  dd$group <- as.factor(dd$group)
  levels(dd$group) <- levels(grps)
  dd$rep <- as.factor(dd$rep)
  levels(dd$rep) <- levels(reps)
  # dd <- dd[,c(1:4,ncol(dd))]
  ## Mean curves ----
  if(basisFunction=='B-Splines')
    basisObj = create.bspline.basis(range(t),
                                    nbasis = n_basis,
                                    norder = n_order)
  if(basisFunction=='Fourier'){
    basisObj = create.fourier.basis(range(t),
                                    nbasis = n_basis)
  }
  tuni <- unique(t)
  basisMtx = predict(basisObj, tuni)
  basisMtx <- Matrix::bdiag(replicate(C,basisMtx,simplify = FALSE))
  beta_lwr <- beta_out - qnorm(.975)*betaSE
  beta_upr <- beta_out + qnorm(.975)*betaSE
  mc_lwr <- basisMtx %*% beta_lwr
  mc_upr <- basisMtx %*% beta_upr
  mc <- basisMtx%*%beta_out
  mc <- data.frame(cluster = rep(paste('cluster',1:B,sep=''),each=C*T),
                   type = rep(rep(unique(market$type),each=T),times=B),
                   time = rep(tuni,times=C*B),
                   mc = mc@x,
                   mc_lwr = mc_lwr@x,
                   mc_upr = mc_upr@x)
  ## return -----
  outList <- list("probTab"=probTab,
                  "betaPar" = beta_out,
                  "betaSE" = betaSE,
                  "piPar" = pi_out,
                  "covPar" = covPar,
                  "covParSE" = tryCatch(sqrt(diag(solve(opt$hessian))), error=function(e) e),
                  'mc' = mc,
                  'lkVec' = lkVec)
  if(returnFitted) outList[['fitted']] <- dd
  class(outList)='aggrmodel_cluster'
  outList
}
