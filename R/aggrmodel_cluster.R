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
  #
  # fitList <- tapply(X=trial_dd$cluster,INDEX = trial_dd$trial,
  #               FUN = function(x){
  #                 ## Merge datasets
  #                 tmp <- data.frame(group=grps,
  #                                   cluster=x)
  #
  #                 mktc <- merge(market, tmp)
  #                 mktsplit <- split(mktc, mktc$cluster)
  #                 ddc <- merge(data,tmp)
  #                 ddsplit <- split(ddc,ddc$cluster)
  #                 ## Apply fit
  #
  #
  #
  #                 fitListOut <- purrr::map2(ddsplit,
  #                                        mktsplit,
  #                                        ~simpleAggrmodel(.x,.y,
  #                                                         n_basis = K,
  #                                                         n_type = nt,
  #                                                         basisFunction = bFunc,
  #                                                         n_order = nord))
  #                 SR_list <- lapply(fitListOut, function(x) sum(abs(resid(x))))
  #                 Reduce('+',SR_list)
  #               })


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
  #
  # ## Return fit objects also
  # mktc <- merge(market, clusterOut)
  # mktsplit <- split(mktc, mktc$cluster)
  # ddc <- merge(data,clusterOut)
  # ddsplit <- split(ddc,ddc$cluster)
  # ## Apply fit
  # fitListOut <- purrr::map2(ddsplit,
  #                           mktsplit,
  #                           ~simpleAggrmodel(.x,.y,
  #                                            n_basis = n_knots,
  #                                            n_type = C,
  #                                            basisFunction = bFunc,
  #                                            n_order = n_order))
   return(list('clusterConfig'=clusterOut,
              'fitList' = fitListOut))
}

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
#' @param returnPred
#' @param optLk
#' @param corPar_init
#'
#' @return An aggrmodel_cluster object
#'
#' @name aggrmodel_cluster
#'
#' @import optimx
#' @importFrom mvtnorm dmvnorm
#' @export
#' @examples library(dplyr)
#'
#' data = simuData %>%
#'   select(group, rep, time, y)
#' df=data
#' mkt = market
#' mkt$Cluster=NULL
#' colnames(mkt) = c('group','type','num')
#'
#' ## Go for a walk after run the code below
#' fitCluster = aggrmodel_cluster(data = df, market=mkt, Y = 'y', timeVar =
#' 'time',groupVar = 'group', repVar = 'rep', n_basis = 7,n_basis_cov = NULL,
#' n_cluster = 2,n_trials = 1000, n_order = 4, corType = 'periodic', verbose=TRUE)
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
                              returnPred = FALSE,
                              optLk = TRUE,
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
  sigPar_init <- unlist(lapply(cl_init[[2]],function(x) summary(x)$sigma))
  sigPar_init <- matrix(sqrt(sigPar_init), ncol=B)
  if(is.null(corPar_init))
    corPar_init <- matrix(20, ncol=B)
  else
    corPar_init <- matrix(corPar_init, ncol=B)
  if(covType == 'Homog_Uniform'){
    cp_in <- c(sigPar_init,corPar_init)
  }
  if(covType == 'Homog'){
    cp_in <- c(rep(sigPar_init,C),rep(corPar_init,3))
  }
  pi_init <- prop.table(table(cluster_init[,2]))
  pi_init <- matrix(pi_init, ncol=B)
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
    probTab <- data.frame(
    reps = rep(1:I, times=J),
    grps =  rep(1:J, each=I)
    )
    probTab <- cbind(probTab, matrix(nrow=I*J,ncol=B))
    for(b in 1:B){
      tmp <- by(data = dd,INDICES = grp_rep,FUN = function(df){
        j <- df$group[1]
        i <- df$rep[1]
        mvtnorm::dmvnorm(x = as.numeric(df$y),
                         mean = as.numeric(c(xbeta[[j]][,b])),
                         sigma = as.matrix(sigMtxList[[b]][[j]]),
                         log=TRUE)
      })
      probTab[,c(2+b)] <- pi_init[b]*exp(as.numeric(tmp))
    }
    denom <- apply(probTab[,-c(1:2)],1,sum)
    for(b in 1:B)  probTab[,c(2+b)] <- probTab[,c(2+b)]/denom
    ## M-Step -----------------------------------------------------------------
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
            prb <- subset(probTab,subset = reps==i & grps==j)[,c(2+b)]
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
                            corWrap=corType
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
                  corWrap=corType
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
      # sigLongList <- lapply(sigMtxList_out,
      #                       function(m){
      #                         mtx <- lapply(m,function(x) Matrix::bdiag(replicate(n=I,x,simplify=FALSE)))
      #                         mtx <- Matrix::bdiag(mtx)
      #                         mtx
      #                       })
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
      # sigLongList <- lapply(sigMtxList_out,
      #                       function(m){
      #                         mtx <- lapply(m,function(x) Matrix::bdiag(replicate(n=I,x,simplify=FALSE)))
      #                         mtx <- Matrix::bdiag(mtx)
      #                         mtx
      #                       })
    }
    # probLong <- list()
    # probTab2 <- probTab[order(probTab[,2], probTab[,1]),]
    #     for(b in 1:B){
    #       probLong[[b]] <- diag(rep(probTab2[,b+2],each=T))
    #     }
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
          pijb <- subset(probTab,subset = reps==i & grps==j)[,c(2+b)]
          yij <- matrix( dd[dd$group==j&dd$rep==i,"y"], ncol=1)
          tmpLeft <- tmpLeft + pijb*mainLeft
          tmpRight <- tmpRight + pijb*(x_sinv %*% yij)

          bij <- pijb*x_sinv
          se <- se + bij %*% sigMtxList_out[[b]][[j]]%*% t(bij)

      # tmpLeft <- t(XLong)%*%probLong[[b]]%*%
      #   Matrix::solve(sigLongList[[b]],XLong)
      # tmpRight <- t(XLong)%*%probLong[[b]]%*%
      #   Matrix::solve(sigLongList[[b]],matrix(y,ncol=1))

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
    pi_out <- apply(probTab[,-c(1:2)],2,mean)
    xbeta_out <- lapply(X, function(x)  x %*% beta_out)

    ## Check convergence & updates
    # if(!optLk){
    #   lk_part <- matrix(nrow=I*J,ncol=B)
    #    for(b in 1:B){
    #      tmp <- by(data = dd,INDICES = grp_rep,FUN = function(df){
    #        j <- df$group[1]
    #        i <- df$rep[1]
    #        mvtnorm::dmvnorm(x = as.numeric(df$y),
    #                         mean = as.numeric(c(xbeta_out[[j]][,b])),
    #                         sigma = as.matrix(sigMtxList_out[[b]][[j]]),
    #                         log=FALSE)
    #      })
    #      lk_part[,b] <- pi_out[b]*as.numeric(tmp)
    #    }
    #   lk <- as.numeric(apply(lk_part,1,sum))
    #   lkOut <- sum(log(lk))
    #   if(verbose) message("lk value:",lkOut)
    # }
    lkDiff <- abs(lkIn - lkOut)
    beta_init <- beta_out
    pi_init <- matrix(pi_out,ncol=B)
    cp_in <- cp_out
    n_it <- n_it+1
    lkIn <- lkOut
  } # end while loop
  ## OUTPUTS ----
  predList <- apply(beta_out,2,function(b) XLong %*% b)
  predList <- unlist(predList)
  predList <- matrix(predList, ncol=B)
  predListNames <- paste("cluster",1:B,sep='')
  colnames(predList) <- predListNames
  ddOut <- cbind(dd,predList)
  probTab <- as.data.frame(probTab)
  colnames(probTab)[1:2] <- c('rep','group')
  ddOut <- merge(ddOut, probTab)
  ec <- ncol(ddOut)
  colnames(ddOut)[(ec-B+1):ec] <- paste('p',1:B,sep='')
  mm<- apply(ddOut[,(ec-B+1):ec],1,which.max)
  pred <- numeric(nrow(ddOut))
  for(k in 1:nrow(ddOut)){
    pred[k]<-ddOut[k,(mm[k]+4)]
  }
  ddOut$pred <- pred
  ddOut <- ddOut[,c(1:4,ncol(ddOut))]
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
  # betaSE <- matrix(nrow=nrow(beta_out),
  #                  ncol=B)
  # for(b in 1:B){
  #   sigmaum <- sigLongList[[b]]
  #   betaMult <- solve(tmpLeft, t(XLong)%*% probLong[[b]])%*%solve(sigmaum)
  #   betaSE[,b] <- sqrt(diag(betaMult%*%sigmaum%*%t(betaMult)))
  # }
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
                  "piPar" = pi_out[-c(1:2)],
                  "covPar" = covPar,
                  'mc' = mc)
  if(returnPred) outList[['predData']] <- ddOut
  class(outList)='aggrmodel_cluster'
  outList
}
