#' Fit a simple aggregated model
#' @name simpleAggrmodel
#'
#' @param market Market data frame. MUST be a 3 column dataframe with the following order: Group, Type and Number of subjects
#' @param n_basis  Number of basis functions for basis expansion
#' @param n_order Order of basis Splines (Default: 4)
#' @param basisFunction Character indicating which basis: 'B-Splines' (default) or 'Fourier'
#' @param data Dataset containing group, replicates (if any), time and aggregated signal
#'
#' @return lm object
#' @export
simpleAggrmodel <- function(data,
                            market,
                            n_basis,
                            n_order,
                            basisFunction){
  XList <- buildX(market=market,
                  timeVec = data$time,
                  n_basis = n_basis,
                  n_order = n_order,
                  basis = basisFunction)
  I <- length(unique(data$rep))
  X <- lapply(XList,
              function(x)
                do.call(rbind,replicate(n=I,
                                        expr=x,
                                        simplify=FALSE)
                )
  )
  X <- do.call(rbind, X)
  X <- cbind(y=data$y, X)
  X <- as.data.frame(X)
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
#' @export
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
  ## TODO
  ## [] Limit max n_trials

  ## Separate dataset in clusters
  C <- length(unique(market[,2]))
  J <- length(unique(data$group))
  smpObl <- rep(1:n_cluster, each=C)
  if((J-length(smpObl)) <0)
    stop('Cannot fit clusters if number of groups at each cluster are smaller than number of clusters')
  smpVec <- replicate(n = n_trials,
                      expr = sample(
                        c(smpObl, sample(1:n_cluster, size=J-length(smpObl),replace=TRUE))
                        ),
                      simplify = TRUE
  )
  grps <- unique(data$group)
  trial_dd <- data.frame(
    group = rep(grps, times = n_trials),
    cluster= c(smpVec),
    trial=rep(1:n_trials, each = length(grps))
  )
  fitList <- tapply(X=trial_dd$cluster,INDEX = trial_dd$trial,
                FUN = function(x,dd,mkt,K,nord,bFunc){
                  ## Merge datasets
                  tmp <- data.frame(group=grps,
                                    cluster=x)
                  mktc <- merge(mkt, tmp)
                  mktsplit <- split(mktc, mktc$cluster)
                  ddc <- merge(dd,tmp)
                  ddsplit <- split(ddc,ddc$cluster)
                  ## Apply fit
                  fitListOut <- purrr::map2(ddsplit,
                                         mktsplit,
                                         ~simpleAggrmodel(.x,.y,
                                                          n_basis = K,
                                                          basisFunction = bFunc,
                                                          n_order = nord))
                  sumList <- lapply(fitListOut, function(x) summary(x)$sigma)
                  Reduce('+',sumList)
                },
                ## tapply args
                dd=data,  mkt=market,  K=n_knots,  nord=n_order,  bFunc=bFunc)
  ## Choose smaller error and output the selected cluster config
  minErr <- as.integer(which.min(fitList))
  clusterOut <- subset(trial_dd, trial == minErr)
  clusterOut$trial <- NULL
  rownames(clusterOut) <- NULL
  ## Return fit objects also
  mktc <- merge(market, clusterOut)
  mktsplit <- split(mktc, mktc$cluster)
  ddc <- merge(data,clusterOut)
  ddsplit <- split(ddc,ddc$cluster)
  ## Apply fit
  fitListOut <- purrr::map2(ddsplit,
                            mktsplit,
                            ~simpleAggrmodel(.x,.y,
                                             n_basis = n_knots,
                                             basisFunction = bFunc,
                                             n_order = n_order))
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
#'
#' @return An aggrmodel_cluster object
#'
#' @name aggrmodel_cluster
#'
#' @importFrom mvtnorm dmvnorm
#' @export
#' @examples library(dplyr)
#'
#' data = simuData %>%
#'   select(group=Group, rep=Rep, time=Time, y=Load)
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
                              ##cicleRep = FALSE,
                              n_order = 4,
                              covType = 'Homog_Uniform',
                              corType = 'periodic',
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
  B = n_cluster
  K = n_basis
  dd <- data.frame(group=grps,
                   rep=reps,
                   time=t,
                   y = y)
  dd <- dd[order(dd[,2], dd[,1], dd[,3]),]
  y <- dd[,4]

  ## Set all initials parameters
  cl_init <- getClusterInitials(data=dd,market=market,
                                n_knots=n_basis,n_order=n_order,
                                n_trials=n_trials,n_cluster=B,
                                bFunc = basisFunction)
  cluster_init <- cl_init[[1]]
  beta_init <- unlist(lapply(cl_init[[2]],coef))
  beta_init <- matrix(beta_init, ncol=B)
  sigma_init <- unlist(lapply(cl_init[[2]],function(x) summary(x)$sigma))
  sigma_init <- matrix(sqrt(sigma_init), ncol=B)
  theta_init <- matrix(10, ncol=B)
  if(covType == 'Homog_Uniform'){
    cp_in <- c(sigma_init,theta_init)
  }
  if(covType == 'Homog'){
    cp_in <- c(rep(sigma_init,C),rep(theta_init,3))
  }
  pi_init <- prop.table(table(cluster_init[,2]))
  pi_init <- matrix(pi_init, ncol=B)
  ## While loop ------
  lkIn <- Inf
  lkDiff = diffTol+1
  n_it = 1
  XList <- buildX(market = market,timeVec = t,n_basis = K,
                  basis = basisFunction,n_order = n_order)
  X <- XList

  while(lkDiff > diffTol & n_it < itMax){
    if(verbose)
      message(paste("\nIteration num ",n_it))
    ## Compute prob for E-Step
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
    # sigMtxList <- lapply(1:B, function(b){
    #   sig <- covMatrix(market = market,group.name = 'group',
    #                    type.name = 'type',mkt.name = 'num',
    #                    timeVec = t,sigPar = sigma_init[1,b],
    #                    tauPar = NULL,corPar = theta_init[1,b],
    #                    funcMtx = NULL,covType = 'Homog_Uniform',
    #                    corType = 'periodic',nKnots = NULL,
    #                    truncateDec = 8)
    #   sig
    # })
    xbeta <- lapply(X, function(x) apply(beta_init, 2, function(bt) x %*% bt))
    probTab_init <- matrix(nrow=I*J, ncol=(B))
    probTab_names <- data.frame(reps = vector(mode=class(reps),length = I*J),
                                grps = vector(mode=class(grps),length = I*J))
    k=1
    for(i in unique(reps)){
      for(j in unique(grps)){
        for(b in 1:B){
          probTab_names$reps[k] = i
          probTab_names$grps[k] = j
          yij <- subset(dd, rep==i & group==j)$y
          probTab_init[k,b] <- mvtnorm::dmvnorm(x = as.numeric(yij),
                                                mean = as.numeric(c(xbeta[[j]][,b])),
                                                sigma = as.matrix(sigMtxList[[b]][[j]]),
                                                log=TRUE)
        }
        k=k+1
      }
    }
    for(k in 1:(I*J)){
      for(b in 1:B){
        probTab_init[k,b] <- probTab_init[k,b] + log(pi_init[1,b])
      }
    }
    denomProb <- log(rowSums(exp(probTab_init)))
    probTab_ratio=probTab_init ## initiate
    for(k in 1:(I*J)){
      for(b in 1:B){
        probTab_ratio[k,b] <- probTab_init[k,b] - denomProb[k]
      }
      if(all(exp(probTab_ratio[k,])==0) | all(is.infinite(exp(probTab_ratio[k,])))){ ## if all probs are zero
        b1 <- which.min(probTab_init[k,])
        probTab_ratio[k,] <- rep(0,times=B)
        probTab_ratio[k,b1] <- 1
      }else{
        probTab_ratio[k,] <- exp(probTab_ratio[k,])
      }
    }
    probTab <- cbind(probTab_names, probTab_ratio)
    ## M-Step
    # if(covType == 'Homog_Uniform')
    #   cp_in <- c(sigma_init,theta_init)
    # if(covType == 'Homog')
    #   cp_in <- c(rep(sigma_init,C),rep(theta_init,3))
    opt <- optim(par = cp_in,
                 fn = Q_wrapper,
                 #            lower=c(rep(-1e-60,B),rep(0,B)),
                 #             upper=c(rep(1e60,2*B)),
                 #             method = "L-BFGS-B",
                 data = dd,
                 market = market,
                 betaPar = beta_init,
                 piPar = pi_init,
                 pTab = probTab,
                 B=B,
                 t=t,
                 K=n_basis,
                 I=I,
                 J=J,
                 C=C,
                 basisFunction=basisFunction,
                 n_order=n_order,
                 covWrap=covType,
                 corWrap=corType
    )
    lkOut <- opt$value
    if(verbose) message(paste("\nlk value:", round(lkOut,6)))
    cp_out <-opt$par
    if(verbose) message(paste("\ncovPar estimates:",paste(round(opt$par,4),collapse=', ')))
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
      sigLongList <- lapply(sigMtxList_out,
                            function(m)
                              Matrix::bdiag(replicate(n=I,Matrix::bdiag(m))))
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
      sigLongList <- lapply(sigMtxList_out,
                            function(m)
                              Matrix::bdiag(replicate(n=I,Matrix::bdiag(m))))
    }
    probLong <- list()
        for(b in 1:B){
          probLong[[b]] <- diag(rep(probTab[,b+2],each=T))
        }
    ## COMPUTE BETAs
    XLong <- do.call(rbind,XList) ## pile by group
    XLong <- do.call(rbind, replicate(n=I,XLong, simplify = FALSE)) ## replicate piling
    beta_out <- beta_init ## initiate
    for(b in 1:B){
      # tmpLeft <- t(XLong)%*%XLong
      # tmpRight <- t(XLong)%*%matrix(y,ncol=1)
      tmpLeft <- t(XLong)%*%probLong[[b]]%*%
        Matrix::solve(sigLongList[[b]],XLong)
      tmpRight <- t(XLong)%*%probLong[[b]]%*%
        Matrix::solve(sigLongList[[b]],matrix(y,ncol=1))
      beta_out[,b] <- as.numeric(solve(tmpLeft,tmpRight))
    }
    pi_out <- apply(probTab,2,mean)
    ## Check convergence & updates
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
  ## Mean curves
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
  mc <- apply(beta_out,2,function(b) basisMtx%*%b)
  mc <- do.call(rbind,mc)
  mc <- data.frame(cluster = rep(paste('cluster',1:B,sep=''),each=C*T),
                   type = rep(rep(unique(market$type),each=T),times=B),
                   time = rep(tuni,times=C*B),
                   mc = as.numeric(mc))
  ## return -----
  outList <- list("probTab"=probTab,
                  "betaPar" = beta_out,
                  "piPar" = pi_out[-c(1:2)],
                  "covPar" = covPar,
                  'predData' = ddOut,
                  'mc' = mc)
  class(outList)='aggrmodel_cluster'
  outList
}
