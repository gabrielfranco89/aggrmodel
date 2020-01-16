#' Fit a simple aggregated model
#' @name simpleAggrmodel
#'
#' @param market
#' @param n_basis
#' @param n_order
#' @param basisFunction
#' @param data
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
#' @param data Must be an ordered data frame by rep > group > time with theses names
#' @param market A data frame with columns group, type and num
#' @param n_knots Number of knots to be fitted on simple model
#' @param n_trials Number of random clustering configurations to be fitted
#' @param bFunc Basis function to be used: "B-Splines" or "Fourier"
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
#' @param formula
#' @param data
#' @param market
#' @param Y
#' @param timeVar
#' @param groupVar
#' @param repVar
#' @param n_basis
#' @param n_basis_cov
#' @param basisFunction
#' @param n_order
#' @param covType
#' @param corType
#' @param diffTol
#' @param n_cluster
#' @param n_trials
#'
#' @name aggrmodel_cluster
#'
#' @importFrom mvtnorm dmvnorm
#' @export
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
  theta_init <- matrix(90, ncol=B)
  pi_init <- prop.table(table(cluster_init[,2]))
  pi_init <- matrix(pi_init, ncol=B)
  ## Get probabilities for first E-Step
  sigMtxList <- lapply(1:B, function(b){
    sig <- covMatrix(market = market,group.name = 'group',
                     type.name = 'type',mkt.name = 'num',
                     timeVec = t,sigPar = sigma_init[1,b],
                     tauPar = NULL,corPar = theta_init[1,b],
                     funcMtx = NULL,covType = 'Homog_Uniform',
                     corType = 'periodic',nKnots = NULL,
                     truncateDec = 8)
    sig
  })
  XList <- buildX(market = market,timeVec = t,n_basis = K,
                  basis = basisFunction,n_order = n_order)
  X <- XList
  xbeta <- lapply(X, function(x) apply(beta_init, 2, function(bt) x %*% bt))
  # yj <- split(y,grps)
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
  for(k in 1:(I*J)){
    for(b in 1:B){
      probTab_init[k,b] <- probTab_init[k,b] - denomProb[k]
    }
    if(all(exp(probTab_init[k,])==0)){ ## if all probs are zero
      b1 <- which.min(probTab_init[k,])
      probTab_init[k,] <- rep(0,times=B)
      probTab_init[k,b1] <- 1
    }else{
      probTab_init[k,] <- exp(probTab_init[k,])
    }
  }
  pTab_init <- cbind(probTab_names, probTab_init)

  ## While loop ------
  lkDiff = diffTol+1
  n_it = 1
  while(lkDiff > diffTol & n_it < itMax){
  cp <- c(sigma_init,theta_init)
  opt <- optim(cp,
               Q_wrapper,
               lower=c(rep(-Inf,B),rep(0,B)),
               data = dd,
               market = market,
               betaPar = beta_init,
               piPar = pi_init,
               B=B,
               t=t,
               K=n_basis,
               I=I,
               J=J,
               basisFunction=basisFunction,
               n_order=n_order
               )
  lkOut <- opt$value
  cp_out <- opt$par

  ## TODO
  # [] Compute Betas
  # [] Compute Pi's
  # [] Check convergence
  }

}
