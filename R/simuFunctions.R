#' Simulate aggregated data
#' @name createSimuData
#'
#' @param B1 Number of groups from cluster 1
#' @param B2 Number of groups from cluster 2
#' @param nRep Number of replicates
#' @param market Optional. Insert the market as a data frame with the columns 'group', 'type' and 'num'
#' @param market_range If market not inserted, the range to be sampled
#' @param sigPar Parameters of sigma: vector of length 3
#' @param tauPar Parameters of tau: vector of length 3
#' @param corPar Parameters of correlation: vector of length 3
#' @param nu1 Functional of variance for cluster 1 (must be of length N)
#' @param nu2 Functional of variance for cluster 2 (must be of length N)
#' @param seed Optional. Seed to be used for simulation
#' @param tempPar Default:1. If you do not want temperature effect, set temPar=0
#' @param tmp_bt1 Beta parameters for winter temperature simulation. Must be of size  5*J*nRep. If NULL,  rnorm(5*J*nRep,sd=2) will be used
#' @param tmp_bt2 Beta parameters for summer temperature simulation. Must be of size  5*J*nRep. If NULL,  rnorm(5*J*nRep,sd=2) will be used
#'
#' @importFrom splines bs
#'
#' @return
#' Return a data frame with nRep replicates and 3 types of consumer observed in 48 times. Parameters and market are returned as attributes.
#'
#' @examples
#' dd = createSimuData()
#' dd = createSimuData(B1=4, B2=8)
#' @export
createSimuData <- function(B1 = 4,
                           B2 = 6,
                           nRep = 20,
                           market=NULL,
                           market_range = 5:25,
                           sigPar = c(2,2,2,4,4,4),
                           tauPar = c(.5,.5,.5, .4,.4,.4),
                           corPar = c(4,4,4,6,6,6),
                           nu1 = NULL,
                           nu2 = NULL,
                           tmp_bt1 = NULL,
                           tmp_bt2 = NULL,
                           beta1 = 0,
                           beta2 = 0,
                           seed=NULL,
                           tempPar=1
                           ){
    ## Preamble
    if(any(B1<3,B2<3))
        stop("B1 and B2 must be greater than 3 for model identifiability!")
    J <- B1+B2
    C <- 3
    sigPar <- matrix(sigPar, ncol=2)
    tauPar <- matrix(tauPar, ncol=2)
    corPar <- matrix(corPar, ncol=2)
    if(!is.null(seed)) set.seed(seed)
    if(is.null(market)){
        mkt <- data.frame(group=rep(1L:J, each=C),
                          type=rep(1L:C,times=J),
                          num=sample(market_range,
                                     size=J*C,
                                     replace=TRUE))
    }else{
        mkt <- market
    }# end if/else market is null

    ## READ MEAN CURVES AND TEMPERATURES
    mc <- simulatedMeanCurves
    mc$Type <- rep(1:3,each=48)
    mc$Time2 <- mc$Time/24
    ## Create base dataset
    mktMtx <- matrix(mkt[,3], nrow=J, byrow=TRUE)
    colnames(mktMtx) <- paste('C',1:C,sep='')
    mktMtx <- as.data.frame(mktMtx)
    mktMtx$group=1:J
    b1 <- matrix(mc$Cluster1, ncol=C)
    b2 <- matrix(mc$Cluster2, ncol=C)
    ## colocar n_c1 e n_c2
    b <- rbind(
        do.call(rbind,replicate(n=B1,expr=b1, simplify=FALSE)),
        do.call(rbind,replicate(n=B2,expr=b2, simplify=FALSE))
    )
    dd <- data.frame(group = rep(1:J, each=48),
                     time=rep(x=(1:48)/48, times=J)
                     )
    dd <- merge(dd,mktMtx)
    dd <- cbind(dd, b)
    ## CREATE REPLICATES ----
    nr <- nrow(dd)
    dd <- do.call(rbind, replicate(n=nRep,
                                   expr=dd,
                                   simplify=FALSE))
    dd$rep <- rep(1:nRep, each=nr)
    if(tempPar != 0){
        tmp_x <- splines::bs(seq(0,1, length.out=nrow(dd)/2),
                             df = 5*J*nRep,
                             intercept = TRUE)
        if(is.null(tmp_bt1)) tmp_bt1 <- rnorm(5*J*nRep,sd=2)
        if(is.null(tmp_bt2)) tmp_bt2 <- rnorm(5*J*nRep,sd=2)
        dd$temp <- c(15+tmp_x %*% tmp_bt1, 24+tmp_x %*% tmp_bt2)
    }
    else
        dd$temp <- 10L
    ## ADD TEMPERATURE EFFECT ----
    dd$`1` <- dd$`1` * log(log(dd$temp))^tempPar
    dd$`2` <- dd$`2` * log(log(dd$temp))^tempPar
    dd$`3` <- dd$`3` * log(log(dd$temp))^tempPar
    ## ADD DUMMY VAR
    ##   added only on cluster 2
    dd$dummy1 <- ifelse(dd$group %in%  (B1+1):(B1+2),1,0)
    dd$dummy2 <- ifelse(dd$group == J, 1, 0)
    ## AGGREGATED DATA ----
    dd$signal <- dd$C1*dd$`1` +
        dd$C2*dd$`2` +
        dd$C3*dd$`3` +
        beta1*dd$dummy1 +
        beta2*dd$dummy2
    dd <- dd[order(dd$group,dd$rep,dd$time),]
    dd$cluster <- ifelse(dd$group %in% 1:B1, 1,2)
    ## COVARIANCE MATRICES FOR BOTH CLUSTES ----
    if(is.null(nu1)) nu1 <- b1
    if(is.null(nu2)) nu2 <- b2
    mkt1 <- mkt[1:(B1*C),]
    mkt2 <- mkt[-c(1:(B1*C)),]
    covMt1 <- covMatrix(market=mkt1,
                       group.name=group,
                       type.name=type,
                       mkt.name=num,
                       timeVec=mc$Time2,
                       sigPar=sigPar[,1],
                       tauPar=tauPar[,1],
                       corPar=corPar[,1],
                       funcMtx=nu1,
                       covType='Heterog',
                       truncateDec=8)
    covMt2 <- covMatrix(market=mkt2,
                       group.name=group,
                       type.name=type,
                       mkt.name=num,
                       timeVec=mc$Time2,
                       sigPar=sigPar[,2],
                       tauPar=tauPar[,2],
                       corPar=corPar[,2],
                       funcMtx=nu2,
                       covType='Heterog',
                       truncateDec=8)
    ## CREATE NOISY DATA
    dd$flag <- rep(1:(J*nRep), each=48)
    ddList <- split(x=dd,f=dd$group)
    obs <- lapply(ddList,
           function(x){
               jflag <- as.character(x$group[1])
               bflag <- x$cluster[1]
               if(bflag==1)
                   sigMtx <- covMt1[[jflag]]
               else
                   sigMtx <- covMt2[[jflag]]
               tmp <- tapply(x$signal,x$flag,
                             function(m){
                                 MASS::mvrnorm(n=1,
                                               mu=m,
                                               Sigma=sigMtx)
                             })
               unlist(tmp)
           })
    obs <- unlist(obs)
    attr(obs,'names') <- NULL
    dd$obs <- obs
    ## ORGANIZE FOR OUTPUT
    dd$flag <- NULL
    colnames(dd) <- c('group', 'time',
                      'c1','c2','c3',
                      'mc1','mc2','mc3',
                      'rep',
                      'temperature',
                      'dummy1',
                      'dummy2',
                      'signal',
                      'cluster',
                      'obs') ## the dependent variable
    attr(dd,"sigPar") = sigPar
    attr(dd,"corPar") = corPar
    attr(dd,"tauPar") = tauPar
    attr(dd,"tempPar") = tempPar
    attr(dd,"seed") = ifelse(is.null(seed), "Not specified", seed )
    attr(dd,"market") = mkt
    dd
}



