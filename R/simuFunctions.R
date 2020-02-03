#' Simulate aggregated data
#' @name createSimuData
#' @param B1 Number of groups from cluster 1
#' @param B2 Number of groups from cluster 2
#' @param nRep Number of replicates
#' @param market Optional. Insert the market as a data frame with the columns 'group', 'type' and 'num'
#' @param market_range If market not inserted, the range to be sampled
#' @param sigPar Parameters of sigma: vector of length 3
#' @param tauPar Parameters of tau: vector of length 3
#' @param corPar Parameters of correlation: vector of length 3
#' @param seed Optional. Seed to be used for simulation
#' @param tempPar Default:1. If you do not want temperature effect, set temPar=0
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
                           market_range = 20:80,
                           sigPar = c(2,2,2,4,4,4),
                           tauPar = c(.5,.5,.5, .4,.4,.4),
                           corPar = c(4,4,4,6,6,6),
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
    temp <- simuTemperature
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
    dd$aggr <- dd$C1*dd$`1` +
        dd$C2*dd$`2` +
        dd$C3*dd$`3`
    ## ADD TEMPERATURES
    nr <- nrow(dd)
    tmp <- do.call(rbind, replicate(n=nRep,
                                   expr=dd,
                                   simplify=FALSE))
    tmp$rep <- rep(1:nRep, each=nr)
    temp$Time <- temp$Time/24
    colnames(temp) <- c('time','rep','temp')
    dd <- merge(tmp,temp)
    rm(tmp)
    dd <- dd[order(dd$group,dd$rep,dd$time),]
    dd$cluster <- ifelse(dd$group %in% 1:B1, 1,2)
    dd$signal <- dd$aggr*( log(log(dd$temp)) )^tempPar
    ## COVARIANCE MATRICES FOR BOTH CLUSTES
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
                       funcMtx=b1,
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
                       funcMtx=b2,
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
    colnames(dd) <- c('time', 'rep','group',
                      'c1','c2','c3',
                      'mc1','mc2','mc3',
                      'real_aggr',
                      'temperature',
                      'cluster',
                      'signal_unnoisy',
                      'obs') ## the dependent variable
    attr(dd,"sigPar") = sigPar
    attr(dd,"corPar") = corPar
    attr(dd,"tauPar") = tauPar
    attr(dd,"tempPar") = tempPar
    attr(dd,"seed") = ifelse(is.null(seed), "Not specified", seed )
    attr(dd,"market") = mkt
    dd
}



