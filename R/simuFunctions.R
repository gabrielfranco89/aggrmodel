#' Simulate aggregated data
#' @name simuData
#' @param nGroups Number of groups to be simulated
#' @param market Optional. Insert the market as a data frame with the columns 'group', 'type' and 'num'
#' @param market_range If market not inserted, the range to be sampled
#' @param sigPar Parameters of sigma: vector of length 3
#' @param tauPar Parameters of tau: vector of length 3
#' @param corPar Parameters of correlation: vector of length 3
#' @param seed Optional. Seed to be used for simulation
#' @param tempPar Default:1. If you do not want temperature effect, set temPar=0
#'
#' @return
#' Return a data frame with 20 replicates and 3 types of consumer observed in 48 times.
#'
#' @examples
#' dd = simuData()
#' dd = simuData(nGroups=4)
simuData <- function(nGroups=10,
                     market=NULL,
                     market_range = 1:100,
                     sigPar = c(2,2,2),
                     tauPar = c(.5,.5,.5),
                     corPar = c(13,13,13),
                     seed=NULL,
                     tempPar=1
                     ){
    ## Preamble
    if(nGroups %% 2 != 0)
        stop("Please use an even number of groups to make it possible to divide into two clusters. \nIf you do not want to use clusters, simply double the number of groups and use just the half of generated data :)")
    J <- nGroups
    C <- 3
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
    b <- rbind(
        do.call(rbind,replicate(n=J/2,expr=b1, simplify=FALSE)),
        do.call(rbind,replicate(n=J/2,expr=b2, simplify=FALSE))
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
    tmp <- do.call(rbind, replicate(n=20,
                                   expr=dd,
                                   simplify=FALSE))
    tmp$rep <- rep(1:20, each=nr)
    temp$Time <- temp$Time/24
    colnames(temp) <- c('time','rep','temp')
    dd <- merge(tmp,temp)
    rm(tmp)
    dd <- dd[order(dd$group,dd$rep,dd$time),]
    dd$cluster <- ifelse(dd$group %in% 1:(J/2), 1,2)
    dd$signal <- dd$aggr*( log(log(dd$temp)) )^tempPar
    ## COVARIANCE MATRICES FOR BOTH CLUSTES
    mkt1 <- mkt[1:(nrow(mkt)/2),]
    mkt2 <- mkt[-c(1:(nrow(mkt)/2)),]
    covMt1 <- covMatrix(market=mkt1,
                       group.name=group,
                       type.name=type,
                       mkt.name=num,
                       timeVec=mc$Time2,
                       sigPar=sigPar,
                       tauPar=tauPar,
                       corPar=corPar,
                       funcMtx=b1,
                       covType='Heterog',
                       truncateDec=8)
    covMt2 <- covMatrix(market=mkt2,
                       group.name=group,
                       type.name=type,
                       mkt.name=num,
                       timeVec=mc$Time2,
                       sigPar=sigPar,
                       tauPar=tauPar,
                       corPar=corPar,
                       funcMtx=b2,
                       covType='Heterog',
                       truncateDec=8)
    ## CREATE NOISY DATA
    dd$flag <- rep(1:(J*20), each=48)
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
    dd
}



