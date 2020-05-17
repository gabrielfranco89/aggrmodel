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
    fit_init <- lm(data$y~X-1)
    beta_init <- coef(fit_init)
    sigma_fit <- sqrt(summary(fit_init)$sigma/(I))

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
             betaCov_init <- rep(sigma_fit, C*n_basis_cov)
        }
        else{
            if(length(betaCov_init)!=C*n_basis_cov) stop("betaCov_init must have the same length as number of basis for covariance")
        } # end if/else
        if(is.null(corPar_init))
            corPar <- rep(1, C)
        else
            corPar <- corPar_init
        parIn <- c(betaCov_init,
                   corPar)
        lowBoundVec <- c(rep(1e-8,times=(C*n_basis_cov)), ##beta
                         rep(1e-4,C)) ## corPar
        ubCor <- ifelse(is.null(truncateDec), Inf, log(10^truncateDec))
        upperBoundVec <- c(rep(Inf,times=(C*n_basis_cov)),
                           rep(ubCor,C))
     }# end if heterog

    output <- list()
    output$X <- X
    output$beta <- beta_init
    output$lb <- lowBoundVec
    output$ub <- upperBoundVec
    output$par <- parIn
    output$n_basis_cov <- n_basis_cov
    output
}
