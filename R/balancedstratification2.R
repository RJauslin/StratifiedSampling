#' @title Balanced Stratification
#' 
#' @description 
#' 
#' Select a strafitifed balanced sample. The function is similar to \code{\link[sampling:balancedstratification]{balancedstratification}} unless it uses internal functions that are more stable.
#' 
#' @param X matrix of auxiliary variables.
#' @param strata vector of integer that specifies the stratification.
#' @param pik vector of inclusion probabilities.
#'
#' @return A vector with elements equal to 0 or 1. The value 1 indicates that the unit is selected while the value 0 is for rejected units.
#'
#' @examples
#' 
#' N <- 1000
#' n <- 100
#' p <- 4
#' X_tmp <- matrix(rgamma(N*p,4,25),ncol = p)
#' strata <- as.matrix(rep(1:n,each = N/n))
#' pik <- rep(n/N,N)
#' X <- cbind(pik,X_tmp)
#' balstrat(X,strata,pik)
#' @export
balstrat <- function (X, strata, pik) 
{
  
  
  ##----------------------------------------------------------------
  ##                          Initializing                         -
  ##----------------------------------------------------------------
  
  
  EPS = 1e-8
  strata = sampling::cleanstrata(strata)
  H = max(strata)
  N = dim(X)[1]
  pikstar = rep(0, times = N)
  
  
  ##----------------------------------------------------------------
  ##                            Step 1                             -
  ##----------------------------------------------------------------
  
  
  for (h in 1:H) {
    pikstar[strata == h] = ffphase(cbind(pik[strata == h],X[strata == h, ]), pik[strata == h])
  }
  
  
  ##----------------------------------------------------------------
  ##                            Step 2                             -
  ##----------------------------------------------------------------
  
  
  XN = cbind(sampling::disjunctive(strata) * pik, X)/pik * pikstar
  if (is.null(colnames(X)) == FALSE) 
    colnames(XN) <- c(paste("Stratum", 1:H, sep = ""), 
                      colnames(X))
  
  i <- which(pikstar > EPS & pikstar < (1-EPS))
  pikstar[i] = ffphase(as.matrix(cbind(X[i,],XN[i,])), pikstar[i])
  
  
  
  ##----------------------------------------------------------------
  ##                            Step 3                             -
  ##----------------------------------------------------------------
  
  
  i <- which(pikstar > EPS & pikstar < (1-EPS))
  pikstar[i] <- landingRM(as.matrix(cbind(X[i,],XN[i,])),
                          pikstar[i],
                          pik[i])
  pikstar <- round(pikstar,8)

  return(pikstar)
}