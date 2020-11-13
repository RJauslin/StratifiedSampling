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
#' X <- matrix(rgamma(N*p,4,25),ncol = p)
#' strata <- as.matrix(rep(1:n,each = N/n))
#' pik <- rep(n/N,N)
#' 
#' balstrat(X,strata,pik)
#' @export
balstrat <- function (X, strata, pik) 
{
  
  H <- as.numeric(ncat(strata))
  pik_tmp <- pik
  EPS <- 1e-7
  
  ##----------------------------------------------------------------
  ##                  Flightphase on each strata                   -
  ##----------------------------------------------------------------
  
  
  for(k in 1:H){
    # pik_tmp[Xcat == k] <- sampling::fastflightcube(cbind(pik[which(Xcat == k)],as.matrix(X[which(Xcat == k),])),
    #                                                pik[Xcat == k],
    #                                                comment = FALSE)
    pik_tmp[strata == k] <- ffphase(as.matrix(cbind(pik[which(strata == k)],as.matrix(X[which(strata == k),]))),
                                    pik[strata == k])
    
  }
  sum(pik_tmp)
  
  ##----------------------------------------------------------------
  ##          Flightphase on the uninon of strata U1 -- Uk         -
  ##----------------------------------------------------------------
  i <- which(pik_tmp > EPS & pik_tmp < (1-EPS))
  if(length(i) != 0){
    
    Xcat_tmp2 <- disjMatrix(as.matrix(strata[i,]))
    Xcat_tmp2 <- Xcat_tmp2*pik_tmp[i]
    
    X_tmp <- as.matrix((X[i,]*pik_tmp[i]/pik[i]))
    pik_tmp[i] <- ffphase(as.matrix(cbind(X_tmp,Xcat_tmp2)),pik_tmp[i])
    
    
    # for(k in 1:H){
    #   # print(k)
    #   
    #   i <- which(Xcat <= k & (pik_tmp > EPS & pik_tmp < (1-EPS)))
    #   
    #   # Xcat_tmp2 <- Xnn[i,1:k]
    #   Xcat_tmp2 <- disjMatrix(as.matrix(Xcat[i,]))
    #   Xcat_tmp2 <- Xcat_tmp2*pik_tmp[i]
    #   
    #   X_tmp <- as.matrix((X[i,]*pik_tmp[i]/pik[i]))
    #   
    #   # pik_tmp[i] <- sampling::fastflightcube(as.matrix(cbind(X_tmp,Xcat_tmp2)),
    #   #                                        pik_tmp[i],
    #   #                                        1,
    #   #                                        comment = FALSE)
    #   pik_tmp[i] <- ffphase(as.matrix(cbind(X_tmp,Xcat_tmp2)),
    #                         pik_tmp[i])
    # }
  }
  sum(pik_tmp)
  
  ##---------------------------------------------------------------
  ##              Landing by suppression of variables             -
  ##---------------------------------------------------------------
  
  i <- which(pik_tmp > EPS & pik_tmp < (1-EPS))
  
  
  if(length(i) != 0){
    Xnn <- disjMatrix(as.matrix(strata))
    Xcat_tmp3 <- as.matrix(Xnn)
    # Xcat_tmp3 <- as.matrix(Xnn[i,])
    pik_tmp <- landingRM(as.matrix(cbind(Xcat_tmp3, X/pik)),
                         pik_tmp)
    # pik_tmp[i] <- landingRM(as.matrix(cbind(Xnn[i,], X[i,])),
    #                         pik_tmp[i],
    #                         pik[i])
    
  }
  return(round(pik_tmp,10))
  
}