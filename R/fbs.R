

#' Title
#'
#' @param X
#' @param Xcat
#' @param pik
#'
#' @return
#' @export
#'
#' @examples
#' rm(list = ls())
#' N <- 10000
#' n <- 1000
#' x1 <- rgamma(N,4,25)
#' x2 <- rgamma(N,4,25)
#'
#' Xcat <- as.matrix(rep(1:n,each = N/n))
#'
#' pik <- inclusionprobastrata(Xcat,rep(1,n))
#' X <- as.matrix(cbind(matrix(c(x1,x2),ncol = 2)))
#'
#' A <- cbind(X,disjunctive(as.matrix(Xcat)))
#'
#' system.time(s <- fbs(X,Xcat,pik))
#' system.time(s <- stratifiedcube(cbind(pik,X),Xcat,pik))
#'
fbs <- function(X,Xcat,pik){
  
  H <- as.numeric(ncat(Xcat))
  pik_tmp <- pik
  EPS <- 1e-7
  
  ##----------------------------------------------------------------
  ##                  Flightphase on each strata                   -
  ##----------------------------------------------------------------
  
  
  for(k in 1:H){
    # pik_tmp[Xcat == k] <- sampling::fastflightcube(cbind(pik[which(Xcat == k)],as.matrix(X[which(Xcat == k),])),
    #                                                pik[Xcat == k],
    #                                                comment = FALSE)
    pik_tmp[Xcat == k] <- ffphase(as.matrix(cbind(pik[which(Xcat == k)],as.matrix(X[which(Xcat == k),]))),
                                                   pik[Xcat == k])
    
  }
  
  
  ##----------------------------------------------------------------
  ##          Flightphase on the uninon of strata U1 -- Uk         -
  ##----------------------------------------------------------------
  Xnn <- disjMatrix(as.matrix(Xcat))
  for(k in 1:H){
    # print(k)
    
    i <- which(Xcat <= k & (pik_tmp > EPS & pik_tmp < (1-EPS)))
    
    # Xcat_tmp2 <- Xnn[i,1:k]
    Xcat_tmp2 <- disjMatrix(as.matrix(Xcat[i,]))
    Xcat_tmp2 <- Xcat_tmp2*pik_tmp[i]
    
    X_tmp <- as.matrix((X[i,]*pik_tmp[i]/pik[i]))
    
    # pik_tmp[i] <- sampling::fastflightcube(as.matrix(cbind(X_tmp,Xcat_tmp2)),
    #                                        pik_tmp[i],
    #                                        1,
    #                                        comment = FALSE)
    pik_tmp[i] <- ffphase(as.matrix(cbind(X_tmp,Xcat_tmp2)),
                                           pik_tmp[i])
  }
  
  
  # sum(pik_tmp)
  
  ##---------------------------------------------------------------
  ##              Landing by suppression of variables             -
  ##---------------------------------------------------------------
  
  i <- which(pik_tmp > EPS & pik_tmp < (1-EPS))
  
  # pik_tmp[i] <- samplecube(as.matrix(cbind(X[i,],Xnn[i,])),
  #                          pik_tmp[i],
  #                          order = 1,
  #                          comment = FALSE,
  #                          method = 2)
  
  pik_tmp[i] <- landingRM(as.matrix(cbind(X[i,],Xnn[i,])),
                          pik_tmp[i],
                          pik[i])
  
  
  return(round(pik_tmp,10))
}
