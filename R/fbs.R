#' @title  Fast balanced Sampling
#' 
#' @description 
#' 
#' This function implements the method proposed by Hasler and Tille (2019). It should be used for selecting a sample for highly stratified population.
#'
#' @param X matrix of auxiliary variables.
#' @param Xcat matrix of categorical variables.
#' @param pik vector of inclusion probabilities.
#'
#' @details 
#' If the number of element selected in each strata is not equal to a integer, the function can be very time-consuming.
#'
#' @return A vector with elements equal to 0 or 1. The value 1 indicates that the unit is selected while the value 0 is for rejected units.
#'
#' @examples
#' 
#' N <- 1000
#' n <- 100
#' x1 <- rgamma(N,4,25)
#' x2 <- rgamma(N,4,25)
#'
#' Xcat <- as.matrix(rep(1:n,each = N/n))
#'
#' pik <- sampling::inclusionprobastrata(Xcat,rep(1,n))
#' X <- as.matrix(cbind(matrix(c(x1,x2),ncol = 2)))
#' A <- cbind(X,sampling::disjunctive(as.matrix(Xcat)))
#' 
#' s <- fbs(X,Xcat,pik)
#' 
#' 
#' 
#' 
#' ##----------------------------------------------------------------
#' ##                              Data                             -
#' ##----------------------------------------------------------------
#' 
#' rm(list = ls())
#' N <- 1000
#' n <- 100
#' x1 <- rgamma(N,4,25)
#' x2 <- rgamma(N,4,25)
#' strata <- as.matrix(rep(1:n,each = N/n))
#' 
#' 
#' 
#' ##---------------------------------------------------------------
#' ##                        different cases                       -
#' ##---------------------------------------------------------------
#' 
#' 
#'  #-------- CASE 0 pik equal and only pik as variable
#'  
#'  pik <- inclusionprobastrata(strata,rep(1,n))
#'  # X <- as.matrix(pik)
#'  
#'  X <- matrix(c(x1,x2),ncol = 2)
#'  Xcat <- strata
#'  system.time(s <- fbs(X,strata,pik))
#'  sum(s)
#'  t(X/pik)%*%s
#'  t(X/pik)%*%pik
#'  
#'  t(disjMatrix(strata))%*%s
#'  t(disjMatrix(strata))%*%pik
#' 
#' 
#'
#' @references 
#' Hasler, C. and Tille Y. (2014). Fast balanced sampling for highly stratified population.
#' \emph{Computational Statistics and Data Analysis}, 74, 81-94
#' 
#' @export
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
  sum(pik_tmp)
  
  ##----------------------------------------------------------------
  ##          Flightphase on the uninon of strata U1 -- Uk         -
  ##----------------------------------------------------------------
  i <- which(pik_tmp > EPS & pik_tmp < (1-EPS))
  if(length(i) != 0){
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
  }
  sum(pik_tmp)

  ##---------------------------------------------------------------
  ##              Landing by suppression of variables             -
  ##---------------------------------------------------------------
  
  i <- which(pik_tmp > EPS & pik_tmp < (1-EPS))

  
  if(length(i) != 0){
  Xcat_tmp3 <- as.matrix(Xnn[i,]*pik_tmp[i])
  pik_tmp[i] <- landingRM(as.matrix(cbind(Xcat_tmp3, X[i,]/pik[i]*pik_tmp[i])),
                          pik_tmp[i])
  # pik_tmp[i] <- landingRM(as.matrix(cbind(Xnn[i,], X[i,])),
  #                         pik_tmp[i],
  #                         pik[i])
  
  }
  return(round(pik_tmp,10))
}
