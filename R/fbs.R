#' @title  Fast balanced Sampling
#' 
#' @description 
#' 
#' This function implements the method proposed by Hasler and Tille (2019). It should be used for selecting a sample for highly stratified population.
#'
#' @param X matrix of auxiliary variables.
#' @param strata matrix of categorical variables.
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
#' strata <- as.matrix(rep(1:n,each = N/n))
#' strata <- rep(1:n,each = N/n)
#'
#' pik <- sampling::inclusionprobastrata(strata,rep(1,n))
#' X <- as.matrix(cbind(matrix(c(x1,x2),ncol = 2)))
#' 
#' s <- fbs(X,strata,pik)
#' 
#'
#' @references 
#' Hasler, C. and Tille Y. (2014). Fast balanced sampling for highly stratified population.
#' \emph{Computational Statistics and Data Analysis}, 74, 81-94
#' 
#' @export
fbs <- function(X,strata,pik){
  
  H <- as.numeric(ncat(as.matrix(strata)))
  pik_tmp <- pik
  EPS <- 1e-8
  
  ##----------------------------------------------------------------
  ##                  Flightphase on each strata                   -
  ##----------------------------------------------------------------
  
  
  for(k in 1:H){
    # pik_tmp[strata == k] <- sampling::fastflightcube(cbind(pik[which(strata == k)],as.matrix(X[which(strata == k),])),
    #                                                pik[strata == k],
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
    Xnn <- disj(strata)
    for(k in 1:H){
      # print(k)
      
      i <- which(strata <= k & (pik_tmp > EPS & pik_tmp < (1-EPS)))
      
      # strata_tmp2 <- Xnn[i,1:k]
      # strata_tmp2 <- disjMatrix(as.matrix(strata[i]))
      
      strata_tmp2 <- disj(strata[i])
      strata_tmp2 <- strata_tmp2*pik_tmp[i]
      
      X_tmp <- as.matrix((X[i,]*pik_tmp[i]/pik[i]))
      
      # pik_tmp[i] <- sampling::fastflightcube(as.matrix(cbind(X_tmp,strata_tmp2)),
      #                                        pik_tmp[i],
      #                                        1,
      #                                        comment = FALSE)
      pik_tmp[i] <- ffphase(as.matrix(cbind(X_tmp,strata_tmp2)),
                                             pik_tmp[i])
    }
  }
  sum(pik_tmp)

  ##---------------------------------------------------------------
  ##              Landing by suppression of variables             -
  ##---------------------------------------------------------------
  
  i <- which(pik_tmp > EPS & pik_tmp < (1-EPS))

  
  if(length(i) != 0){
  strata_tmp3 <- as.matrix(Xnn)
  # strata_tmp3 <- as.matrix(Xnn[i,])
  pik_tmp <- landingRM(as.matrix(cbind(strata_tmp3, X/pik)),
                          pik_tmp)
  # pik_tmp[i] <- landingRM(as.matrix(cbind(Xnn[i,], X[i,])),
  #                         pik_tmp[i],
  #                         pik[i])
  
  }
  return(round(pik_tmp,10))
}
