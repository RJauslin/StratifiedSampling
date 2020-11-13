#' @title Fast flight phase of the cube method
#'
#' @description 
#' 
#' This function is computing the flight phase of the cube method. (Chauvet and Tille 2006)
#'
#' @param X a matrix of size (\eqn{N} x \eqn{p}) of auxiliary variables on which the sample must be balanced.
#' @param pik vector of inclusion probabilities.
#'
#' @details 
#' This function implements the method proposed by (Chauvet and Tille 2006). It recursively transform the vector of inclusion probabilities \code{pik} into a
#' sample that respect the balancing equations. The algorithm stops when the null space of the sub-matrix \eqn{B} is empty.
#' For more information see (Chauvet and Tille 2006).
#' 
#' The function uses the function \code{\link[MASS:Null]{Null}} to find the null space of the sub-matrix \eqn{B}.
#'
#' @return Updated vector of \code{pik} that contains 0 and 1 for unit that are rejected or selected.
#'
#' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
#'
#'
#' @references
#' Chauvet, G. and Tille, Y. (2006). A fast algorithm of balanced sampling. Computational Statistics, 21/1:53-62
#'
#'
#' @seealso \code{\link[sampling:samplecube]{fastflightphase}}, \code{\link[BalancedSampling:flightphase]{flightphase}}. 
#'
#' @examples
#' 
#' N <- 100
#' n <- 10
#' p <- 4
#' 
#' pik <- rep(n/N,N)
#' X <- cbind(pik,matrix(rgamma(N*p,4,25),ncol= p))
#' 
#' pikstar <- ffphase(X,pik) 
#' t(X/pik)%*%pikstar
#' t(X/pik)%*%pik
#' pikstar
#' 
#' 
#' @export
ffphase <- function(X, pik){

  ##----------------------------------------------------------------
  ##                        Initialization                         -
  ##----------------------------------------------------------------
  EPS = 1e-8
  A <- X/pik
  N <- length(pik)  
  p = ncol(X)
  
  ##----------------------------------------------------------------
  ##                Number of non 0-1 inclusion prob               -
  ##----------------------------------------------------------------
  
  i <- which(pik > EPS & pik < (1-EPS))
  i_size = length(i)
  
  ##----------------------------------------------------------------
  ##                          flight phase                         -
  ##----------------------------------------------------------------
  
  
  while(i_size > 0){
    
    ##------ Find B
    if(i_size >= (p+1)){
      i_tmp <- i[1:(p+1)]
      pik_tmp <- pik[i_tmp]
    }else{
      i_tmp <- i
      pik_tmp <- pik[i_tmp]
    }
    
    B <- as.matrix(A[i_tmp,])

    ##------ onestep and check if null
    tmp <-  onestep(B,pik_tmp,EPS)
    if(is.null(tmp)){
      break;
    }else{
      pik[i_tmp] <- tmp
    }
    
    ##------ update i
    i <- which(pik > EPS & pik < (1-EPS))
    i_size = length(i)
    
    if(i_size == 1){
      break;
    }
  
  }
  
  return(pik)
}








