#' @title Fast flight phase
#'
#' @description 
#' 
#' This function is computing the flight phase of the cube method. 
#'
#' @param X a matrix of size (N x p) of auxiliary variables on which the sample must be balanced.
#' @param pik a vector of size N of inclusion probabilities.
#'
#' @details 
#' This function implements the method proposed by (Chauvet and Tille 2006). It progressively transform the vector of inclusion probabilities \code{pik} into a sample while respecting the balancing equations.
#' The algorithm stops when the null space of the matrix \eqn{B} is empty. For more information see (Chauvet and Tille 2006).
#' 
#' The function uses the function \code{\link[MASS:Null]{Null}} to find the null space of the matrix \eqn{B}.
#'
#' @return the updated vector of \code{pik} that contains 0 and 1 for unit that are rejected or selected.
#'
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
#' #pik <- rep(n/N,N)
#' pik <- inclusionprobabilities(runif(N),n)
#' X <- cbind(pik,matrix(rnorm(N*p),ncol= p))
#' 
#' pikstar <- ffphase(X,pik) 
#' 
#' t(X/pik)%*%pikstar
#' t(X/pik)%*%pik
#' pikstar
#' 
#' 
#' 
#' 
#' X <- cbind(matrix(rep(0,2*10),nrow = 2),c(0.45,0.55))
#' X <- matrix(c(1,1),ncol = 1)
#' pik <- c(0.1260076, 0.8739924)
#' 
#' fastflightcube(X,pik)
#' s <- ffphase(X*pik,pik)
#' 
#' t(X/pik)%*%pik
#' 
#' @export
ffphase <- function(X, pik){

  ##----------------------------------------------------------------
  ##                        Initialization                         -
  ##----------------------------------------------------------------
  EPS = 1e-8
  pikInit <- pik
  # AInit <- X/pik
  A <- X/pikInit
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
  
  # print(i_size)
  # print(p)
  while(i_size > 0){
  
    # print(i_size)
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
    # print(tmp)
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
    # A <- X*(pik/pikInit)
  }
  
  return(pik)
}








