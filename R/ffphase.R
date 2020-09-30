#' @title Cube method with reduction of the auxiliary variables matrix
#'
#' @description Modified cube method.
#' This function reduces considerably the execution time when the matrix of auxiliary variables \code{X} contains lot of 0s.
#' It is based on the function \code{\link[sampling:samplecube]{samplecube}} from the package \code{sampling}.
#'
#'
#' @param X a matrix of size (N x p) of auxiliary variables on which the sample must be balanced.
#' @param pik a vector of size N of inclusion probabilities.
#'
#' @details In case where the number of auxiliary variables is great (i.e. p very large), even if we use the fast implementation proposed by
#' (Chauvet and Tille 2005), the problem is time consuming.
#' This function reduces considerably the execution time when the matrix of auxiliary variables \code{X} contains lot of 0s.
#' It considers a reduced matrix \code{X} by removing columns and rows that sum to 0 in the flight phase of the method (see  \code{\link{ReducedMatrix}} and \code{\link{ReducedFlightphase}}).
#'
#'
#' @return the updated vector of \code{pik} that contains only 0s and 1s that indicates if a unit is selected or not at each wave.
#'
#'
#' @author Esther Eustache \email{esther.eustache@@unine.ch}
#'
#'
#' @references
#' Chauvet, G. and Tille, Y. (2006). A fast algorithm of balanced sampling. Computational Statistics, 21/1:53-62
#'
#'
#' @seealso \code{\link[sampling:samplecube]{samplecube}}, \code{\link[sampling:landingcube]{landingcube}}, \code{\link{ReducedFlightphase}}, \code{\link{ReducedMatrix}}
#'
#'
#' @export
#' @examples
#' \dontrun{
#' 
#' N <- 100
#' n <- 10
#' p <- 4
#' 
#' pik <- sampling::inclusionprobabilities(runif(N),n)
#' pik <- rep(n/N,N)
#' 
#' X <- cbind(pik,matrix(rnorm(N*p),ncol= p))
#' s <- ffphase(X,pik) 
#' 
#' 
#' sum(s)
#' 
#' }
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
    
    B <- A[i_tmp,]

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
  }
  
  
  
  return(pik)
}








