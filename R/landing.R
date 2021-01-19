#' @title Landing by suppression of variables
#'
#' @param X matrix of auxiliary variables on which the sample must be balanced. (The matrix should be divided by the original inclusion probabilities.)
#' @param pikstar vector of updated inclusion probabilities by the flight phase. See \code{\link{ffphase}}
#'
#' @return A vector with elements equal to 0 or 1. The value 1 indicates that the unit is selected while the value 0 is for rejected units.
#' 
#' @description
#' This function performs the landing phase of the cube method using suppression of variables proposed by Chauvet and Tillé (2006).
#'
#' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
#'
#' @seealso \code{\link{fbs}}, \code{\link{balstrat}}.
#' 
#' @references
#' Chauvet, G. and Tillé, Y. (2006). A fast algorithm of balanced sampling. \emph{Computational Statistics}, 21/1:53-62
#'
#' @export
#' @examples
#' 
#' N <- 1000
#' n <- 10
#' p <- 4
#' pik <- rep(n/N,N)
#' X <- cbind(pik,matrix(rgamma(N*p,4,25),ncol= p))
#' pikstar <- ffphase(X,pik) 
#' s <- landingRM(X/pik,pikstar)
#' sum(s)
#' t(X/pik)%*%pik
#' t(X/pik)%*%pikstar
#' t(X/pik)%*%s
landingRM <- function(X,pikstar){


  ##----------------------------------------------------------------
  ##                          Initializing                         -
  ##----------------------------------------------------------------


  EPS = 1e-8
  N = nrow(X)
  i = which(pikstar > EPS & pikstar < (1 - EPS))
  i_size = length(i)
  
  pikland <- pikstar[i]

  Xland <- X[i,]
  Nland = length(pikland)
  nland = sum(pikland)
  p <- ncol(Xland)
  
  
  j <-  which(pikland > EPS & pikland < (1 - EPS))
  j_size <- length(j)
  

  ##---------------------------------------------------------------
  ##                          Main loop                           -
  ##---------------------------------------------------------------
  
  
  for(k in 0:(p-1)){

    Bland <- Xland[j,]*pikland[j]
    Bland <- Bland[,1:(p-k)]

    kern <- MASS::Null(Bland)
    if(length(kern)!=0){
      
      pikland[j] <- ffphase(as.matrix(Bland),pikland[j])
      
      j = which(pikland > EPS & pikland < (1 - EPS))
      j_size <- length(j)
      # print(i_size)
    }
    if(j_size <= 1){
      break;
    }

  }
  
  pikstar[i] = pikland
  i <- which(pikstar > EPS & pikstar < (1 - EPS))
  
  if(length(i) != 0){
    pikstar[i] <- stats::rbinom(1,1,pikstar[i])
    # cat("it remains",length(i), "units that are not put to 0 or 1")
  }

  return(round(pikstar,10))
}
