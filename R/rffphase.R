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
#' @examples
#' set.seed(1)
#' ## Matrix of 8 auxilary variables and 10 units with lot of 0s ##
#'
#' ## Inclusion probabilities with 10 units ##
#' pik <- sampling::inclusionprobabilities(runif(40),10)
#' X   <- matrix(sample(c(0,0,0,1),240,replace=TRUE), nrow = 40, ncol =  8)
#' X <- cbind(pik,X)
#' ## Cube method ##
#' s <- ReducedSamplecube(X, pik)
#' # s <- sampling::samplecube(X,pik)
#' s
#' sum(pik)
#' sum(s)
#'
#' @export
rffphase <- function(X, pik, redux = TRUE){

  ##----------------------------------------------------------------
  ##                        Initialization                         -
  ##----------------------------------------------------------------


  EPS = 1e-8
  A = X/pik
  J = ncol(X)


  ##----------------------------------------------------------------
  ##                Number of non 0-1 inclusion prob               -
  ##----------------------------------------------------------------




  i <- which(pik > EPS & pik < (1-EPS))
  i_size = length(i)
  if(i_size >= (J+1)){
    i <- i[1:(J+1)]
    B = A[i,]
  }else{
    B = A[i,]
  }

  # i <- which(pik > EPS & pik < (1-EPS))[1:(J+1)]
  # i_size = length(i)



  ##---------------------------------------------------------------
  ##                          Main loop                           -
  ##---------------------------------------------------------------

  while(i_size > 0){


    ##-----------------------
    ##  if redux is desired
    ##-----------------------

    if(redux == TRUE){
      # pik_tmp <- pik[i]
      # tmp <- ReducedMatrix(B)
      # B_tmp <- tmp$B
      # pik_tmp[tmp$ind_row]<- onestep(B_tmp, pik_tmp[tmp$ind_row],EPS)
      # pik[i] <- pik_tmp

      pik_tmp <- pik[i]
      tmp <- ReducedMatrix(B)
      B_tmp <- tmp$B

      print(dim(B_tmp))
      check <- onestep(B_tmp, pik_tmp[tmp$ind_row],EPS)

      if(is.null(check)){
        break;
      }else{
        pik_tmp[tmp$ind_row]<- check
        pik[i] <- pik_tmp
      }


    }else{

      check <- onestep(B,pik[i],EPS)
      if(is.null(check)){
        break;
      }else{
        pik[i] <- check
      }

    }


    ##------------
    ##  update i
    ##------------

    i <- which(pik > EPS & pik < (1-EPS))
    i_size = length(i)



    ##-----------------------------------
    ##  Depending if we have enough row
    ##-----------------------------------

    if(i_size >= (J+1)){
      i <- i[1:(J+1)]
    }
    B <-A[i,]

  }


  return(pik)
}


