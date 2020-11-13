#' @title Find best sub-matrix B in stratifiedcube
#'
#' @description
#' This function is computing a sub-matrix used in stratifiedcube.
#' 
#' @param X A matrix of auxiliary variables.
#' @param strata A vector of categories.
#'
#' @details
#'
#' The function find the smallest matrix B such that is contains only one more row than the number of column.
#' It consecutively add the right number of row depending of the number of categories that is added.
#'
#' @return A list of two components. The sub-matrix of X and the corresponding disjunctive matrix.
#'  If we \code{cbind} the two matrix, the resulting matrix have only one more row than the number of column. 
#'
#' @export
#'
#' @examples
#' N <- 1000
#' strata <-  sample(x = 1:6, size = N, replace = TRUE)
#'
#' p <- 3
#' X <- matrix(rnorm(N*p),ncol = 3)
#' findB(X,strata)
#'
findB <- function(X,
                  strata){

  strata <- as.matrix(strata)
  X <- as.matrix(X)
  eps <- 1e-9
  N <- nrow(X)
  pInit <- ncol(X)

  if(pInit > nrow(strata)){
    strata_tmp <- as.matrix(strata[1:(pInit+1),])
  }else{
    strata_tmp <- as.matrix(strata[1:pInit,])
  }
  
  
  nstrata <- sum(ncat(strata_tmp))
  nstrata_tmp <- 0

  while(nstrata != nstrata_tmp){
    nstrata_tmp = nstrata
    p =  pInit  + nstrata
    if(p >= nrow(X)){
      p <- nrow(X)-1
      strata_tmp <- as.matrix(strata[1:(p+1),])
      nstrata <- sum(ncat(strata_tmp))
      break;
    }
    strata_tmp <- as.matrix(strata[1:(p+1),])
    nstrata <- sum(ncat(strata_tmp))
  }

  strata_tmp <- as.matrix(strata_tmp)
  disj_strata <- disjMatrix(strata_tmp)

  return(list(X = as.matrix(X[1:(p+1),]),Xcat = disj_strata))

}
