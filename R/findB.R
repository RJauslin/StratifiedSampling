#' @title Find best sub-matrix B in stratifiedcube
#'
#' @description
#' This function is computing a sub-matrix used in \code{\link{stratifiedcube}}.
#' 
#' @param X A matrix of size (\eqn{N} x \eqn{p}) of auxiliary variables on which the sample must be balanced.
#' @param strata A vector of integer that specifies the stratification.
#'
#' @details
#'
#' The function find the smallest matrix B such that is contains only one more row than the number of column.
#' It consecutively add the right number of row depending of the number of categories that is added.
#'
#' @return A list of two components. The sub-matrix of \code{X} and the corresponding disjunctive matrix.
#'  If we use the function \code{cbind} to combine the two matrix, the resulting matrix have only one more row than the number of column. 
#'
#' @export
#'
#' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
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
