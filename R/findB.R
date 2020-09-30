#' find the sub-matrix B in the flight phase
#'
#' @description
#' 
#'
#' @param X matrix of auxiliary variables.
#' @param Xcat matrix of categorical variables.
#'
#' @details
#'
#' The function find the smallest matrix B such that is contains only one additional row.
#'  It consecutively add the right number of row depending of the number of categroies that is added.
#'
#' @return matrix B
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
#'
#' N <- 100
#' Xcat <-data.frame(cat1 = rep(1:50,each = 2))
#' Xcat <- as.matrix(Xcat[-2,])
#' Xcat <- as.matrix(Xcat[-4,])
#' X <- matrix(rnorm(N),ncol = 1)
#'
#' findB(X,Xcat)
#'
#'
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
