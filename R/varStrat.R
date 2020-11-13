#' @title Approximated variance
#' 
#' @description 
#' 
#' @param X A matrix of size (\eqn{N} x \eqn{p}) of auxiliary variables on which the sample must be balanced.
#' @param strata A vector of integer that represents the categories.
#' @param pik vector of inclusion probabilities.
#' @param s a sample (vector of 0 and 1, if rejected or selected).
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
#' rm(list = ls())
#' N <- 1000
#' n <- 400
#' x1 <- rgamma(N,4,25)
#' x2 <- rgamma(N,4,25)
#' 
#' strata <- as.matrix(rep(1:40,each = 25)) # 25 strata
#' Xcat <- disjMatrix(strata)
#' 
#' #-------- CASE 1 pik equal and integer in each strata 
#' 
#' pik <- inclusionprobastrata(strata,rep(1,n))
#' X <- as.matrix(matrix(c(x1,x2),ncol = 2))
#'  
#' s <- stratifiedcube(cbind(pik,X),strata,pik)
#'  
#'  
#'  y <- 500 + 5*x1 + 5*x2 + rnorm(1000,0,270)
#'  
#'  y_ht <- sum(y[which(s==1)]/pik[which(s == 1)])
#'  
#'  sum(y)
#'  sum(y_ht)
#'  
#'  varApp(X,strata,pik,s,y)
#'  
varApp <- function(X,strata,pik,y){
  
  Xcat <- disjMatrix(strata)
  X_tmp <- cbind(Xcat,X)
  N <- length(pik)
  H <- ncol(Xcat)
  q <- ncol(X)
  
  beta <- matrix(rep(0,(H+q)^2),ncol = (H+q),nrow = (H+q))
  
  #compute beta
  for(k in 1:N){
    z_k <- X_tmp[k,]
    b_k <- pik[k]*(1-pik[k])*(N/(N-(H+q)))
    beta <- beta + b_k*(z_k/pik[k])%*%t(z_k/pik[k])
  } 
  beta <- solve(beta)
  
  #compute tmp
  tmp <- rep(0,H+q)
  for(k in 1:N){
    z_k <- X_tmp[k,]
    b_k <- pik[k]*(1-pik[k])*(N/(N-(H+q)))
    tmp <- tmp + b_k*(y[k]/pik[k])*(z_k/pik[k])
  } 
  
  beta <- beta%*%tmp
  
  #compute v
  v <- 0
  for(k in 1:N){
    z_k <- X_tmp[k,]
    b_k <- pik[k]*(1-pik[k])*(N/(N-(H+q)))
    
    v <- v + b_k*( (y[k]/pik[k]) - t(beta)%*%(z_k/pik[k]) )^2
  } 
  
  v
  return(v)
  
  
}



#' Title
#'
#' @param X 
#' @param strata 
#' @param pik 
#' @param s 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
#' rm(list = ls())
#' N <- 1000
#' n <- 400
#' x1 <- rgamma(N,4,25)
#' x2 <- rgamma(N,4,25)
#' 
#' strata <- as.matrix(rep(1:40,each = 25)) # 25 strata
#' Xcat <- disjMatrix(strata)
#' 
#' #-------- CASE 1 pik equal and integer in each strata 
#' 
#' pik <- rep(n/N,N)
#' X <- as.matrix(matrix(c(x1,x2),ncol = 2))
#'  
#' s <- stratifiedcube(cbind(pik,X),strata,pik)
#' s <- fbs(X,strata,pik)
#'  
#'  y <- 20*strata + rnorm(1000,120)
#'  y2 <- 500 + 5*x1 + 5*x2 + rnorm(1000,0,270)
#'  
#'  y_ht <- sum(y[which(s==1)]/pik[which(s == 1)])
#'  
#'  (sum(y_ht) - sum(y))^2
#'  varEst(X,strata,pik,s,y)
#'  varApp(X,strata,pik,s,y)
#'  
varEst <- function(X,strata,pik,s,y){
  
  Xcat <- disjMatrix(strata)
  X_tmp <- cbind(Xcat,X)
  # N <- length(pik)
  n <- sum(s)
  index <- which(s == 1)
  H <- ncol(Xcat)
  q <- ncol(X)
  
  beta <- matrix(rep(0,(H+q)^2),ncol = (H+q),nrow = (H+q))
  #compute beta
  for(k in index){
    # print(k)
    z_k <- X_tmp[k,]
    c_k <- (1-pik[k])*(n/(n-(H+q)))
    beta <- beta + c_k*(z_k/pik[k])%*%t(z_k/pik[k])
  } 
  beta <- solve(beta)
  
  #compute tmp
  tmp <- rep(0,H+q)
  for(k in index){
    z_k <- X_tmp[k,]
    c_k <- (1-pik[k])*(n/(n-(H+q)))
    tmp <- tmp + c_k*y[k]/pik[k]*(z_k/pik[k])
  } 
  
  beta <- beta%*%tmp
  
  #compute v
  v <- 0
  for(k in index){
    z_k <- X_tmp[k,]
    c_k <- (1-pik[k])*(n/(n-(H+q)))
    v <- v + c_k*( (y[k]/pik[k]) - t(beta)%*%(z_k/pik[k]) )^2
  } 
  
  v
  return(v)
  
  
}
