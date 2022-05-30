#' @title Approximated variance for balanced sampling
#' 
#' @param X A matrix of size (\eqn{N} x \eqn{p}) of auxiliary variables on which the sample must be balanced.
#' @param strata A vector of integers that represents the categories.
#' @param pik A vector of inclusion probabilities.
#' @param y A variable of interest.
#'
#' @return a scalar, the value of the approximated variance.
#' @export
#' 
#' @details
#' 
#' This function gives an approximation of the variance of the Horvitz-Thompson total estimator presented by Hasler and Tillé (2014). 
#'
#' @references 
#' Hasler, C. and Tillé, Y. (2014). Fast balanced sampling for highly stratified population. \emph{Computational Statistics and Data Analysis}, 74:81-94.
#' 
#' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
#'
#' @seealso \code{\link{varEst}} 
#' 
#' @examples
#' 
#' N <- 100
#' n <- 10
#' p <- 5
#' 
#' Xaux <- matrix(rnorm(N*p),ncol = p)
#' 
#' s <- BalancedSampling::cube(pik,Xaux)
#' s_01 <- rep(0,N)
#' s_01[s] <- 1
#' pik <- rep(n/N,N)
#' y <- Xaux%*%c(1,1,1,1,1) + rnorm(N)
#' y <- rnorm(N)
#' varB(Xaux[s,],pik[s],y[s]) 
#' 
varB <- function(Xaux,pik,y){
  
  A <- Xaux/pik
  n <- length(pik)
  p <- ncol(Xaux)
  beta <- matrix(rep(0,p^2),ncol = p,nrow = p)
  
  #compute beta
  for(k in 1:n){
    c_k <- (1-pik[k])
    beta <- beta + c_k*(A[k,])%*%t(A[k,])
  } 
  # beta <- solve(beta)
  beta <- MASS::ginv(beta)
  #compute tmp
  tmp <- rep(0,p)
  for(k in 1:n){
    c_k <- (1-pik[k])
    tmp <- tmp + c_k*(y[k]/pik[k])*A[k,]
  } 
  beta <- beta%*%tmp
  
  #compute v
  v <- 0
  for(k in 1:n){
    c_k <- (1-pik[k])
    v <- v + c_k*( (y[k]/pik[k] - t(beta)%*%A[k,] ))^2
  } 
  
  v <- v*(n/(n-p))
  return(v)
  
}

#' @title Approximated variance for balanced sampling
#' 
#' @param X A matrix of size (\eqn{N} x \eqn{p}) of auxiliary variables on which the sample must be balanced.
#' @param strata A vector of integers that represents the categories.
#' @param pik A vector of inclusion probabilities.
#' @param y A variable of interest.
#'
#' @return a scalar, the value of the approximated variance.
#' @export
#' 
#' @details
#' 
#' This function gives an approximation of the variance of the Horvitz-Thompson total estimator presented by Hasler and Tillé (2014). 
#'
#' @references 
#' Hasler, C. and Tillé, Y. (2014). Fast balanced sampling for highly stratified population. \emph{Computational Statistics and Data Analysis}, 74:81-94.
#' 
#' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
#'
#' @seealso \code{\link{varEst}} 
#' 
#' @examples
#' 
#' rm(list = ls())
#' N <- 100
#' n <- 10
#' p <- 5
#' k <- 1
#' set.seed(k)
#' 
#' pik <- rep(n/N,N)
#' Xaux <- matrix(rnorm(N*p),ncol = p)
#' 
#' s <- BalancedSampling::cube(pik,Xaux)
#' s_01 <- rep(0,N)
#' s_01[s] <- 1
#' 
#' y <- Xaux%*%c(1,1,1,1,1) + rnorm(N)
#' 
#' vEst(Xaux[s,],pik[s],y[s])
#' varB(Xaux[s,],pik[s],y[s]) 
#' varB2(Xaux[s,],pik[s],y[s])
#' k = k + 1
#' 
varB2 <- function(Xaux,pik,y){
  # Xaux <- Xaux[s,]
  # pik <- pik[s]
  # y <- y[s]
  
  A <- Xaux/pik
  n <- length(pik)
  p <- ncol(Xaux)
  
  c <- (n/(n-p))*(1-pik)
  D <- diag(c)
  XX <- MASS::ginv(t(A)%*%D%*%A)
  
  proj <- A%*%XX%*%t(A)%*%D
  pred <- proj%*%(y/pik)
  diff <- y/pik - pred
  out <- t(diff)%*%D%*%diff 
  out
  return(out)
}


#' @title Approximated variance for balanced sampling
#' 
#' @param X A matrix of size (\eqn{N} x \eqn{p}) of auxiliary variables on which the sample must be balanced.
#' @param strata A vector of integers that represents the categories.
#' @param pik A vector of inclusion probabilities.
#' @param y A variable of interest.
#'
#' @return a scalar, the value of the approximated variance.
#' @export
#' 
#' @details
#' 
#' This function gives an approximation of the variance of the Horvitz-Thompson total estimator presented by Hasler and Tillé (2014). 
#'
#' @references 
#' Hasler, C. and Tillé, Y. (2014). Fast balanced sampling for highly stratified population. \emph{Computational Statistics and Data Analysis}, 74:81-94.
#' 
#' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
#'
#' @seealso \code{\link{varEst}} 
#' 
#' @examples
#' 
#' N <- 1000
#' n <- 100
#' p <- 5
#' 
#' Xaux <- matrix(rnorm(N*p),ncol = p)
#' Xspread <- as.matrix(cbind(runif(N),runif(N)))
#' pik <- rep(n/N,N)
#' Xaux <- cbind(pik,Xaux)
#' s <- BalancedSampling::cube(pik,Xaux)
#' y <- Xaux%*%c(1,1,1,1,1,1) + rnorm(N)
#' 
#' 
#' varB(Xaux[s,],pik[s],y[s]) 
#' varB2(Xaux[s,],pik[s],y[s])
#' varDBS(Xaux[s,],Xspread[s,],pik[s],y[s])
#' vDBS(Xaux[s,],Xspread[s,],pik[s],y[s])
#' 
varDBS <- function(Xaux,Xspread,pik,y){
  
  A <- Xaux/pik
  n <- length(pik)
  p <- ncol(Xaux)
  beta <- matrix(rep(0,p^2),ncol = p,nrow = p)
  # D <- as.matrix(proxy::dist(Xspread))
  
  
  #compute beta
  for(k in 1:n){
    b_k <- (1-pik[k])
    beta <- beta + b_k*(A[k,])%*%t(A[k,])
  } 
  beta <- MASS::ginv(beta) 
  
  #compute tmp
  tmp <- rep(0,p)
  for(k in 1:n){
    b_k <- (1-pik[k])
    tmp <- tmp + b_k*(y[k]/pik[k])*A[k,]
  } 
  beta <- beta%*%tmp
  
  e <- y - Xaux%*%beta
  e_bar <- rep(0,n)
  for(i in 1:n){
    dd <- distUnitk(Xspread,i,tore = FALSE,toreBound = -1)
    ii <- which(base::rank(dd,ties.method = "max") < (p + 1 + 1e-7 ))
    e_bar[i] <- sum((1-pik[ii])*e[ii]/pik[ii])/ sum((1-pik[ii]))
  }
  
  #compute v
  v <- 0
  for(k in 1:n){
    b_k <- (1-pik[k])
    v <- v + b_k*( e[k]/pik[k] - e_bar[k] )^2
  } 
  
  v <- (n/(n-p)) * ((p+1)/p) * v
  return(v)
  
}


