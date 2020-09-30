
#' @title Landing by linear programming
#'
#' @param X matrix of auxiliary variables on which the sample must be balanced.
#' @param pikstar vector of updated inclusion probabilities by the flight phase.
#' @param pik vector of inclusion probabilities.
#' @param Xcat matrix of categorical variable (not in its disjunctive form).
#'
#' @description
#' This function does the landing phase of the cube method using linear programming.
#'  It allows you to put some categorical variable using the \code{Xcat} variable.
#'
#' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
#'
#'
#' @references
#' Chauvet, G. and Tille, Y. (2006). A fast algorithm of balanced sampling. Computational Statistics, 21/1:53-62
#'
#'
#' @export
#' @examples
#' \dontrun{
#' rm(list = ls())
#' N <- 100
#' n <- 10
#' p <- 4
#' 
#' pik <- sampling::inclusionprobabilities(runif(N),n)
#' 
#' X <- cbind(pik,matrix(rnorm(N*p),ncol= p))
#' pikstar <- ffphase(X,pik) 
#' s <- landingLP(X,pikstar,pik)
#' sum(s)
#'
#' }
landingLP <- function(X,pikstar,pik){

  ##----------------------------------------------------------------
  ##                          Initializing                         -
  ##----------------------------------------------------------------

  EPS = 1e-7
  N = nrow(X)
  i = which(pikstar > EPS & pikstar < (1 - EPS))
  pikland = pikstar[i]
  Xland <- X[i,]
  Nland = length(pikland)
  nland = sum(pikland)


  ##---------------------------------------------------------------
  ##            Calculate sampleSet and sampleSetSize             -
  ##---------------------------------------------------------------

  FLAGI = (abs(nland - round(nland)) < EPS)
  if(FLAGI){
    pikland = round(nland) * pikland/nland
    nland = round(nland)
    sampleSet = samplen(Nland,nland)
  }else{
    sampleSet = cbind(samplen(Nland,trunc(nland)), samplen(Nland,trunc(nland) + 1))
  }
  sampleSetSize = ncol(sampleSet)

  ##----------------------------------------------------------------
  ##                        Calculate cost                         -
  ##----------------------------------------------------------------


  Astar <-  t(Xland/pik[i]) %*%(sampleSet-pikland)
  A = Xland/pik[i]
  cost = apply(Astar,
               MARGIN = 2,
               FUN = function(x,H){return(t(x)%*%H%*%x)},
               H =  MASS::ginv(t(A) %*% A))

  V = sampleSet
  b = pikland
  constdir = rep("==", times = (Nland))
  # x = lpSolve::lp("min", rep(1,length(cost)), V, constdir, b)
  x = lpSolve::lp("min", cost, V, constdir, b)
  if(x$status == 2){
    stop("Error: no feasible solution reached.")
  }else{
    x <- x$solution
  }


  ##---------------------------------------------------------------
  ##                  Choose a sampleSet randomly                 -
  ##---------------------------------------------------------------


  u = runif(1, 0, 1)
  s = 0
  ccc = 0
  while (ccc < u) {
    s = s + 1
    ccc = ccc + x[s]
  }


  ##----------------------------------------------------------------
  ##                        update pikstar                         -
  ##----------------------------------------------------------------

  pikfin = pikstar
  pikfin[i] = sampleSet[,s]

  return(round(pikfin,10))

}




#' @title Landing by suppression of variables
#'
#' @param X matrix of auxiliary variables on which the sample must be balanced.
#' @param pikstar vector of updated inclusion probabilities by the flight phase.
#' @param pik vector of inclusion probabilities.
#' @param Xcat matrix of categorical variable (not in its disjunctive form).
#'
#' @description
#' This function does the landing phase of the cube method using suppression of variables.
#'  It allows you to put some categorical variable using the \code{Xcat} variable.
#'
#' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
#'
#'
#' @references
#' Chauvet, G. and Tille, Y. (2006). A fast algorithm of balanced sampling. Computational Statistics, 21/1:53-62
#'
#' @export
#' @examples
#' \dontrun{
#' rm(list = ls())
#' N <- 100
#' n <- 10
#' p <- 4
#' 
#' pik <- sampling::inclusionprobabilities(runif(N),n)
#' 
#' X <- cbind(pik,matrix(rnorm(N*p),ncol= p))
#' pikstar <- ffphase(X,pik) 
#' s <- landingRM(X,pikstar,pik)
#' sum(s)
#' t(X/pik)%*%pik
#' t(X/pik)%*%pikstar
#' t(X/pik)%*%s
#' }
landingRM <- function(X,pikstar,pik){


  ##----------------------------------------------------------------
  ##                          Initializing                         -
  ##----------------------------------------------------------------

  EPS = 1e-11
  N = nrow(X)
  i = which(pikstar > EPS & pikstar < (1 - EPS))
  i_size = length(i)
  Xland <- X[i,]
  pikland = pikstar[i]
  Nland = length(pikland)
  nland = sum(pikland)
  p <- ncol(Xland)
  
  ##---------------------------------------------------------------
  ##                          Main loop                           -
  ##---------------------------------------------------------------
  
  
  for(k in 0:(p-1)){

    Bland <- X[i,]
    Bland <- Bland[,1:(p-k)]

    
    kern <- MASS::Null(Bland)

    if(length(kern)!=0){
      pikstar[i] <- onestep(Bland/pik[i]*pikland,pikland,EPS)
      i = which(pikstar > EPS & pikstar < (1 - EPS))
      pikland = pikstar[i]
      Nland = length(pikland)
      i_size <- length(i)
      # print(i_size)
    }
    if(i_size == 1){
      break;
    }

  }

  if(length(i) > 1){
    stop("error not possible")
  }else{
    pikstar[i] <- rbinom(1,1,pikstar[i])
  }

  return(round(pikstar,10))
}
