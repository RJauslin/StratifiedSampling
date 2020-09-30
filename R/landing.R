
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
#'
#' rm(list = ls())
#' library(sampling)
#'
#' N <- 100
#' n <- 10
#' Xcat <-as.matrix(data.frame(cat1 = rep(1:4,each = N/4),
#'                     cat2 = rep(1:5,each = N/5),
#'                     cat2 = rep(1:10,each = N/10)))
#'
#' p <- 3
#' X <- cbind(rep(1,N),matrix(rnorm(N*p),ncol = 3))
#' Xcat_tmp <- disjMatrix(Xcat)
#' Xred <- as.matrix(cbind(X,Xcat_tmp))
#'
#' pik <- rep(n/N,N)
#' A <- Xred/pik
#'
#' pikstar <- fastcube(X,pik,Xcat)
#' s <- landingLP(X,pikstar,pik,Xcat)
#' sum(s)
#' }
landingLP <- function(X,pikstar,pik,Xcat){

  ##----------------------------------------------------------------
  ##                          Initializing                         -
  ##----------------------------------------------------------------

  EPS = 1e-11
  N = nrow(X)
  i = which(pikstar > EPS & pikstar < (1 - EPS))
  Xdev <- disjMatrix(as.matrix(Xcat[i,]))

  pikland = pikstar[i]
  Xland <- cbind(A[i,],Xdev)
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
#'
#' rm(list = ls())
#'
#' library(sampling)
#'
#' N <- 100
#' Xcat <-as.matrix(data.frame(cat1 = rep(1:4,each = N/4),
#'                     cat2 = rep(1:5,each = N/5),
#'                     cat2 = rep(1:10,each = N/10)))
#' p <- 3
#' X <- cbind(rep(1,N),matrix(rnorm(N*p),ncol = 3))
#'
#' Xcat_tmp <- disjMatrix(Xcat)
#' Xcat_tmp <- do.call(cbind,apply(Xcat,MARGIN = 2,disjunctive))
#' Xred <- as.matrix(cbind(X,Xcat_tmp))
#'
#' pik <- rep(10/N,N)
#' A <- Xred/pik
#'
#'
#' # WITH Xcat
#'
#' pikstar <- fastcube(X,pik,Xcat)
#' s <- landingRM(X,pikstar,pik,Xcat)
#' s
#' sum(s)
#' t(t(Xred/pik)%*%pik)
#' t(t(Xred/pik)%*%pikstar)
#' t(t(Xred/pik)%*%s)
#' t(t(Xred/pik)%*%samplecube(Xred,pik,comment = FALSE,method = 2))
#'
#' # WITHOUT Xcat
#' pikstar <- fastcube(cbind(rep(1,nrow(X)),X),pik)
#' s <- landingRM(X,pikstar,pik)
#' s
#' sum(s)
#' t(cbind(rep(1,nrow(X)),X)/pik)%*%pik
#' t(cbind(rep(1,nrow(X)),X)/pik)%*%pikstar
#' t(cbind(rep(1,nrow(X)),X)/pik)%*%s
#' t(cbind(rep(1,nrow(X)),X)/pik)%*%samplecube(X,pik,comment = FALSE,method = 2)
#' }
landingRM <- function(X,pikstar,pik,Xcat){


  ##----------------------------------------------------------------
  ##                          Initializing                         -
  ##----------------------------------------------------------------

  EPS = 1e-11
  N = nrow(X)
  i = which(pikstar > EPS & pikstar < (1 - EPS))
  i_size = length(i)


  if(!missing(Xcat)){
    CAT = TRUE
  }else{
    CAT = FALSE
  }

  if(CAT == TRUE){
    Xdev <- disjMatrix(Xcat[i,])
    Xland <- cbind(X[i,],Xdev)
  }else{
    Xland <- X[i,]
  }


  pikland = pikstar[i]
  Nland = length(pikland)
  nland = sum(pikland)
  p <- ncol(Xland)

  for(k in 0:(p-1)){
    if(CAT == TRUE){
      Xdev <- disjMatrix(Xcat[i,])
      Bland <- cbind(X[i,],Xdev)
    }else{
      Bland <- X[i,]
    }
    Bland <- Bland[,1:(p-k)]

    # pikstar[i] <- fastcube(Bland/pik[i]*pikland,pikland)
    kern <- MASS::Null(Bland)
    # kern
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
