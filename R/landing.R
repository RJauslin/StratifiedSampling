#' @title Landing by linear programming
#'
#' @param X matrix of auxiliary variables on which the sample must be balanced.
#' @param pikstar vector of updated inclusion probabilities by the flight phase. See \code{\link{ffphase}}.
#' @param pik vector of inclusion probabilities.
#'
#' @description
#' This function does the landing phase of the cube method using linear programming.
#'
#' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
#'
#' @references
#' Chauvet, G. and Tille, Y. (2006). A fast algorithm of balanced sampling. Computational Statistics, 21/1:53-62
#'
#'
#' @export
#' @examples
#' N <- 1000
#' n <- 10
#' p <- 4
#' pik <- sampling::inclusionprobabilities(runif(N),n)
#' X <- cbind(pik,matrix(rnorm(N*p),ncol= p))
#' pikstar <- ffphase(X,pik) 
#' s <- landingLP(X,pikstar,pik)
#' sum(s)
#' t(X/pik)%*%pik
#' t(X/pik)%*%pikstar
#' t(X/pik)%*%s
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
#' @param pikstar vector of updated inclusion probabilities by the flight phase. See \code{\link{ffphase}}
#' @param pik vector of inclusion probabilities.
#'
#' @description
#' This function does the landing phase of the cube method using suppression of variables.
#'
#' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
#'
#'
#' @references
#' Chauvet, G. and Tille, Y. (2006). A fast algorithm of balanced sampling. Computational Statistics, 21/1:53-62
#'
#' @export
#' @examples
#' 
#' N <- 1000
#' n <- 10
#' p <- 4
#' pik <- sampling::inclusionprobabilities(runif(N),n)
#' X <- cbind(pik,matrix(rnorm(N*p),ncol= p))
#' pikstar <- ffphase(X,pik) 
#' s <- landingRM(X/pik,pikstar)
#' sum(s)
#' t(X/pik)%*%pik
#' t(X/pik)%*%pikstar
#' t(X/pik)%*%s
#' t(X/pik)%*%samplecube(X,pik)
landingRM <- function(X,pikstar){


  
  # Xcat_tmp3 <- as.matrix(Xnn[i,]*pik_tmp[i])
  # X <- as.matrix(cbind(Xcat_tmp3, X[i,]/pik[i]*pik_tmp[i]))
  # pikstar <- pik_tmp[i]
  # pik <- pik[i]
  ##----------------------------------------------------------------
  ##                          Initializing                         -
  ##----------------------------------------------------------------

  EPS = 1e-11
  N = nrow(X)
  i = which(pikstar > EPS & pikstar < (1 - EPS))
  i_size = length(i)
  
  pikland <- pikstar[i]
  # pikland = pikstar[i]
  # Xland <- X[i,]/pik[i]
  Xland <- X[i,]
  Nland = length(pikland)
  nland = sum(pikland)
  p <- ncol(Xland)
  
  
  j <-  which(pikland > EPS & pikland < (1 - EPS))
  j_size <- length(j)
  
  # t(Xland)%*%pikland
  # t(Xland)%*%pikstar[i]
  ##---------------------------------------------------------------
  ##                          Main loop                           -
  ##---------------------------------------------------------------
  
  
  for(k in 0:(p-1)){

    Bland <- Xland[j,]*pikland[j]
    Bland <- Bland[,1:(p-k)]

    kern <- MASS::Null(Bland)
    if(length(kern)!=0){
      
      pikland[j] <- ffphase(as.matrix(Bland),pikland[j])
      
      # pikstar[i] <- onestep(Bland/pik[i]*pikland,pikland,EPS)
      j = which(pikland > EPS & pikland < (1 - EPS))
      
      
      j_size <- length(j)
      # print(i_size)
    }
    if(j_size == 1){
      break;
    }

  }
  
  t(Xland)%*%pikland
  t(Xland)%*%pikstar[i]
  
  pikstar[i] = pikland
  i <- which(pikstar > EPS & pikstar < (1 - EPS))

  if(length(i) > 1){
    stop("error not possible")
  }else{
    pikstar[i] <- rbinom(1,1,pikstar[i])
  }

  return(round(pikstar,10))
}
