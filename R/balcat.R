


#' Title
#'
#' @param X
#' @param Xcat
#' @param pik
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#' rm(list = ls())
#' N <- 100
#' pik <- rep(10/N,N)
#'
#' Xcat <-data.frame(cat1 = rep(1:50,each = 2))
#' Xcat_tmp <- disjMatrix(as.matrix(Xcat))
#' X <- cbind(pik,matrix(rnorm(N),ncol = 1))
#' A <- cbind(X,Xcat_tmp)/pik
#'
#' s <- balcat(X,Xcat,pik)
#' t(A)%*%s
#' t(A)%*%pik
#'
#'
#' rm(list = ls())
#' N <- 10000
#' pik <- rep(0.25,N)
#'
#' Xcat <- as.matrix(rep(1:1250,each = N/1250))
#' Xcat_tmp <- disjMatrix(Xcat)
#' X <- cbind(pik,matrix(rnorm(N),ncol = 1))
#' A <- cbind(X,Xcat_tmp)/pik
#'
#' system.time(s <- balcat(X,Xcat,pik))
#'
#'
#'
#' rm(list = ls())
#' N <- 1000
#' pik <- rep(0.25,N)
#'
#' Xcat <- as.matrix(rep(1:250,each = N/250))
#' Xcat_tmp <- disjMatrix(Xcat)
#' X <- cbind(pik,matrix(rnorm(N),ncol = 1))
#' A <- cbind(X,Xcat_tmp)/pik
#'
#'
#' system.time(s <- balcat(X,Xcat,pik))
#' system.time(s <- BalancedSampling::flightphase(pik,cbind(X,Xcat_tmp)))
#' system.time(s <- balancedstratification(X,Xcat,pik,comment=FALSE,method=2))
#'
#' microbenchmark( balcat(X,Xcat,pik),
#' BalancedSampling::flightphase(pik,cbind(X,Xcat_tmp)),
#' balancedstratification(X,Xcat,pik,comment=FALSE,method=2),times = 10)
#'
#'
#'
#' rm(list = ls())
#' N <- 120
#' pik <- rep(0.4,N)
#' pikInit <- pik
#' Xcat <-data.frame(cat1 = rep(1:40,each = 3))
#' Xcat_tmp <- disjMatrix(as.matrix(Xcat))
#' X <- cbind(pik,matrix(rnorm(N),ncol = 1))
#' X <- as.matrix(pik)
#' A <- cbind(X,Xcat_tmp)/pik
#'
#' s <- balcat(X,Xcat,pik)
#'
#' t(A)%*%s
#' t(A)%*%pik
#'
#'
#'
#' rm(list = ls())
#' N <- 10000
#' n <- 1000
#' x1 <- rgamma(N,4,25)
#' x2 <- rgamma(N,4,25)
#'
#' Xcat <- as.matrix(rep(1:n,each = N/n))
#'
#' pik <- inclusionprobastrata(Xcat,rep(1,n))
#' #pik <- rep(inclusionprobabilities(runif(400),1),25)
#'
#' X <- as.matrix(cbind(matrix(c(x1,x2),ncol = 2)))
#'
#' system.time(s <- balcat(X,Xcat,pik))
#'
balcat <- function(X,Xcat,pik){

  ##----------------------------------------------------------------
  ##                        Initialization                         -
  ##----------------------------------------------------------------
  Xcat <- as.matrix(Xcat)
  EPS = 1e-8
  A = X/pik

  p = ncol(X)
  pikInit <- pik
  n_all_cat <- sum(ncat(Xcat))


  ##----------------------------------------------------------------
  ##                  Flightphase on each strata                   -
  ##----------------------------------------------------------------

  H <- as.numeric(ncat(Xcat))
  EPS <- 1e-7

  ##----------------------------------------------------------------
  ##                  Flightphase on each strata                   -
  ##----------------------------------------------------------------

  for(k in 1:H){
    pik[Xcat == k] <- sampling::fastflightcube(as.matrix(X[Xcat == k,]),
                                                   pik[Xcat == k],
                                                   order = 1,
                                                   comment = FALSE)
  }


  ##----------------------------------------------------------------
  ##                Number of non 0-1 inclusion prob               -
  ##----------------------------------------------------------------

  i <- which(pik > EPS & pik < (1-EPS))
  i_size = length(i)


  ##----------------------------------------------------------------
  ##            flight phase with categorical variable             -
  ##----------------------------------------------------------------


  # while(i_size > p + n_all_cat){
  # step = 1
  while(i_size > 0){
    # print(step)
    ##------ Remove unique category

    uCat <- i[duplicated(Xcat[i,]) | duplicated(Xcat[i,], fromLast = TRUE)]
    if(length(uCat) == 0){
      break;
    }

    ##------ Find B
    A_tmp <- as.matrix(X[uCat,]/pikInit[uCat])
    # A_tmp <- as.matrix(X[uCat,]*pik[uCat]/pikInit[uCat])

    B <- findB(A_tmp,as.matrix(Xcat[uCat,]))
    B <- cbind(B$X,B$Xcat/pikInit[uCat[1:nrow(B$X)]])
    # print(dim(B))

    # B <- cbind(B$X,B$Xcat*pik[uCat[1:nrow(B$X)]])$

    # B <- B*pik[uCat[1:nrow(B)]]/pikInit[uCat[1:nrow(B)]]
    ##------ onestep and check if null

    tmp <-  onestep(B,pik[uCat[1:nrow(B)]],EPS)
    if(is.null(tmp)){
      break;
    }else{
      pik[uCat[1:nrow(B)]] <- tmp
    }

    ##------ update i

    i <- which(pik > EPS & pik < (1-EPS))
    i_size = length(i)
    # step = step + 1
  }

  # image(as(disjMatrix(as.matrix(Xcat[i,])),"sparseMatrix"))
  # A <- cbind(X,Xcat_tmp)/pikInit
  # Xcat_tmp <- disjMatrix(as.matrix(Xcat))
  # Atest <- cbind(X,Xcat_tmp)/pikInit
  # t(t(Atest)%*%pik)
  # t(t(Atest)%*%pikInit)

  ##----------------------------------------------------------------
  ##            end of flight phase without categories             -
  ##----------------------------------------------------------------


  while(i_size > 0){


    ##------ Find B

    if(i_size <= p){
      B <- as.matrix(A[i,]) # if not enough row
    }else{
      B <- as.matrix(A[i[1:(p+1)],]) # A has p columns
    }

    ##------ onestep and check if null

    tmp <-  onestep(B,pik[i[1:nrow(B)]],EPS)
    if(is.null(tmp)){
      break;
    }else{
      pik[i[1:nrow(B)]] <- tmp
    }

    ##------ update i

    i <- which(pik > EPS & pik < (1-EPS))
    i_size = length(i)

  }

  # head(t(AInit)%*%pik)
  # head(t(AInit)%*%pikInit)



  ##----------------------------------------------------------------
  ##                 LANDING by droping step by step               -
  ##----------------------------------------------------------------


  if(i_size > 1){
    # Xdev <- disjMatrix(as.matrix(Xcat[i,]))
    # Xland <- cbind(X[i,],Xdev)
    Xland <- X[i,]

    pikland = pik[i]
    Nland = length(pikland)
    nland = sum(pikland)
    p <- ncol(Xland)

    for(k in 0:(p-1)){
      # Xdev <- disjMatrix(as.matrix(Xcat[i,]))
      # Bland <- cbind(X[i,],Xdev)
      Bland <- X[i,]
      Bland <- Bland[,1:(p-k)]

      # pikstar[i] <- fastcube(Bland/pik[i]*pikland,pikland)
      kern <- MASS::Null(Bland)
      # kern
      if(length(kern)!=0){
        pik[i] <- onestep(Bland/pikInit[i]*pikland,pikland,EPS)
        i = which(pik > EPS & pik < (1 - EPS))
        pikland = pik[i]
        Nland = length(pikland)
        i_size <- length(i)
        # print(i_size)
      }
      if(i_size == 1){
        break;
      }
    }
  }

  if(length(i) > 1){
      stop("error not possible")
  }
  if(length(i) == 1){
      pik[i] <- rbinom(1,1,pik[i])
  }



  pik <- round(pik,10)
  # Xcat_tmp <- disjMatrix(as.matrix(Xcat))
  # Atest <- cbind(X,Xcat_tmp)/pikInit
  # t(t(Atest)%*%pik)
  # t(t(Atest)%*%pikInit)
  return(pik)
}
