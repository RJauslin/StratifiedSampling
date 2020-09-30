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
#' @export
#' @examples
#' \dontrun{
#'
#'
#' ########################################################
#'
#'
#' rm(list = ls())
#'
#' library(sampling)
#'
#' N <- 1000
#'
#' Xcat <-as.matrix(data.frame(cat1 = rep(1:40,each = N/40),
#'                     cat2 = rep(1:50,each = N/50),
#'                     cat2 = rep(1:100,each = N/100)))
#'
#'
#' p <- 30
#' X <- matrix(rnorm(N*p),ncol = 30)
#'
#'
#' Xcat_tmp <- disjMatrix(Xcat)
#' Xcat_tmp <- do.call(cbind,apply(Xcat,MARGIN = 2,disjunctive))
#' Xred <- as.matrix(cbind(X,Xcat_tmp))
#'
#' pik <- rep(300/N,N)
#' A <- Xred/pik
#'
#'
#' system.time(s1 <- fastcube(Xred,pik))
#' system.time(s1 <- ffphase(X,pik,Xcat))
#' system.time(s1 <- fastcube(X,pik))
#' as.vector(t(A)%*%s1)
#' as.vector(t(A)%*%pik)
#'
#' system.time(s2 <- ReducedSamplecube(Xred,pik,redux = TRUE))
#' A <- Xred/pik
#' as.vector(t(A)%*%s2)
#' as.vector(t(A)%*%pik)
#'
#' system.time(s3 <- fastflightcube(Xred,pik,order = 2,comment = FALSE))
#' as.vector(t(A)%*%s3)
#' as.vector(t(A)%*%pik)
#'
#'
#'
#' system.time(s4 <- BalancedSampling::flightphase(pik,Xred))
#' as.vector(t(A)%*%s4)
#' as.vector(t(A)%*%pik)
#'
#'
#'
#' rm(list = ls())
#' N <- 100
#' Xcat <-data.frame(cat1 = rep(1:50,each = 2))
#' pik <- rep(10/N,N)
#' X <- cbind(pik,matrix(rnorm(N),ncol = 1))
#' sum(ffphase(X,pik,Xcat))
#'
#'
#' rm(list = ls())
#' N <- 10000
#' Xcat <- as.matrix(rep(1:1250,each = N/1250))
#' Xcat_tmp <- disjMatrix(Xcat)
#' sum(ncat(Xcat))
#' pikInit <- rep(0.25,N)
#' X <- cbind(pikInit,matrix(rnorm(N),ncol = 1))
#' Xconc <- as.matrix(cbind(X,Xcat_tmp))
#'
#'
#' AInit <- Xconc/pikInit
#' pik <- pikInit
#' system.time(s <- ffphase(X,pik,Xcat))
#'
#' system.time(sampling::balancedstratification(X,Xcat,pik,comment=TRUE,method=1))
#'
#' head(t(AInit)%*%pikInit)
#' head(t(AInit)%*%s)
#' system.time(s <- BalancedSampling::flightphase(pik,Xconc))
#'
#'
#' }
ffphase <- function(X, pik, Xcat){

  ##----------------------------------------------------------------
  ##                        Initialization                         -
  ##----------------------------------------------------------------

  Xcat <- as.matrix(Xcat)

  EPS = 1e-8
  A = X/pik
  p = ncol(X)
  if(!missing(Xcat)){
    CAT = TRUE
    n_all_cat <- sum(ncat(Xcat))
  }else{
    CAT = FALSE
    n_all_cat <- 0
  }



  ##----------------------------------------------------------------
  ##                Number of non 0-1 inclusion prob               -
  ##----------------------------------------------------------------

  i <- which(pik > EPS & pik < (1-EPS))
  i_size = length(i)

  ##---------------------------------------------------------------
  ##                          Main loop                           -
  ##---------------------------------------------------------------


  while(i_size > p + n_all_cat){
    uCat <- i[duplicated(Xcat[i,]) | duplicated(Xcat[i,], fromLast = TRUE)]


    # cat("i_size : ",i_size,"\n\n")

    ##-----------------------------------
    ##  Find B
    ##-----------------------------------

    # if(i_size <= p + n_all_cat){
    #   Xdev <- disjMatrix(as.matrix(Xcat[i,]))
    #   B <- cbind(A[i,],Xdev)
    # }else{
      B <- findB(as.matrix(A[uCat,]),as.matrix(Xcat[uCat,]))
      # Sys.sleep(3)
      # print(image(as(B,"sparseMatrix")))
      # print(dim(B))
    # }
    tmp <-  onestep(B,pik[uCat[1:nrow(B)]],EPS)
    if(is.null(tmp)){
      break;
    }else{
      pik[uCat[1:nrow(B)]] <- tmp
    }



    # cat("i_size : ",i_size,"\n\n")
    # ##-----------------------------------
    # ##  Find B
    # ##-----------------------------------
    #
    # if(i_size <= p + n_all_cat){
    #   if(CAT == TRUE){
    #     Xdev <- disjMatrix(as.matrix(Xcat[i,]))
    #     B <- cbind(A[i,],Xdev)
    #   }else{
    #     B <- A[i,]
    #   }
    # }else{
    #   if(CAT == TRUE){
    #     B <- findB(A[i,],Xcat[i,])
    #
    #     # Sys.sleep(3)
    #     # print(image(as(B,"sparseMatrix")))
    #     print(dim(B))
    #   }else{
    #     B <- A[i[1:(p+1)],]
    #   }
    # }
    # tmp <-  onestep(B,pik[i[1:nrow(B)]],EPS)
    # if(is.null(tmp)){
    #   break;
    # }else{
    #   pik[i[1:nrow(B)]] <- tmp
    # }


    ##------------
    ##  update i
    ##------------

    i <- which(pik > EPS & pik < (1-EPS))
    i_size = length(i)

  }



  image(as(disjMatrix(as.matrix(Xcat[i,])),"sparseMatrix"))
  head(t(AInit)%*%pik)
  head(t(AInit)%*%pikInit)

  while(i_size > 0){

    ##-----------------------------------
    ##  Find B
    ##-----------------------------------


    if(i_size <= p){
      B <- as.matrix(A[i,])
    }else{
      B <- as.matrix(A[i[1:(p+1)],])
    }
    tmp <-  onestep(B,pik[i[1:nrow(B)]],EPS)
    if(is.null(tmp)){
      break;
    }else{
      pik[i[1:nrow(B)]] <- tmp
    }

    ##------------
    ##  update i
    ##------------

    i <- which(pik > EPS & pik < (1-EPS))
    i_size = length(i)

  }

  head(t(AInit)%*%pik)
  head(t(AInit)%*%pikInit)


  return(pik)
}








