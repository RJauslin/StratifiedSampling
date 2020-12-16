#' @title Stratified Sampling
#' 
#' @description  
#' sndgj
#'
#' @param X matrix of auxiliary variables.
#' @param strata matrix of categorical variables.
#' @param pik vector of inclusion probabilities.
#'
#' @return A vector with elements equal to 0 or 1. The value 1 indicates that the unit is selected while the value 0 is for rejected units.
#' @export
#'
#' @examples
#'
#' ##----------------------------------------------------------------
#' ##                              Data                             -
#' ##----------------------------------------------------------------
#' 
#' rm(list = ls())
#' N <- 1000
#' n <- 100
#' x1 <- rgamma(N,4,25)
#' x2 <- rgamma(N,4,25)
#' strata <- as.matrix(rep(1:n,each = N/n))
#' Xcat <- disjMatrix(strata)
#' 
#' 
#' 
#' ##---------------------------------------------------------------
#' ##                        different cases                       -
#' ##---------------------------------------------------------------
#'   
#' 
#' #-------- CASE 1 pik equal and integer in each strata 
#' 
#'  pik <- sampling::inclusionprobastrata(strata,rep(1,n))
#'  X <- as.matrix(cbind(pik,matrix(c(x1,x2),ncol = 2)))
#'  system.time(s <- stratifiedcube(X,strata,pik))
#'  sum(s)
#'  t(X/pik)%*%s
#'  t(X/pik)%*%pik
#'  
#'  t(Xcat)%*%s
#'  t(Xcat)%*%pik
#'  
#' 
#' #-------- CASE 2 pik unequal and integer in each strata  
#' 
#'  pik <- rep(sampling::inclusionprobabilities(runif(N/n),1),n)
#'  X <- as.matrix(cbind(pik,matrix(c(x1,x2),ncol = 2)))
#'  system.time(s <- stratifiedcube(X,strata,pik))
#'  sum(s)
#'  t(X/pik)%*%s
#'  t(X/pik)%*%pik
#'  
#'  t(Xcat)%*%s
#'  t(Xcat)%*%pik  
#'
#' #-------- CASE 3 pik equal and NOT integer in each strata 
#'  
#'  pik <- rep(0.25,N)
#'  X <- as.matrix(cbind(pik,matrix(c(x1,x2),ncol = 2)))
#'  system.time(s <- stratifiedcube(X,strata,pik))
#'  sum(s)
#'  t(X/pik)%*%s
#'  t(X/pik)%*%pik
#'  
#'  t(Xcat)%*%s
#'  t(Xcat)%*%pik
#'  
#' #-------- CASE 4 pik unequal and NOT integer in each strata 
#'  
#' pik <- rep(sampling::inclusionprobabilities(runif(N/n),2.5),n)
#' X <- as.matrix(cbind(pik,matrix(c(x1,x2),ncol = 2)))
#' system.time(s <- stratifiedcube(X,strata,pik))
#' sum(s)
#' t(X/pik)%*%s
#' t(X/pik)%*%pik
#' 
#' t(Xcat)%*%s
#' t(Xcat)%*%pik
stratifiedcube <- function(X,strata,pik){
  
  ##----------------------------------------------------------------
  ##                        Initialization                         -
  ##----------------------------------------------------------------
  strata <- as.matrix(strata)
  EPS = 1e-8
  A = X/pik
  pikstar <- pik
  p = ncol(X)
  nstrata <- length(unique(strata))
  
  ##----------------------------------------------------------------
  ##                  Flightphase on each strata                   -
  ##----------------------------------------------------------------
  
  for(k in 1:nstrata){
    # pik[Xcat == k] <- sampling::fastflightcube(as.matrix(X[Xcat == k,]),
    #                                            pik[Xcat == k],
    #                                            order = 1,
    #                                            comment = FALSE)
    #
    # pikstar[strata == k] <-ffphase(cbind(pik[which(strata == k)],as.matrix(X[which(strata == k),])),
    #        pikstar[strata == k])
    pikstar[strata == k] <-ffphase(as.matrix(cbind(pikstar[strata == k],X[which(strata == k),])),
                                   pikstar[strata == k])
  }
  
  ###################### CHECK
  # t(X/pik)%*%pikstar
  # t(X/pik)%*%pik
  # t(Xcat)%*%pik
  # t(Xcat)%*%pikstar
  
  ##----------------------------------------------------------------
  ##                Number of non 0-1 inclusion prob               -
  ##----------------------------------------------------------------
  
  i <- which(pikstar > EPS & pikstar < (1-EPS))
  i_size = length(i)
  i_size
  
  
  ##----------------------------------------------------------------
  ##            flight phase with categorical variable             -
  ##----------------------------------------------------------------
  
  if(i_size > 0 ){
    while(i_size > 0){
    
      ##------ Remove unique category
      
      uCat <- i[duplicated(strata[i,]) | duplicated(strata[i,], fromLast = TRUE)]
      if(length(uCat) == 0){
        break;
      }
      
      ##------ Find B
      A_tmp <- as.matrix(X[uCat,]/pik[uCat])
      B <- findB(A_tmp,as.matrix(strata[uCat,]))
      # B <- cbind(B$X,B$Xcat/pik[uCat[1:nrow(B$X)]])
      B <- cbind(B$X,B$Xcat)
      # print(B)
    
      ##------ onestep and check if null
      tmp <-  onestep(B,pikstar[uCat[1:nrow(B)]],EPS)
      if(is.null(tmp)){
        break;
      }else{
        pikstar[uCat[1:nrow(B)]] <- tmp
      }
      
      ##------ update i
      
      i <- which(pikstar > EPS & pikstar < (1-EPS))
      i_size = length(i)
      # print(sum(pikstar))
    }
    
    
    # ##----------------------------------------------------------------
    # ##          end of flight phase on strata categories             -
    # ##----------------------------------------------------------------
    # Sometimes some stata does could not be balanced at the end and so we 
    # drop some auxiliary variable such that we could have only one unit that
    # are not put equal to 0 or 1 within each strata
    #
    
    p <- ncol(X)
    k = 1
    while(length(uCat) != 0){
      ##------ Find B
      if(k == p){
        B <- as.matrix(pikstar[uCat])/pikstar[uCat]
      }else{
        A_tmp <- as.matrix(as.matrix(X[uCat,1:(p-k)])/pik[uCat])
        B <- findB(A_tmp,as.matrix(strata[uCat,]))
        B <- cbind(B$X,B$Xcat)
      }
      
      # print(pikstar[uCat[1:nrow(B)]])
      tmp <-  onestep(B,pikstar[uCat[1:nrow(B)]],EPS)
      # print(tmp)
      if(!is.null(tmp)){
        pikstar[uCat[1:nrow(B)]] <- tmp  
      }
      i <- which(pikstar > EPS & pikstar < (1-EPS))
      i_size = length(i)
      uCat <- i[duplicated(strata[i,]) | duplicated(strata[i,], fromLast = TRUE)]
      k = k + 1
    }
  }
  
  #---------------- check
  # Xcat <- disjMatrix(strata)
  # 
  # t(X/pik)%*%pik
  # t(X/pik)%*%pikstar
  # t(Xcat)%*%pik
  # t(Xcat)%*%pikstar
  # length(i)
  # sum(pikstar)
  
  
  # ##----------------------------------------------------------------
  # ##            Landing on unit that are alone in the strata       -
  # ##----------------------------------------------------------------
  if(length(i) > 0){
    pikstar[i] <- landingRM(cbind(pikstar[i],X[i,]/pik[i]),pikstar[i])    
  }

  
  return(round(pikstar,10))
}

