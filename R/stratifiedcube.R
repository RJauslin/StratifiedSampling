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
#'
#' @importFrom sampling inclusionprobastrata
#' @examples
#' N <- 100
#' n <- 10
#' p <- 4
#' X <- matrix(rgamma(N*p,4,25),ncol = p)
#' strata <- as.matrix(rep(1:n,each = N/n))
#' pik <- rep(n/N,N)
#' 
#' s <- stratifiedcube(X,strata,pik)
#' 
#' t(X/pik)%*%s
#' t(X/pik)%*%pik
#' 
#' Xcat <- disj(strata)
#' 
#' t(Xcat)%*%s
#' t(Xcat)%*%pik
#' 
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
      B <- cbind(B$X,B$Xcat)
    
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

