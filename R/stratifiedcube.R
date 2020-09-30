


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
#' ##----------------------------------------------------------------
#' ##                              Data                             -
#' ##----------------------------------------------------------------
#' 
#' rm(list = ls())
#' N <- 10000
#' n <- 1000
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
#'  #-------- CASE 0 pik equal and only pik as variable
#'  
#'  pik <- inclusionprobastrata(strata,rep(1,n))
#'  X <- as.matrix(pik)
#'  system.time(s <- stratifiedcube(X,strata,pik))
#'  sum(s)
#'  t(X/pik)%*%s
#'  t(X/pik)%*%pik
#'  
#'  t(Xcat)%*%s
#'  t(Xcat)%*%pik
#'  
#' #-------- CASE 0.1 pik equal and only pik as variable
#'  
#'  pik <- rep(inclusionprobabilities(runif(N/n),1),n)
#'  X <- as.matrix(pik)
#'  system.time(s <- stratifiedcube(X,strata,pik))
#'  sum(s)
#'  t(X/pik)%*%s
#'  t(X/pik)%*%pik
#'  
#'  t(Xcat)%*%s
#'  t(Xcat)%*%pik
#'   
#' 
#' #-------- CASE 1 pik equal and integer in each strata 
#' 
#'  pik <- inclusionprobastrata(strata,rep(1,n))
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
#'  pik <- rep(inclusionprobabilities(runif(N/n),1),n)
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
#'  pik <- rep(inclusionprobabilities(runif(N/n),2.5),n)
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
    pikstar[strata == k] <-ffphase(as.matrix(X[which(strata == k),]),
                                   pikstar[strata == k])
  }
  
  
  ##----------------------------------------------------------------
  ##                Number of non 0-1 inclusion prob               -
  ##----------------------------------------------------------------
  
  i <- which(pikstar > EPS & pikstar < (1-EPS))
  i_size = length(i)
  
  
  ##----------------------------------------------------------------
  ##            flight phase with categorical variable             -
  ##----------------------------------------------------------------
  
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
  
  #---------------- check
  # Xcat <- disjMatrix(strata)
  # 
  # t(X/pik)%*%pik
  # t(X/pik)%*%pikstar
  # 
  # t(Xcat)%*%pik
  # t(Xcat)%*%pikstar
  
  ##----------------------------------------------------------------
  ##            end of flight phase without categories             -
  ##----------------------------------------------------------------

  if(i_size > 0){
    pikstar[i] <- ffphase(cbind(pikstar[i],X[i,]),pikstar[i])  
    i <- which(pikstar > EPS & pikstar < (1-EPS))
    i_size <- length(i)
  }
  
  
  
  ##----------------------------------------------------------------
  ##                 LANDING by droping step by step               -
  ##----------------------------------------------------------------
  if(i_size > 0){
  pikstar[i] <- landingRM(X[i,],
            pikstar[i],
            pik[i])
  }
  
  pikstar <- round(pikstar,10)
 
  return(pikstar)
}

