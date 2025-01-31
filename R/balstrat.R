#' @title Balanced Stratification
#' @name balstrat
#' @description 
#' 
#' Select a stratified balanced sample. The function is similar to \code{\link[sampling:balancedstratification]{balancedstratification}} of the package sampling.
#' 
#' @param X A matrix of size (\eqn{N} x \eqn{p}) of auxiliary variables on which the sample must be balanced.
#' @param strata A vector of integers that specifies the stratification.
#' @param pik A vector of inclusion probabilities.
#' @param rand if TRUE, the data are randomly arranged. Default TRUE
#' @param landing if TRUE, landing by linear programming otherwise supression of variables. Default TRUE
#'
#' @return A vector with elements equal to 0 or 1. The value 1 indicates that the unit is selected while the value 0 is for rejected units.
#'
#' @details The function implements the method proposed by Chauvet (2009). Firstly, a flight phase is performed on each strata. Secondly, a flight phase is applied on the whole population by aggregating the strata. Finally, a landing phase is applied by suppression of variables.
#' 
#' @references Chauvet, G. (2009). Stratified balanced sampling. \emph{Survey Methodology}, 35:115-119.
#'
#' @import Matrix
#' @useDynLib StratifiedSampling, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom sampling landingcube
#' 
#' @seealso \code{\link{ffphase}}, \code{\link{landingRM}}
#' 
#' 
#' @examples
#' 
#' N <- 100
#' n <- 10
#' p <- 4
#' X <- matrix(rgamma(N*p,4,25),ncol = p)
#' strata <- as.matrix(rep(1:n,each = N/n))
#' pik <- rep(n/N,N)
#' 
#' s <- balstrat(X,strata,pik)
#' 
#' t(X/pik)%*%s
#' t(X/pik)%*%pik
#' 
#' Xcat <- disj(strata)
#' 
#' t(Xcat)%*%s
#' t(Xcat)%*%pik
#' @export
balstrat <- function (X, strata, pik,rand = TRUE,landing = TRUE) 
{
 
  
  if(rand == TRUE){
    strataInit <- strata
    ## cleanstrata
    a = sort(unique(as.vector(strata)))
    b = 1:length(a)
    names(b) = a
    strata <- as.vector(b[as.character(strata)])
    
    ## RANDOMIZATION 
    o <- order(strata)
    # strata[o] # this order the vector
    
    o_split <- split(o, f = strata[o] ) # split and randomize in each category
    for( i in 1:length(o_split)){
      o_tmp <- sample(1:length(o_split[[i]]))
      o_split[[i]] <- o_split[[i]][o_tmp]
    }
    o_out <- unlist(o_split[sample(1:length(o_split))]) # unlist and randomize each category
    
    XInit <- X
    pikInit <- pik  
    
    X <- as.matrix(X[o_out,])
    strata <- as.matrix(strata[o_out])
    pik <- pik[o_out]
  }
  
  H <- as.numeric(ncat(as.matrix(strata)))
  pik_tmp <- pik
  EPS <- 1e-7
  Xnn <- disj(strata)
  
  ##----------------------------------------------------------------
  ##                  Flightphase on each strata                   -
  ##----------------------------------------------------------------
  
  
  for(k in 1:H){
    pik_tmp[strata == k] <- ffphase(as.matrix(cbind(pik[which(strata == k)],as.matrix(X[which(strata == k),]))),
                                    pik[strata == k])
    
  }
  
  # ###################### CHECK
  # sum(pik_tmp)
  # t(X/pik)%*%pik_tmp
  # t(X/pik)%*%pik
  # 
  # t(Xnn)%*%pik
  # t(Xnn)%*%pik_tmp
  
  ##----------------------------------------------------------------
  ##          Flightphase on the uninon of strata U1 -- Uk         -
  ##----------------------------------------------------------------
  i <- which(pik_tmp > EPS & pik_tmp < (1-EPS))
  if(length(i) != 0){
    
    # Xcat_tmp2 <- disjMatrix(as.matrix(strata[i,]))
    Xcat_tmp2 <- disj(strata[i])
    Xcat_tmp2 <- Xcat_tmp2*pik_tmp[i]
    
    X_tmp <- as.matrix((X[i,]*pik_tmp[i]/pik[i]))
    pik_tmp[i] <- ffphase(as.matrix(cbind(X_tmp,Xcat_tmp2)),pik_tmp[i])
    
  }
  
  
  # ###################### CHECK
  # sum(pik_tmp)
  # t(X/pik)%*%pik_tmp
  # t(X/pik)%*%pik
  # # 
  # t(Xnn)%*%pik
  # t(Xnn)%*%pik_tmp
  
  ##---------------------------------------------------------------
  ##              Landing by suppression of variables             -
  ##---------------------------------------------------------------
  
  i <- which(pik_tmp > EPS & pik_tmp < (1-EPS))
  
  
  if(length(i) > 0){
    if(landing == TRUE){
      if(length(i) > 20){
        warnings("The landing by using linear programming might be very time consuming. Think about landing by using drop variables.")
      }
      pik_tmp <- sampling::landingcube(cbind(pik,Xnn,X),pik_tmp,pik,comment = FALSE) # pass the whole matrix to compute t(A)%*%A
    }else{
      pik_tmp[i] <- landingRM(cbind(pik[i],Xnn[i,],X[i,])/pik[i]*pik_tmp[i],pik_tmp[i],EPS)
    }
  }
  
  
  if(rand == TRUE){
    s_01 <- rep(0,length(pik_tmp))
    s_01[o_out[which(pik_tmp > (1-EPS))]] <- 1
  }else{
    s_01 <- round(pik_tmp,10)
  }
  
  
  return(s_01)
  
}