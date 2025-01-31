#' @title Stratified Sampling
#' @name stratifiedcube
#' @description  
#' 
#' This function implements a method for selecting a stratified sample. It really improves the performance of the function \code{\link{fbs}} and \code{\link{balstrat}}.
#'
#' @param X A matrix of size (\eqn{N} x \eqn{p}) of auxiliary variables on which the sample must be balanced.
#' @param strata A vector of integers that specifies the stratification..
#' @param pik A vector of inclusion probabilities.
#' @param EPS epsilon value
#' @param landing if FALSE, no landing phase is done.
#' @param rand if TRUE, the data are randomly arranged. Default TRUE
#' @param lp if TRUE, landing by linear programming otherwise supression of variables. Default TRUE
#'
#' @details 
#' 
#' The function is selecting a balanced sample very quickly even if the sum of inclusion probabilities within strata are non-integer. The function should be used in preference. Firstly, a flight phase is performed on each strata. Secondly, the function \code{\link{findB}} is used to find a particular matrix to apply a flight phase by using the cube method proposed by Chauvet, G. and Tillé, Y. (2006). Finally, a landing phase is applied by suppression of variables.
#' 
#' @seealso \code{\link{fbs}}, \code{\link{balstrat}}, \code{\link{landingRM}}, \code{\link{ffphase}}
#'
#' @return A vector with elements equal to 0 or 1. The value 1 indicates that the unit is selected while the value 0 is for rejected units.
#' @export
#'
#' @references 
#' Chauvet, G. and Tillé, Y. (2006). A fast algorithm of balanced sampling. \emph{Computational Statistics}, 21/1:53-62
#'
#' @examples
#' 
#' # EXAMPLE WITH EQUAL INCLUSION PROBABILITES AND SUM IN EACH STRATA INTEGER
#' N <- 100
#' n <- 10
#' p <- 4
#' X <- matrix(rgamma(N*p,4,25),ncol = p)
#' strata <- rep(1:n,each = N/n)
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
#' 
#' # EXAMPLE WITH UNEQUAL INCLUSION PROBABILITES AND SUM IN EACH STRATA INTEGER
#' N <- 100
#' n <- 10
#' X <- cbind(rgamma(N,4,25),rbinom(N,20,0.1),rlnorm(N,9,0.1),runif(N))
#' colSums(X)
#' strata <- rbinom(N,10,0.7)
#' strata <- sampling::cleanstrata(strata)
#' pik <- as.vector(sampling::inclusionprobastrata(strata,ceiling(table(strata)*0.10)))
#' EPS = 1e-7
#' 
#' s <- stratifiedcube(X,strata,pik)
#' test <- stratifiedcube(X,strata,pik,landing = FALSE)
#' 
#' t(X/pik)%*%s
#' t(X/pik)%*%test
#' t(X/pik)%*%pik
#' 
#' Xcat <- disj(strata)
#' 
#' t(Xcat)%*%s
#' t(Xcat)%*%test
#' t(Xcat)%*%pik
#' 
#' 
#' # EXAMPLE WITH UNEQUAL INCLUSION PROBABILITES AND SUM IN EACH STRATA NOT INTEGER
#' set.seed(3)
#' N <- 100
#' n <- 10
#' X <- cbind(rgamma(N,4,25),rbinom(N,20,0.1),rlnorm(N,9,0.1),runif(N))
#' strata <- rbinom(N,10,0.7)
#' strata <- sampling::cleanstrata(strata)
#' pik <- runif(N)
#' EPS = 1e-7
#' tapply(pik,strata,sum)
#' table(strata)
#' 
#' 
#' s <- stratifiedcube(X,strata,pik,landing = TRUE)
#' test <- stratifiedcube(X,strata,pik,landing = FALSE)
#' 
#' 
#' t(X/pik)%*%s
#' t(X/pik)%*%test
#' t(X/pik)%*%pik
#' 
#' Xcat <- disj(strata)
#' 
#' t(Xcat)%*%s
#' t(Xcat)%*%pik
#' t(Xcat)%*%test
#' 
#' 
stratifiedcube <- function(X,
                           strata,
                           pik,
                           EPS = 1e-7,
                           rand = TRUE,
                           landing = TRUE,
                           lp = TRUE){
  
  if(rand == TRUE){
    strataInit <- strata
    ## cleanstrata
    a = sort(unique(strata))
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
  
  ##----------------------------------------------------------------
  ##                        Initialization                         -
  ##----------------------------------------------------------------
  # strata <- as.matrix(strata)
  A = X/pik
  pikstar <- pik
  p = ncol(X)
  nstrata <- length(unique(strata))
  # TOT <- t(A)%*%pik
  # XX <- MASS::ginv(t(A) %*% A)
  
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
  # Xcat <- disj(strata)
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
      uCat <- i[duplicated(strata[i]) | duplicated(strata[i], fromLast = TRUE)]
      if(length(uCat) == 0){
        break;
      }
      
      ##------ Find B
      A_tmp <- X[uCat,]/pik[uCat]
      B <- findB(A_tmp,strata[uCat])
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
    
      if(landing == TRUE){
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
            A_tmp <- X[uCat,1:(p-k)]/pik[uCat]
            B <- findB(A_tmp,strata[uCat])
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
          uCat <- i[duplicated(strata[i]) | duplicated(strata[i], fromLast = TRUE)]
          k = k + 1
        }
      }
  }
    
    #---------------- check
    # Xcat <- disjMatrix(strata)
    # 
    # print(t(X/pik)%*%pik)
    # print(t(X/pik)%*%pikstar)
    # print(t(Xcat)%*%pik)
    # print(t(Xcat)%*%pikstar)
    # print(length(i))
    # print(sum(pikstar))
  
  
    # ##----------------------------------------------------------------
    # ##            Landing on unit that are alone in the strata       -
    # ##----------------------------------------------------------------
    if(landing == TRUE){
      i <- which(pikstar > EPS & pikstar < (1-EPS))
      if(length(i) > 0){
        if(lp == TRUE){
          if(length(i) > 20){
            warning("The landing by using linear programming might be very time consuming. Think about landing by using drop variables.")
          }
          pikstar <- sampling::landingcube(X,pikstar,pik,comment = FALSE) # pass the whole matrix to compute t(A)%*%A
        }else{
          pikstar[i] <- landingRM(cbind(pik[i],X[i,])/pik[i]*pikstar[i],pikstar[i],EPS)  
        }
      }
  }
  #---------------- check
  # XcatInit <- disj(strataInit)
  # 
  # print(t(XInit/pikInit)%*%pikInit)
  # print(t(X/pik)%*%pikstar)
  # print(t(XInit/pikInit)%*%s_01)
  # 
  # 
  # print(t(Xcat)%*%pik)
  # print(t(Xcat)%*%pikstar)
  # print(t(XcatInit)%*%s_01)
  # 
  # 
  # print(sum(pikstar))
  # print(sum(s_01))
  
  if(rand == TRUE){
    s_01 <- rep(0,length(pikstar))
    if(landing == TRUE){
      # put to 1 if > than 0
      s_01[o_out[which(pikstar > (1-EPS))]] <- 1
    }else{
      s_01[o_out[which(pikstar > EPS)]] <- pikstar[which(pikstar > EPS)]
    }
  }else{
    s_01 <- round(pikstar,10)
  }
  
  
  # if(rand == TRUE){
  #   out <- rep(0,length(pikstar))
  #   # out[o_out[which(pikstar > (1-EPS))]] <- 1
  #   out[o_out[which(pikstar > EPS)]] <- pikstar[which(pikstar > EPS)]
  # }else{
  #   out <- round(pikstar,10)
  # }
  # return(out)
  
  return(s_01)
  
}


#' 
#' #' @noRd
#' stratifiedcube_wo_landing <- function(X,strata,pik,EPS = 1e-7,rand = TRUE,landing = TRUE){
#'   
#'   if(rand == TRUE){
#'     strataInit <- strata
#'     ## cleanstrata
#'     a = sort(unique(strata))
#'     b = 1:length(a)
#'     names(b) = a
#'     strata <- as.vector(b[as.character(strata)])
#'     
#'     ## RANDOMIZATION 
#'     o <- order(strata)
#'     # strata[o] # this order the vector
#'     
#'     o_split <- split(o, f = strata[o] ) # split and randomize in each category
#'     for( i in 1:length(o_split)){
#'       o_tmp <- sample(1:length(o_split[[i]]))
#'       o_split[[i]] <- o_split[[i]][o_tmp]
#'     }
#'     o_out <- unlist(o_split[sample(1:length(o_split))]) # unlist and randomize each category
#'     
#'     XInit <- X
#'     pikInit <- pik  
#'     
#'     X <- X[o_out,]
#'     strata <- strata[o_out]
#'     pik <- pik[o_out]
#'   }
#'   
#'   ##----------------------------------------------------------------
#'   ##                        Initialization                         -
#'   ##----------------------------------------------------------------
#'   # strata <- as.matrix(strata)
#'   A = X/pik
#'   pikstar <- pik
#'   p = ncol(X)
#'   nstrata <- length(unique(strata))
#'   # TOT <- t(A)%*%pik
#'   # XX <- MASS::ginv(t(A) %*% A)
#'   
#'   ##----------------------------------------------------------------
#'   ##                  Flightphase on each strata                   -
#'   ##----------------------------------------------------------------
#'   
#'   for(k in 1:nstrata){
#'     pikstar[strata == k] <-ffphase(as.matrix(cbind(pikstar[strata == k],X[which(strata == k),])),
#'                                    pikstar[strata == k])
#'   }
#'   
#'   ###################### CHECK
#'   # t(X/pik)%*%pikstar
#'   # t(X/pik)%*%pik
#'   # Xcat <- disj(strata)
#'   # t(Xcat)%*%pik
#'   # t(Xcat)%*%pikstar
#'   
#'   ##----------------------------------------------------------------
#'   ##                Number of non 0-1 inclusion prob               -
#'   ##----------------------------------------------------------------
#'   
#'   i <- which(pikstar > EPS & pikstar < (1-EPS))
#'   i_size = length(i)
#'   i_size
#'   
#'   
#'   ##----------------------------------------------------------------
#'   ##            flight phase with categorical variable             -
#'   ##----------------------------------------------------------------
#'   
#'   if(i_size > 0 ){
#'     while(i_size > 0){
#'       
#'       ##------ Remove unique category
#'       uCat <- i[duplicated(strata[i]) | duplicated(strata[i], fromLast = TRUE)]
#'       if(length(uCat) == 0){
#'         break;
#'       }
#'       
#'       ##------ Find B
#'       A_tmp <- X[uCat,]/pik[uCat]
#'       B <- findB(A_tmp,strata[uCat])
#'       B <- cbind(B$X,B$Xcat)
#'       ##------ onestep and check if null
#'       tmp <-  onestep(B,pikstar[uCat[1:nrow(B)]],EPS)
#'       if(is.null(tmp)){
#'         break;
#'       }else{
#'         pikstar[uCat[1:nrow(B)]] <- tmp
#'       }
#'       
#'       ##------ update i
#'       
#'       i <- which(pikstar > EPS & pikstar < (1-EPS))
#'       i_size = length(i)
#'       # print(sum(pikstar))
#'     }
#'     
#'     
#'     # ##----------------------------------------------------------------
#'     # ##          end of flight phase on strata categories             -
#'     # ##----------------------------------------------------------------
#'     # Sometimes some stata does could not be balanced at the end and so we 
#'     # drop some auxiliary variable such that we could have only one unit that
#'     # are not put equal to 0 or 1 within each strata
#'     #
#'     
#'     # p <- ncol(X)
#'     # k = 1
#'     # while(length(uCat) != 0){
#'     #   ##------ Find B
#'     #   if(k == p){
#'     #     B <- as.matrix(pikstar[uCat])/pikstar[uCat]
#'     #   }else{
#'     #     A_tmp <- X[uCat,1:(p-k)]/pik[uCat]
#'     #     B <- findB(A_tmp,strata[uCat])
#'     #     B <- cbind(B$X,B$Xcat)
#'     #   }
#'     #   
#'     #   # print(pikstar[uCat[1:nrow(B)]])
#'     #   tmp <-  onestep(B,pikstar[uCat[1:nrow(B)]],EPS)
#'     #   # print(tmp)
#'     #   if(!is.null(tmp)){
#'     #     pikstar[uCat[1:nrow(B)]] <- tmp  
#'     #   }
#'     #   i <- which(pikstar > EPS & pikstar < (1-EPS))
#'     #   i_size = length(i)
#'     #   uCat <- i[duplicated(strata[i]) | duplicated(strata[i], fromLast = TRUE)]
#'     #   k = k + 1
#'     # }
#'   }
#'   
#'   #---------------- check
#'   # Xcat <- disjMatrix(as.matrix(strata))
#'   # 
#'   # print(t(X/pik)%*%pik)
#'   # print(t(X/pik)%*%pikstar)
#'   # print(t(Xcat)%*%pik)
#'   # print(t(Xcat)%*%pikstar)
#'   # print(length(i))
#'   # print(sum(pikstar))
#' 
#'     
#'   if(rand == TRUE){
#'     out <- rep(0,length(pikstar))
#'     # out[o_out[which(pikstar > (1-EPS))]] <- 1
#'     out[o_out[which(pikstar > EPS)]] <- pikstar[which(pikstar > EPS)]
#'   }else{
#'     out <- round(pikstar,10)
#'   }
#'   return(out)
#'   
#'   
#' }
