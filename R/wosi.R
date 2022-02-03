


#' @noRd
split_pik <- function(pik,H){
  
  
  # initialize
  V0 <- 0
  cont <- 0
  N <- length(pik)
  
  # split
  n_int <- round(sum(pik)/H) # number of interval
  abmat <- matrix(0,ncol = n_int - 1,nrow = 3)
  
  
  for(j in 1:(N-1)){
    Vm10 <- V0
    V0 <- V0 + pik[j]
    if((ceiling(V0)-floor(Vm10)) == 2 && (round(V0)%%H)==0){
      cont<-cont+1
      abmat[2,cont]<-ceiling(Vm10)-Vm10
      abmat[3,cont]<-V0-floor(V0)
      abmat[1,cont]<-j
    }    
    
  }
  
  abmat<-cbind(abmat,c(N,pik[N],0))
  abmat<-round(abmat,8)
  return(abmat)
}



#' Windows of size Integer
#'
#' @param pik A vector of inclusion probabilities.
#' @param H An integer, the considered windows size.
#'
#' @return A vector with elements equal to 0 or 1. The value 1 indicates that the unit is selected while the value 0 is for rejected units.
#' @export
#'
#' @examples
#' pik<-c(.4,.7,.1,.2,.4,.2,.2,.2,.1,.6,.2,.1,.4,.5,.5,.6,.1,.2)
#' H <- 2
#' s <- wosi(pik,2)
#' 
wosi <- function(pik,H)
{
  ### Phantom unit ###
  phan<-0
  if(sum(pik)%%1!=0){
    pik<-c(pik,ceiling(sum(pik))-sum(pik))
    phan<-1
  }
  
  N<-length(pik)
  a <- rep(0,N)
  
  
  
  # Problem general, the function is not exactly 
  # in a stream here, the implementation is splitting the full vector first.
  # We should implement a version such that you never apply a function on the full vector pik.
  
  
  
  ##############################
  ### If H is more than one ####
  ##############################
  if(H>1){
    
    
    ####  Indicating Cross-border #####
    
    
    abmat <- split_pik(pik,H) 
    id<-c(abmat[1,])  
    aa<-c(abmat[2,])
    bb<-c(abmat[3,])
    
    
    ### Selecting the sample ####
    lab<-length(abmat[1,])
    pikkk<-c(pik[1:(id[1]-1)],aa[1])
    fsam<-osod(pikkk)
    a[1:(id[1]-1)]<-fsam[c(1:(length(fsam)-1))]
    for(k in 1:(lab-1)){
      pikkk <- c(pik[c((id[k]+1):(id[k+1]-1))],aa[k+1])
      print(sum(pikkk))
      
      
      #----- why ceiling here... possible to have bad behaviour with 3.0000001 that turns out to be transofmr to 4
      pikkk2 <- inclusionprobabilities(pikkk,H)
      # pikkk2 <- inclusionprobabilities(pikkk,sum(pikkk))
      
      
      if(fsam[length(fsam)]==1){
        a[id[k]]<-1
        pikkk<-pikkk2
        print(sum(pikkk))
        fsam<-osod(pikkk)
        a[(id[k]+1):(id[k+1]-1)]<-fsam[c(1:(length(fsam)-1))]
      }else{
        
        print(aa[k])
        print(bb[k])
        print(bb[k] + aa[k])
        print(pik[id[k]])
        cat("\n\n")

        
        bb[k]<-(bb[k])/(1-aa[k])
        
        
        print(sum(pikkk))
        
        pikkk<-(pikkk - pikkk2*aa[k])/(1-aa[k])
        pikkk<-c(bb[k],pikkk)
        
        
        
        fsam<-osod(pikkk)
        a[(id[k]):(id[k+1]-1)]<-fsam[c(1:(length(fsam)-1))]
      }
    } 
    
  }
  
  ####################
  ### For Size One ###
  ####################
  if(H==1){
    V<-0
    for(k in 1:(N-1))
    {
      Vm1=V
      V=V+pik[k]
      
      flo<-ceiling(V)-floor(Vm1)-1
      clo<-ceiling(V)-sum(a)
      
      if(  runif(1) < (min(clo,1))* (pik[k]-(flo*(2-clo)*(1-Vm1%%1)))/
           ((1-flo)*(1-Vm1%%1)+(flo*(Vm1%%1*(2-clo)+pik[k]*(clo-1)))) ) a[k]=1
    }
    
  }
  
  
  a[N] = sum(pik)-sum(a)
  if(phan==1){
    a <- a[-N]
  } 
  
  return(a)
}















#' wosi
#'
#' @param pik 
#' @param H 
#'
#' @return
#' @export
#'
#' @examples
#' rm(list = ls())
#' 
#' library(sampling)
#' pik<-c(.4,.7,.1,.2,.4,.2,.2,.2,.1,.6,.2,.1,.4,.5,.5,.6,.1,.2)
#' pik <- inclusionprobabilities(runif(100),10)
#' H <- 2
#' s <- wosi2(pik,H)
#' s
#' sum(s)
wosi2 <- function(pik,H){
  
  
  ### Phantom unit ###
  phan<-0
  if(sum(pik)%%1!=0){
    pik<-c(pik,ceiling(sum(pik))-sum(pik))
    phan<-1
  }
  
  N <- length(pik)
  out <- rep(0,N)
  
  
  # FIRST STEP
  
  k <- 1
  m <- pik[k]
  while(m < H){
    k = k + 1
    m = m + pik[k]
  }
  
  pika <- H - (m-pik[k])
  pikb <- pik[k] - pika
  s_tmp <- osod(c(pik[1:(k-1)],pika)) # osod on [1:k-1 ; pika]
  out[1:(k-1)] <- s_tmp[1:(length(s_tmp)-1)] # update out only on [1:k]
  
  
  # NEXT STEPS
  m <- pikb
  for(i in (k+1):N){
    
    m = m + pik[i]
    
    if(m >= H){
      
      pikai <- H - (m - pik[i])
      pikstar <- c(pik[(k+1):(i-1)],pikai)
      pikstar2 <- inclusionprobabilities(pikstar,H)
      
      
      if(s_tmp[length(s_tmp)] == 1){
        
        out[k] <- 1
        s_tmp <- osod(pikstar2)  
        out[(k+1):(i-1)] <- s_tmp[1:(length(s_tmp)-1)]
        
      }else{
        
      
        pikb <- pikb/(1-pika)
        pikstar2 <- (pikstar - pikstar2*pika)/(1-pika)
        pikstar2 <- c(pikb,pikstar2)
        s_tmp <- osod(pikstar2)  
        out[k:(i-1)] <- s_tmp[1:(length(s_tmp)-1)]
        
      }
      
      pika <- pikai
      pikb <- pik[i] - pika # update pikb
      m <- pikb
      k = i
      
    }
    
  }
  
  out[N] = sum(pik)-sum(out) # to take the decision for the last unit
  if(phan==1){
    out <- out[-N]
  } 
  
  return(out)
}









