


#' @noRd
split_pik <- function(pik){
  
  # initialize
  V0 <- 0
  cont <- 0
  
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

  ##############################
  ### If H is more than one ####
  ##############################
  if(H>1){
    
  
    ####  Indicating Cross-border #####
    abmat <- split_pik(pik)
    id<-c(abmat[1,])  
    aa<-c(abmat[2,])
    bb<-c(abmat[3,])
    
    
    ### Selecting the sample ####
    lab<-length(abmat[1,])
    pikkk<-c(pik[1:(id[1]-1)],aa[1])
    fsam<-osod(pikkk)
    a[1:(id[1]-1)]<-fsam[c(1:(length(fsam)-1))]
    for(k in 1:(lab-1)){
      pikkk<-c(pik[c((id[k]+1):(id[k+1]-1))],aa[k+1])
      
      #----- why ceiling here... possible to have bad behaviour with 3.0000001 that turns out to be transofmr to 4
      # pikkk2<-inclusionprobabilities(pikkk,ceiling(sum(pikkk)))
      pikkk2 <- inclusionprobabilities(pikkk,sum(pikkk))
  
      
      if(fsam[length(fsam)]==1){
        a[id[k]]<-1
        pikkk<-pikkk2
        fsam<-osod(pikkk)
        a[(id[k]+1):(id[k+1]-1)]<-fsam[c(1:(length(fsam)-1))]
      }else{
        bb[k]<-(bb[k])/(1-aa[k])
        pikkk<-(pikkk-pikkk2*aa[k])/(1-aa[k])
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
    a < -a[-N]
  } 
  
  return(a)
}

