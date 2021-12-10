

rm(list = ls())
library(sampling)


fT <- function(w,i,j){
  tmp <- w^i
  res <- w[j]^i
  return(sum(tmp) - res)
}


R <- function(w,k,j){

  T_tmp <- c()
  for(i in 1:k){
    T_tmp[i] <- fT(w,i,j)
  }
  
  R_tmp <- rep(0,k)
  for(i in 1:k){
    for(l in 1:i){
      beta <-  (1/i)*(-1)^(l+1)*fT(w,l,j)
      if((i-l) == 0){
        R_tmp[i] <- R_tmp[i] + beta
      }else{
        R_tmp[i] <- R_tmp[i] + beta * R_tmp[i-l]      
      }
      
    }
  }
  return(R_tmp)
}

# verification
# R(w,3,N-1)
# fT(w,1,N-1)
# 
# 0.5*fT(w,1,N-1)^2 - 0.5*fT(w,2,N-1) # second entry
# (1/6)*fT(w,1,N-1)^3 - (1/6)*fT(w,1,N-1)*fT(w,2,N-1) - (1/3)*fT(w,2,N-1)*fT(w,1,N-1) + (1/3)*fT(w,3,N-1) # thrid entry



updatew <- function(pik,n,N){
  w <- pik
  
  for(t in 1:100){
    for(i in 1:N){
      num <- R(w,n-1,N)
      num <- num[length(num)]
      
      den <- R(w,n-1,i)
      den <- den[length(den)]
      
      
      w[i] <- pik[i]*num/den
    }  
  }
  return(w)
}





pik <- c(0.1,0.4,0.7,0.8)
w <- updatew(pik,2,4)
w <- w - mean(w)



sum(pik)
sum(w)




pik <- inclusionprobabilities(seq(1,100,by = 1),10)

N <- length(pik)
n <- sum(pik)
w <- updatew(pik,n,N)
w <- w - mean(w)

q <- qfromw(w,n)
s <- sfromq(q)
  









