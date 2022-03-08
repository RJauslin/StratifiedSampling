

rm(list = ls())
library(sampling)
N <- 200
n <- 50
# CSR
pp <- spatstat.random::rpoispp(N)
X_pp <- matrix(cbind(pp$x,pp$y),nrow = pp$n,ncol = 2)
N_pp <- nrow(X_pp)
while(N_pp != N){
  pp <- spatstat.random::rpoispp(N)
  X_pp <- matrix(cbind(pp$x,pp$y),nrow = pp$n,ncol = 2)
  N_pp <- nrow(X_pp)
}


swap <- function(s,i,j){
  tmp <- s[i]
  s[i] <- s[j]
  s[j] <- tmp
  return(s)
}

D <- as.matrix(proxy::dist(X_pp))

s <- srswor(n,N)

s0 <- which(s == 0)
s1 <- which(s == 1)
sj <- s0[sample(length(s0),1)]
si <- s1[sample(length(s1),1)]

s[si]
s[sj]

se <- swap(s,si,sj)



MD <-  function(D,s){
  N <- nrow(D)
  p <- ncol(D)
  out <- 0
  for(i in 1:N){
    if(s[i] == 1){
      for(j in 1:p){
        if(s[j] == 1){
          out <- out + D[i,j]  
        }
      }
    }
  }
  return(out)
}


MD(D,s)
MD(D,se)  
  

p <- function(s,se,D,beta){
  out <- min(1,(MD(D,se)/MD(D,s))^beta)
  return(out)
}

prob <- p(s,se,D,3)





