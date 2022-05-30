


######### Package loading

library(Spbsampling)

######### set seed and rm environment

rm(list = ls())

N <- 300

# CSR
pp <- spatstat.random::rpoispp(N)
X_pp <- matrix(cbind(pp$x,pp$y),nrow = pp$n,ncol = 2)
N_pp <- nrow(X_pp)
while(N_pp != N){
  pp <- spatstat.random::rpoispp(N)
  X_pp <- matrix(cbind(pp$x,pp$y),nrow = pp$n,ncol = 2)
  N_pp <- nrow(X_pp)
}

# distance
D <- as.matrix(proxy::dist(X_pp))
con <- rep(0,nrow(D))
stand_D <- stprod(mat = D, con = con)$mat
n <- 60

s <- pwd(stand_D,n)$s

plot(X_pp)
points(X_pp[s,],pch = 16)
