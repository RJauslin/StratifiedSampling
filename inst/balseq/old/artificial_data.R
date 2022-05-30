
library(sampling)
library(Spbsampling)
rm(list = ls())
N <- 300
set.seed(1234)
# CSR
pp <- spatstat.random::rpoispp(N)
X_pp <- matrix(cbind(pp$x,pp$y),nrow = pp$n,ncol = 2)
N_pp <- nrow(X_pp)
while(N_pp != N){
  pp <- spatstat.random::rpoispp(N)
  X_pp <- matrix(cbind(pp$x,pp$y),nrow = pp$n,ncol = 2)
  N_pp <- nrow(X_pp)
}

# Neyman-Scottt
# nclust <- function(x0, y0, radius, n) {
# return(spatstat.random::runifdisc(n, radius, centre=c(x0, y0)))
# }
# N_ns <- 0
# while(N_ns != N){
#   ns <- spatstat.random::rNeymanScott(17, 1, nclust, radius=0.08, n=17)
#   X_ns <- matrix(cbind(ns$x,ns$y),nrow = ns$n,ncol = 2)
#   N_ns <- nrow(X_ns)
#   print(N_ns)
# }
# saveRDS(X_ns,"C:/Users/Raphael/switchdrive/StratifiedSampling/StratifiedSampling/inst/balseq/X_ns.rds")

X_ns <- readRDS("C:/Users/Raphael/switchdrive/StratifiedSampling/StratifiedSampling/inst/balseq/X_ns.rds")

plot(X_ns)
p <- 5 
Xaux <- matrix(rgamma(N*p,1,1/5),ncol = p)

n <- 10
pik <- rep(n/N,N)

s_pp <- balseq(pik,as.matrix(pik),X_pp)

o <- sample(N,N)
pik <- pik[o]
X_ns <- X_ns[o,]
s_ns <- balseq(pik,as.matrix(pik),X_ns)

D_ns <- as.matrix(proxy::dist(X_ns))
D_pp <- as.matrix(proxy::dist(X_pp))

s_pwd_pp <- as.vector(pwd(D_pp,n)$s)
s_pwd_ns <- as.vector(pwd(D_ns,n)$s)


plot(X_pp)
lines(X_pp[s_pwd_pp,],pch = 16,type = "p")

plot(X_pp)
lines(X_pp[s_pwd_pp,],pch = 16,type = "p")
plot(X_ns)
lines(X_ns[s_pwd_ns,],pch = 16,type = "p")





