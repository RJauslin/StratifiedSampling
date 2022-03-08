
######### Package loading

library(sampling)
library(Spbsampling)
library(BalancedSampling)
library(WaveSampling)
devtools::load_all(".")


######### set seed and rm environment

rm(list = ls())
source("C:/Users/Raphael/switchdrive/StratifiedSampling/StratifiedSampling/inst/balseq/sim.R")
set.seed(1234)

######### Data generation loading

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


######### Auxiliary variables and variable of interest

p <- 5 
Xaux <- matrix(rgamma(N*p,1,1/5),ncol = p)
n <- 50
pik <- rep(n/N,N)
beta <- c(1,2,3,4,5)
y <- list(Xaux%*%beta + rnorm(N,0,20))



Xaux <- list(eq = Xaux)
pik <- list(eq = pik)
Xspread <- list(X_pp)
 


######### Distance matrix

D_ns <- as.matrix(proxy::dist(X_ns))
D_pp <- as.matrix(proxy::dist(X_pp))

######### Sampling functions 


f <- list(balseq,
          lcube,
          lpm1,
          srswor,
          pwd)

names(f) <- c("stream",
              "lcube",
              "lpm1",
              "srswor",
              "pwd")

######### Simulation


sim(1,f,Xaux,Xspread,D_pp,pik,y)



