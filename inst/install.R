library(devtools)
install_github("Rjauslin/fastcube@master")
library(fastcube)

rm(list = ls())
# set.seed(1)
eps <- 1e-13
library(Matrix)
N <- 500
Pik <- matrix(c(sampling::inclusionprobabilities(runif(N),50),
                sampling::inclusionprobabilities(runif(N),100),
                sampling::inclusionprobabilities(runif(N),150)),ncol = 3)
X <- PM(Pik)$PM
pik <- PM(Pik)$P
image(as(X,"sparseMatrix"))
# rrefArma(X)
dim(X)
order = 2
EPS = 1e-11
A=X/pik
system.time(test <- flightphase_arma(X,pik,redux = TRUE))
system.time(test1 <- ReducedSamplecube(X,pik,t = 3))
t(A)%*%test1
t(A)%*%pik




rm(list = ls())
N <- 100000
strata <- sample(x = 1:3400, size = N, replace = TRUE)
system.time(M <- disj(strata,3400))
system.time(M <- model.matrix(~as.factor(strata)-1))
system.time(M <- sampling::disjunctive(strata))



#'
#' #' rm(list = ls())
#' data(swissmunicipalities)
#' swiss=swissmunicipalities
#' X=cbind(swiss$HApoly,
#'         swiss$Surfacesbois,
#'         swiss$P00BMTOT,
#'         swiss$P00BWTOT,
#'         swiss$POPTOT,
#'         swiss$Pop020,
#'         swiss$Pop2040,
#'         swiss$Pop4065,
#'         swiss$Pop65P,
#'         swiss$H00PTOT )
#' pik=inclusionprobabilities(swiss$POPTOT,400)
#'
#'
#' Xcat <-data.frame(cat1 = swiss$REG)
#'
#' system.time(s <- fastcube(X,Xcat,pik))
#' system.time(s2 <-balancedstratification(X,swiss$REG,pik,comment=FALSE))
#' as.character(swiss$Nom[s2==1])
#' t(X/pik)%*%s2
#'
#' t(X/pik)%*%s
#' t(X/pik)%*%pik


rm(list = ls())

library(sampling)
N <- 10000

Xcat <-as.matrix(data.frame(cat1 = rep(1:40,each = N/40),
                            cat2 = rep(1:50,each = N/50),
                            cat2 = rep(1:100,each = N/100)))

p <- 30
X <- matrix(rnorm(N*p),ncol = 30)


Xcat_tmp <- disjMatrix(Xcat)
Xcat_tmp <- do.call(cbind,apply(Xcat,MARGIN = 2,disjunctive))
Xred <- as.matrix(cbind(X,Xcat_tmp))

pik <- rep(300/N,N)
A <- Xred/pik
system.time(s1 <- fastcubeArma(X,Xcat,pik))
system.time(s1 <- fastcube(X,Xcat,pik))
as.vector(t(A)%*%as.vector(s1))
as.vector(t(A)%*%pik)








rm(list = ls())
# set.seed(1)
library(sampling)
N <- 1000

Xcat <-as.matrix(data.frame(cat1 = rep(1:40,each = N/40),
                            cat2 = rep(1:50,each = N/50),
                            cat2 = rep(1:100,each = N/100)))

p <- 30
X <- matrix(rnorm(N*p),ncol = 30)


Xcat_tmp <- disjMatrix(Xcat)
# Xcat_tmp <- do.call(cbind,apply(Xcat,MARGIN = 2,disjunctive))
Xred <- as.matrix(cbind(X,Xcat_tmp))
dim(Xred)
pik <- rep(300/N,N)
A <- Xred/pik
system.time(s1 <- fastcubeArma(X,Xcat,pik))
sum(s1)
as.vector(t(A)%*%as.vector(s1))
as.vector(t(A)%*%pik)
system.time(s1 <- fastcube(X,Xcat,pik))
as.vector(t(A)%*%as.vector(s1))
as.vector(t(A)%*%pik)
