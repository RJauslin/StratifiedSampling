

rm(list = ls())

library(sampling)

N <- 1000

Xcat <-as.matrix(data.frame(cat1 = rep(1:40,each = N/40),
                    cat2 = rep(1:50,each = N/50),
                    # cat3 = rep(1:500,times = 2)))
                    cat3 = rep(1:100,each = N/100)))

pik <- inclusionprobastrata(Xcat[,1],rep(1,40))
p <- 30
X <- matrix(rnorm(N*p),ncol = 30)
Xdev <- disjMatrix(Xcat)

Xall <- cbind(X,Xdev)
J <- ncol(Xall)
B1 <- Xall[1:(J+1),]
B1 <- ReducedMatrix(B1)
dim(B1$B)

B2 <- findB(X,Xcat)
dim(B2)

image(as(B1$B,"sparseMatrix"))
image(as(B2,"sparseMatrix"))
