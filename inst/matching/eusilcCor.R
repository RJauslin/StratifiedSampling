# install.packages("laeken")
rm(list = ls())
library(laeken)
library(sampling)
library(parallel)
devtools::load_all("C:/Users/jauslinr/switchdrive/matching_optimal_transport/transportMatching")
source("C:/Users/jauslinr/switchdrive/matching_optimal_transport/simucode/simu.R")


data("eusilc")
data("ses")
eusilc
eusilc <- na.omit(eusilc)


# Neu <- nrow(eusilc)
# eusilc <- eusilc[sample(1:Neu,5000),]


# Xm <- eusilc[,-c(1,5,26,27,28,4)]
# Xmcat <- do.call(cbind,apply(Xm[,1:5],MARGIN = 2,FUN = disjunctive))
# Xm <- cbind(Xmcat,Xm[,-(1:5)])


Xm <- eusilc[,c("hsize","db040","age","rb090","pb220a")]
Xmcat <- do.call(cbind,apply(Xm[,c(2,4,5)],MARGIN = 2,FUN = disjunctive))
Xm <- cbind(Xmcat,Xm[,-c(2,4,5)])
id <- eusilc$rb030


# Y <- data.frame(c_age = c_age)
Y <- data.frame(hy040n = eusilc$hy080n)
Z <- data.frame(eqIncome = eusilc$eqIncome)

rownames(Xm) <- rownames(Y) <- rownames(Z)<- id

# YZ <- table(cbind(Y,Z))
c_YZ <- cor(Y,Z)
c_YZ





mean((Y$hy040n - mean(Y$hy040n))*(Z$eqIncome - mean(Z$eqIncome)))/ sqrt(mean((Y$hy040n - mean(Y$hy040n))^2)*mean((Z$eqIncome - mean(Z$eqIncome))^2))

## RENSSEN

N <- nrow(eusilc)
n1 <- 5000
n2 <- 3000



SIM <- 5
