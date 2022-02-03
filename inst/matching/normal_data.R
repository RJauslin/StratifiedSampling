

rm(list = ls())
library(MASS)
library(BalancedSampling)

Posdef <- function (n, ev = runif(n, 0, 10)) 
{
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

# test
set.seed(1)
Sigma <- Posdef(3)
N <- 10000
Xm <- mvrnorm(N,mu = c(0,0,0),Sigma)



# b1 <- matrix(c(2,-3,1),ncol = 1)
b1 <- matrix(c(2,-3,1,1.2,4,-5),ncol = 2)
Y <- Xm%*%b1


# b2 <- matrix(c(-4,1,-3),ncol = 1)
b2 <- matrix(c(-4,1,-3,-1.4,2.3,-5.6),ncol = 2)
Z <- Xm%*%b2


ALL <- cbind(Xm,Y,Z)
S <- cov(ALL)


# S_YX%*%(S_XX^{-1})%*%S_XZ -------------> CIA
# S_YX <- S[4,1:3]
# S_XZ <- S[1:3,5]

S_YX <- S[4:5,1:3]
S_XZ <- S[1:3,6:7]

S_YX%*%solve(S[1:3,1:3])%*%S_XZ
S
cov(Y,Z)

#------------------------------------ CIA HOLDS 

n1 <- 600
n2 <- 3000
id <- seq(1,N,1)
rownames(Xm) <- id

# simu_c <- function(Xm,Y,Z,id,n1,n2,totals = FALSE){
# 
#   # SET UP
#   y <- colnames(Y)
#   z <- colnames(Z)
#   N <- nrow(Xm)
#   
#   # samples
#   pik1 <- rep(n1/N,N)
#   pik2 <- rep(n2/N,N)
#   
#   Xm1 <- cbind(pik1,Xm)
#   Xm2 <- cbind(pik2,Xm)
# 
#   s1 <- rep(0,N)
#   s2 <- rep(0,N)
#   
#   s1[cube(pik1,as.matrix(Xm1))] <- 1
#   s2[cube(pik2,as.matrix(Xm2))] <- 1
#   
#   
#   X1 <- Xm[s1 == 1,]
#   X2 <- Xm[s2 == 1,]
#   
#   
#   # Y1 observed Y for S1
#   # Z2 observed Z for S2
#   
#   Y1 <- data.frame(Y[s1 == 1,])
#   colnames(Y1) <- y
#   Z2 <- data.frame(Z[s2 == 1,])
#   colnames(Z2) <- z
#   
#   id1 <- id[s1 == 1]
#   id2 <- id[s2 == 1]
#   
#   rownames(Y1) <- id1
#   rownames(Z2) <- id2
#   
#   d1 <- rep(N/n1,n1)
#   d2 <- rep(N/n2,n2)
# 
#   # harmonization
#   if(totals == TRUE){
#     re <- harmonize(X1,d1,id1,X2,d2,id2,totals = c(N,colSums(Xm)))
#   }else{
#     re <- harmonize(X1,d1,id1,X2,d2,id2)
#   }
#   w1 = re$w1
#   w2 = re$w2
#   
#   cat("Harmonization done \n\n")
#   if(abs(sum(w1) - sum(w2)) > 1e-7){
#     w2 <- w2/sum(w2)*sum(w1)
#   }
#   
#   
#   
#   #----------------- RENSSEN
#   
# 
#   fit1 <- lm(as.matrix(Y1)~X1)
#   fit2 <- lm(as.matrix(Z2)~X2)
#   
#   
#   # predict Z for S1
#   Y_pred_S1 <- fit1$fitted.values
#   Z_pred_S1 <- cbind(rep(1,nrow(X1)),X1)%*%fit2$coefficients
#   
#   # predict Y for S2
#   Z_pred_S2 <- fit2$fitted.values
#   Y_pred_S2 <- cbind(rep(1,nrow(X2)),X2)%*%fit1$coefficients
#   
#   
#   m1_ren_S1 <- colSums(w1*Y_pred_S1)/sum(w1)
#   m2_ren_S1 <- colSums(w1*Z_pred_S1)/sum(w1)
#   
#   m1_ren_S2 <- colSums(w2*Y_pred_S2)/sum(w2)
#   m2_ren_S2 <- colSums(w2*Z_pred_S2)/sum(w2)
#   
#   
#   c_ren1 <- t(w1*(Y_pred_S1 - m1_ren_S1))%*%(Z_pred_S1 - m2_ren_S1)/sum(w1)
#   c_ren2 <- t(w2*(Y_pred_S2 - m1_ren_S2))%*%(Z_pred_S2 - m2_ren_S2)/sum(w2)
#   
#   
#   #------------------- OPTIMAL
#   object = otmatch(X1,id1,X2,id2,w1,w2)
#   
#   w_opt <- object[,3]
#   
#   y1_opt <- as.matrix(Y1[as.character(object$id1),])
#   z2_opt <-  as.matrix(Z2[as.character(object$id2),])
#   
#   m1_opt <- colSums(w_opt*y1_opt)/sum(w_opt)
#   m2_opt <- colSums(w_opt*z2_opt)/sum(w_opt)
#   
#   c_opt <- t(w_opt*(y1_opt - m1_opt))%*%(z2_opt - m2_opt)/sum(w_opt)
#   
#   
#   cat("Optimal done \n\n")  
#   #------------------- RANDOM
#   
#   out <- bsmatch(object)
#   
#   w_ran <- out$object[,3]
#   y1_ran <- as.matrix(Y1[as.character(out$object$id1),])
#   z2_ran <-  as.matrix(Z2[as.character(out$object$id2),])
# 
#   
#   m1_ran <- colSums(w_ran*y1_ran)/sum(w_ran)
#   m2_ran <- colSums(w_ran*z2_ran)/sum(w_ran)
#   
#   c_ran <- t(w_ran*(y1_ran - m1_ran))%*%(z2_ran - m2_ran)/sum(w_ran)
#   
#   
#   cat("Random done \n\n")
#   
#   
#   return(list(c_ren1 = c_ren1,
#               c_ren2 = c_ren2,
#               c_opt = c_opt,
#               c_ran = c_ran))
#   
# }



source("C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/matching/simu_c.R")

simu_c(Xm,Y,Z,id,n1,n2,totals = TRUE)
S[4:5,6:7]



SIM <- 100
library(parallel)

##------------- simu parrallel
cl <- makeCluster(detectCores())
clusterEvalQ(cl,{
  library(devtools)
  library(BalancedSampling)
  devtools::load_all("C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling")
  source("C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/matching/simu_c.R")
  
})

f <- function(n,Xm,Y,Z,id,n1,n2){
  l <- simu_c(Xm,Y,Z,id,n1,n2,totals = FALSE)
  return(l)
}


start <- Sys.time()
l_sim <- parLapply(cl = cl,
                   X = 1:SIM,
                   fun = f,
                   Xm = Xm,
                   Y = Y,
                   Z = Z,
                   id = id,
                   n1 = n1,
                   n2 = n2)

Sys.time() - start



c_opt <- lapply(l_sim,function(x){return(x$c_opt)})
c_ren1 <- lapply(l_sim,function(x){return(x$c_ren1)})
c_ren2 <- lapply(l_sim,function(x){return(x$c_ren2)})
c_ran <- lapply(l_sim,function(x){return(x$c_ran)})



cov_opt <- c_opt[[1]]/SIM
cov_ren1 <- c_ren1[[1]]/SIM
cov_ren2 <- c_ren2[[1]]/SIM
cov_ran <- c_ran[[1]]/SIM

for(i in 2:SIM){
  cov_opt <- cov_opt + c_opt[[i]]/SIM
  cov_ren1 <- cov_ren1 + c_ren1[[i]]/SIM
  cov_ren2 <- cov_ren2 + c_ren2[[i]]/SIM
  cov_ran <- cov_ran + c_ran[[i]]/SIM
}


cov_opt
cov_ren1
cov_ren2
cov_ran
S[4:5,6:7]


mse_opt <- (c_opt[[1]] - S[4:5,6:7])^2/SIM
mse_ren1 <- (c_ren1[[1]] - S[4:5,6:7])^2/SIM
mse_ren2 <- (c_ren2[[1]] - S[4:5,6:7])^2/SIM
mse_ran <- (c_ran[[1]] - S[4:5,6:7])^2/SIM

for(i in 2:SIM){
  mse_opt <- mse_opt + (c_opt[[i]] - S[4:5,6:7])^2/SIM
  mse_ren1 <- mse_ren1 + (c_ren1[[i]] - S[4:5,6:7])^2/SIM
  mse_ren2 <- mse_ren2 + (c_ren2[[i]] - S[4:5,6:7])^2/SIM
  mse_ran <- mse_ran + (c_ran[[i]] - S[4:5,6:7])^2/SIM
}



# 
# Xm <- X
# 
# n1 <- 600
# n2 <- 3000
# id <- seq(1,N,1)
# rownames(Xm) <- id
# # simu(Xm,Y,Z,id,n1,n2,totals = FALSE)
# 
# 
# # SET UP
# y <- colnames(Y)
# z <- colnames(Z)
# N <- nrow(Xm)
# 
# 
# # samples
# pik1 <- rep(n1/N,N)
# pik2 <- rep(n2/N,N)
# 
# Xm1 <- cbind(pik1,Xm)
# Xm2 <- cbind(pik2,Xm)
# # s1 <- landingRM(Xm1/pik1,ffphase(Xm1,pik1))
# # s2 <- landingRM(Xm2/pik2,ffphase(Xm2,pik2))
# 
# s1 <- rep(0,N)
# s2 <- rep(0,N)
# 
# s1[cube(pik1,as.matrix(Xm1))] <- 1
# s2[cube(pik2,as.matrix(Xm2))] <- 1
# 
# # print(sum(s1))
# # print(sum(s2))
# 
# X1 <- Xm[s1 == 1,]
# X2 <- Xm[s2 == 1,]
# 
# 
# Y1 <- data.frame(Y[s1 == 1,])
# colnames(Y1) <- y
# Z2 <- data.frame(Z[s2 == 1,])
# colnames(Z2) <- z
# 
# id1 <- id[s1 == 1]
# id2 <- id[s2 == 1]
# 
# rownames(Y1) <- id1
# rownames(Z2) <- id2
# 
# d1 <- rep(N/n1,n1)
# d2 <- rep(N/n2,n2)
# 
# # disjunctive form
# # Y_dis <- sampling::disjunctive(as.matrix(Y))
# # Z_dis <- sampling::disjunctive(as.matrix(Z))
# 
# # Y1_dis <- sampling::disjunctive(as.matrix(Y1))
# # Z2_dis <- sampling::disjunctive(as.matrix(Z2))
# 
# # Y1_dis <- Y_dis[s1 ==1,]
# # Z2_dis <- Z_dis[s2 ==1,]
# 
# # harmonization
# # if(totals == TRUE){
#   re <- harmonize(X1,d1,id1,X2,d2,id2,totals = c(N,colSums(Xm)))
# # }else{
#   # re <- harmonize(X1,d1,id1,X2,d2,id2)  
# # }
# w1 = re$w1
# w2 = re$w2
# 
# cat("Harmonization done \n\n")
# if(abs(sum(w1) - sum(w2)) > 1e-7){
#   w2 <- w2/sum(w2)*sum(w1)
# }
# 
# 
# 
# 
# #----------------- RENSSEN
# 
# 
# 
# 
# fit1 <- lm(as.matrix(Y1)~X1)
# fit2 <- lm(as.matrix(Z2)~X2)
# 
# 
# # predict Z for S1
# Y_pred_S1 <- fit1$fitted.values
# Z_pred_S1 <- cbind(rep(1,nrow(X1)),X1)%*%fit2$coefficients
# 
# # predict Y for S2
# Z_pred_S2 <- fit2$fitted.values
# Y_pred_S2 <- cbind(rep(1,nrow(X2)),X2)%*%fit1$coefficients
# 
# 
# m1_ren_S1 <- colSums(w1*Y_pred_S1)/sum(w1)
# m2_ren_S1 <- colSums(w1*Z_pred_S1)/sum(w1)
# 
# m1_ren_S2 <- colSums(w2*Y_pred_S2)/sum(w2)
# m2_ren_S2 <- colSums(w2*Z_pred_S2)/sum(w2)
# 
# 
# c_ren1 <- t(w1*(Y_pred_S1 - m1_ren_S1))%*%(Z_pred_S1 - m2_ren_S1)/sum(w1)
# c_ren2 <- t(w2*(Y_pred_S2 - m1_ren_S2))%*%(Z_pred_S2 - m2_ren_S2)/sum(w2)
# # c_ren
# c_ren1
# c_ren2
# 
# 
# 
# #------------------- OPTIMAL
# object = otmatch(X1,id1,X2,id2,w1,w2)
# 
# Y1_optimal <- cbind(X1[as.character(object$id1),],y = Y1[as.character(object$id1),])
# Z2_optimal <- cbind(X2[as.character(object$id2),],z = Z2[as.character(object$id2),])
# 
# w <- object[,3]
# y1_opt <- as.matrix(Y1[as.character(object$id1),])
# z2_opt <-  as.matrix(Z2[as.character(object$id2),])
# 
# length(z2_opt)
# length(y1_opt)
# 
# 
# m1_opt <- colSums(w*y1_opt)/sum(w)
# m2_opt <- colSums(w*z2_opt)/sum(w)
# 
# 
# 
# 
# c_opt <- t(w*(y1_opt - m1_opt))%*%(z2_opt - m2_opt)/sum(w)
# 
# # sum((y1_opt - m1_opt)*(z2_opt - m2_opt)*w)/sum(w)
# S[4:5,6:7]
# # S[4,5]
# 
# cat("Optimal done \n\n")  
# #------------------- RANDOM
# 
# out <- bsmatch(object)
# 
# w <- out$object[,3]
# y1_opt <- as.matrix(Y1[as.character(out$object$id1),])
# z2_opt <-  as.matrix(Z2[as.character(out$object$id2),])
# 
# length(z2_opt)
# length(y1_opt)
# 
# 
# m1_opt <- colSums(w*y1_opt)/sum(w)
# m2_opt <- colSums(w*z2_opt)/sum(w)
# 
# 
# 
# 
# c_ran <- t(w*(y1_opt - m1_opt))%*%(z2_opt - m2_opt)/sum(w)
# 
# # sum((y1_opt - m1_opt)*(z2_opt - m2_opt)*w)/sum(w)
# S[4:5,6:7]
# # S[4,5]
# 
# 
# cat("Random done \n\n")
# 
# 
# 
# c_ren1
# c_ren2 
# c_opt
# c_ran
# S[4:5,6:7]
# 
