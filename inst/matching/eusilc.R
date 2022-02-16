# install.packages("laeken")
rm(list = ls())
library(laeken)
library(sampling)
library(parallel)
devtools::load_all("C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling")
source("C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/matching/simu.R")


data("eusilc")
# data("ses")
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




# categorical age

# c_age <- eusilc$age
# c_age[which(eusilc$age <= 20)] <- "(0,20]" 
# c_age[which(20 < eusilc$age & eusilc$age <= 40)] <- "(20,40]" 
# c_age[which(40 < eusilc$age & eusilc$age <= 60)] <- "(40,60]" 
# c_age[which(eusilc$age > 60)] <- "(60,Inf)"

# categorical income

c_income  <- eusilc$eqIncome
q <- quantile(eusilc$eqIncome, probs = seq(0, 1, 0.15))
c_income[which(eusilc$eqIncome <= q[2])] <- "(0,15]"
c_income[which(q[2] < eusilc$eqIncome & eusilc$eqIncome <= q[3])] <- "(15,30]"
c_income[which(q[3] < eusilc$eqIncome & eusilc$eqIncome <= q[4])] <- "(30,45]"
c_income[which(q[4] < eusilc$eqIncome & eusilc$eqIncome <= q[5])] <- "(45,60]"
c_income[which(q[5] < eusilc$eqIncome & eusilc$eqIncome <= q[6])] <- "(60,75]"
c_income[which(q[6] < eusilc$eqIncome & eusilc$eqIncome <= q[7])] <- "(75,90]"
c_income[which(  eusilc$eqIncome > q[7] )] <- "(90,100]"



# Y <- data.frame(c_age = c_age)
Y <- data.frame(ecostat = eusilc$pl030)
Z <- data.frame(c_income = c_income)

rownames(Xm) <- rownames(Y) <- rownames(Z)<- id

YZ <- table(cbind(Y,Z))
YZ
ncat1 <- nlevels(factor(Y[,1]))
ncat2 <- nlevels(factor(Z[,1]))



## RENSSEN

N <- nrow(eusilc)
n1 <- 600
n2 <- 3000



######################### SRSWOR
# s1 <- srswor(n1,N)
# s2 <- srswor(n2,N)



SIM <- 10000

##------------- simu
# res_ren <- res_opt <- res_ran <- list()
# l_save <- list()
# for(i in 1:SIM){
#   print(i)
#   # i = 3078
#   # set.seed(i+3077)
#   l <- simu(Xm,Y,Z,id,n1,n2)
#   # l_save[[i]] <- l
#   # saveRDS(l_save, file = "C:/Users/jauslinr/Desktop/l_save.rds")
# 
#   res_ren[[i]] <- l$YZ_ren
#   res_opt[[i]] <- l$YZ_opt
#   res_ran[[i]] <- l$YZ_ran
# 
# }
# 
# 
# # l <- readRDS("C:/Users/jauslinr/Desktop/l_save10000v2.rds")
# # l <- readRDS(file = "C:/Users/jauslinr/switchdrive/matching_optimal_transport/l_save10000.rds")
# 
# res_ren <- lapply(l,FUN = function(x){return(x$YZ_ren)})
# res_opt <- lapply(l,FUN = function(x){return(x$YZ_opt)})
# res_ran <- lapply(l,FUN = function(x){return(x$YZ_ran)})
# SIM <- length(l)
# # l <- readRDS("C:/Users/jauslinr/Desktop/l_save10000v2.rds")
# res_ren <- lapply(l,FUN = function(x){return(x$YZ_ren)})
# res_opt <- lapply(l,FUN = function(x){return(x$YZ_opt)})
# res_ran <- lapply(l,FUN = function(x){return(x$YZ_ran)})
# SIM <- length(l)

##------------- simu parrallel
cl <- makeCluster(detectCores())
clusterEvalQ(cl,{
  library(devtools)
  library(BalancedSampling)
  devtools::load_all("C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling")
  source("C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/matching/simu.R")

})

f1 <- function(n,Xm,Y,Z,id,n1,n2){
  l <- simu(Xm,Y,Z,id,n1,n2,totals = FALSE)
  return(list(YZ_ren =  l$YZ_ren,YZ_opt = l$YZ_opt, YZ_ran = l$YZ_ran))
}


# f1(1,Xm,Y,Z,id,n1,n2)


start <- Sys.time ()
l_sim <- parLapply(cl = cl,
                   X = 1:SIM,
                   fun = f1,
                   Xm = Xm,
                   Y = Y,
                   Z = Z,
                   id = id,
                   n1 = n1,
                   n2 = n2)
print(Sys.time () - start)

saveRDS(l_sim, file = "C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/matching/l_simeusilc.rds")

stopCluster(cl)


res_ren <- lapply(l_sim,function(x){x$YZ_ren})
res_opt <- lapply(l_sim,function(x){x$YZ_opt})
res_ran <- lapply(l_sim,function(x){x$YZ_ran})


do.call(rbind,lapply(res_ren,function(x){return(dim(x))} ))

##------------- Remove NA

res_opt <- lapply(res_opt,function(x){x[is.na(x)] <- 0; return(x)} )
res_ran <- lapply(res_ran,function(x){x[is.na(x)] <- 0; return(x)} )
res_ren <- lapply(res_ren,function(x){x[is.na(x)] <- 0; return(x)} )

##------------- MSE



MSE_ren <- lapply(res_ren,FUN = function(x){ (x - YZ)^2})
tmp <- MSE_ren[[1]]
for(j in 2:SIM){
  tmp <- tmp + MSE_ren[[j]]
}
MSE_ren <- tmp/SIM

MSE_opt <- lapply(res_opt,FUN = function(x){ (x - YZ)^2})
tmp <- MSE_opt[[1]]
for(j in 2:SIM){
  tmp <- tmp + MSE_opt[[j]]
}
MSE_opt <- tmp/SIM

MSE_ran <- lapply(res_ran,FUN = function(x){ (x - YZ)^2})
tmp <- MSE_ran[[1]]
for(j in 2:SIM){
  tmp <- tmp + MSE_ran[[j]]
}
MSE_ran <- tmp/SIM


##------------- Biais
b_ren <- lapply(res_ren,FUN = function(x){ (x - YZ)})
tmp <- b_ren[[1]]
for(j in 2:SIM){
  tmp <- tmp + b_ren[[j]]
}
b_ren <- tmp/SIM

b_opt <- lapply(res_opt,FUN = function(x){ (x - YZ)})
tmp <- b_opt[[1]]
for(j in 2:SIM){
  tmp <- tmp + b_opt[[j]]
}
b_opt <- tmp/SIM

b_ran <- lapply(res_ran,FUN = function(x){ (x - YZ)})
tmp <- b_ran[[1]]
for(j in 2:SIM){
  tmp <- tmp + b_ran[[j]]
}
b_ran <- tmp/SIM


##------------- var

m_ren <- res_ren[[1]]
for(j in 2:SIM){
  m_ren <- m_ren + res_ren[[j]]
}
m_ren <- m_ren/SIM

v_ren <- (res_ren[[1]]-m_ren)^2
for(j in 2:SIM){
  v_ren <- v_ren + (res_ren[[j]]-m_ren)^2
}
v_ren <- v_ren/SIM


m_opt <- res_opt[[1]]
for(j in 2:SIM){
  m_opt <- m_opt + res_opt[[j]]
}
m_opt <- m_opt/SIM

v_opt <- (res_opt[[1]]-m_opt)^2
for(j in 2:SIM){
  v_opt <- v_opt + (res_opt[[j]]-m_opt)^2
}
v_opt <- v_opt/SIM




m_ran <- res_ran[[1]]
for(j in 2:SIM){
  m_ran <- m_ran + res_ran[[j]]
}
m_ran <- m_ran/SIM

v_ran <- (res_ran[[1]]-m_ran)^2
for(j in 2:SIM){
  v_ran <- v_ran + (res_ran[[j]]-m_ran)^2
}
v_ran <- v_ran/SIM



saveRDS(res_opt, file = "C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/matching/res_opt.rds")
saveRDS(res_ren, file = "C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/matching/res_ren.rds")
saveRDS(res_ran, file = "C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/matching/res_ran.rds")


saveRDS(MSE_opt, file = "C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/matching/MSE_opt.rds")
saveRDS(MSE_ren, file = "C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/matching/MSE_ren.rds")
saveRDS(MSE_ran, file = "C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/matching/MSE_ran.rds")

saveRDS(b_opt, file = "C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/matching/b_opt.rds")
saveRDS(b_ren, file = "C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/matching/b_ren.rds")
saveRDS(b_ran, file = "C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/matching/b_ran.rds")

saveRDS(v_opt, file = "C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/matching/v_opt.rds")
saveRDS(v_ren, file = "C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/matching/v_ren.rds")
saveRDS(v_ran, file = "C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/matching/v_ran.rds")


##------------- khi2



phi <- function(YZ){
  E <- rowSums(YZ)%*%t(colSums(YZ))/sum(YZ)
  k <- sum((YZ - E)^2/E)
  return(sqrt(k/sum(YZ)))
}


phi_true <- phi(YZ)
phi_opt <- mean(do.call(rbind,lapply(res_opt,FUN = phi)),na.rm = TRUE)
phi_ren <- mean(do.call(rbind,lapply(res_ren,FUN = phi)),na.rm = TRUE)
phi_ran <- mean(do.call(rbind,lapply(res_ran,FUN = phi)),na.rm = TRUE)


MSE_phi_opt <- mean((do.call(rbind,lapply(res_opt,FUN = phi)) - phi_true)^2,na.rm = TRUE)

MSE_phi_ren <-  mean((do.call(rbind,lapply(res_ren,FUN = phi)) - phi_true)^2,na.rm = TRUE)

MSE_phi_ran <- mean((do.call(rbind,lapply(res_ran,FUN = phi)) - phi_true)^2,na.rm = TRUE)

b_phi_opt <- mean((do.call(rbind,lapply(res_opt,FUN = phi)) - phi_true),na.rm = TRUE)
b_phi_ren <- mean((do.call(rbind,lapply(res_ren,FUN = phi)) - phi_true),na.rm = TRUE)
b_phi_ran <- mean((do.call(rbind,lapply(res_ran,FUN = phi)) - phi_true),na.rm = TRUE)


var_phi_opt <- mean((do.call(rbind,lapply(res_opt,FUN = phi)) - phi_opt)^2,na.rm = TRUE)
var_phi_ren <- mean((do.call(rbind,lapply(res_ren,FUN = phi)) - phi_ren)^2,na.rm = TRUE)
var_phi_ran <- mean((do.call(rbind,lapply(res_ran,FUN = phi)) - phi_ran)^2,na.rm = TRUE)


phi_tab <- data.frame(phi = c(phi_true,phi_ren,phi_opt,phi_ran),
                      mse = c(NA,MSE_phi_ren,MSE_phi_opt,MSE_phi_ran),
                      b = c(NA,b_phi_ren,b_phi_opt,b_phi_ran) ,
                      v = c(NA,var_phi_ren,var_phi_opt,var_phi_ran))


phi_tab <- round(phi_tab,5)
phi_tab[1,2:4] <- "-"
phi_tab
saveRDS(phi_tab, file = "C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/matching/phi_tab2.rds")

##------------- results

MSE_ren
v_ren + b_ren^2
v_ren
round(b_ren^2,4)



MSE_opt
v_opt + b_opt^2
round(b_opt^2,4)
v_opt



MSE_ran
v_ran + b_ran^2
v_ran
round(b_ran^2,4)




sum(MSE_ren)
sum(MSE_opt)
sum(MSE_ran)

sum(b_ren^2)
sum(b_opt^2)
sum(b_ran^2)

sum(v_ren)
sum(v_opt)
sum(v_ran)




