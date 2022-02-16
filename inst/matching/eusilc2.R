# install.packages("laeken")
rm(list = ls())
library(laeken)
library(sampling)
library(parallel)
devtools::load_all("C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling")
source("C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/matching/simu.R")

# k = 1
# set.seed(k)

data("eusilc")
# data("ses")
eusilc
eusilc <- na.omit(eusilc)

Xm <- eusilc[,c("hsize","db040","age","rb090","pb220a")]
Xmcat <- do.call(cbind,apply(Xm[,c(2,4,5)],MARGIN = 2,FUN = disjunctive))
Xm <- cbind(Xmcat,Xm[,-c(2,4,5)])
id <- eusilc$rb030

Y <- data.frame(ecostat = eusilc$pl030)
Z <- data.frame(income = eusilc$eqIncome)

rownames(Xm) <- rownames(Y) <- rownames(Z)<- id

YZ <- cbind(Y,Z)
ll <- split(YZ,f = YZ$ecostat)
ZbyY <- tapply(YZ$income,list(YZ$ecostat),mean)
tapply(object[,3]*Z2_optimal$z,list(Y1_optimal$y),mean)

library(ggplot2)
ggplot(YZ,aes(x = ecostat,y = income))+
  geom_boxplot()

ncat1 <- nlevels(factor(Y[,1]))



N <- nrow(eusilc)
n1 <- 600
n2 <- 3000


# SET UP
y <- colnames(Y)
z <- colnames(Z)
N <- nrow(Xm)
ZbyY <- t(do.call(rbind,lapply(ll,function(x){mean(x$income)})))  
  
# samples
pik1 <- rep(n1/N,N)
pik2 <- rep(n2/N,N)

# s1 <- stratifiedcube(as.matrix(Xm),strata = as.numeric(as.character(Y$ecostat)),pik1)
# s2 <- stratifiedcube(as.matrix(Xm),strata = as.numeric(as.character(Y$ecostat)),pik2)
# print(sum(s2))
# k = k + 1

s1 <- srswor(n1,N)
s2 <- srswor(n2,N)


X1 <- Xm[s1 == 1,]
X2 <- Xm[s2 == 1,]
  
Y1 <- data.frame(Y[s1 == 1,])
colnames(Y1) <- y
Z2 <- as.data.frame(Z[s2 == 1,])
colnames(Z2) <- z
  
id1 <- id[s1 == 1]
id2 <- id[s2 == 1]
  
rownames(Y1) <- id1
rownames(Z2) <- id2
  
d1 <- rep(N/n1,n1)
d2 <- rep(N/n2,n2)
  
# disjunctive form
Y_dis <- sampling::disjunctive(as.matrix(Y))
# Z_dis <- sampling::disjunctive(as.matrix(Z))
  
  # Y1_dis <- sampling::disjunctive(as.matrix(Y1))
  # Z2_dis <- sampling::disjunctive(as.matrix(Z2))
  
Y1_dis <- Y_dis[s1 ==1,]
# Z2_dis <- Z_dis[s2 ==1,]
  
re <- harmonize(X1,d1,id1,X2,d2,id2,totals = c(N,colSums(Xm)))
    
w1 = re$w1
w2 = re$w2
 
cat("Harmonization done \n\n")
if(abs(sum(w1) - sum(w2)) > 1e-7){
    w2 <- w2/sum(w2)*sum(w1)
}

#----------------- optimal transport

object <- otmatch(X1,id1,X2,id2,w1,w2)
# Y1_optimal <- cbind(X1[as.character(object$id1),],y = Y1[as.character(object$id1),])
# Z2_optimal <- cbind(X2[as.character(object$id2),],z = Z2[as.character(object$id2),])

res_opt <- tapply(object$weight*Z2[as.character(object$id2),],list(Y1[as.character(object$id1),]),sum)/tapply(object$weight,list(Y1[as.character(object$id1),]),sum)
ZbyY <- tapply(YZ$income,list(YZ$ecostat),mean)
ZbyY
res_opt
sum((res_opt - ZbyY)^2)

 #------------------- balanced sampling match

out <- bsmatch(object)

res_ran <- tapply( (out$object$weight/out$q)*Z2[as.character(out$object$id2),],list(Y1[as.character(out$object$id1),]),sum)/tapply(out$object$weight/out$q,list(Y1[as.character(out$object$id1),]),sum)
res_ran
sum((res_ran - ZbyY)^2)




#----------------- RENSSEN

QR.1 <- qr(X1 * sqrt(w1))
beta.yx.1 <- qr.coef(QR.1, Y1_dis * sqrt(w1))
beta.yx.1[is.na(beta.yx.1)] <- 0
QR.2 <- qr(X2 * sqrt(w2))
beta.zx.2 <- qr.coef(QR.2, Z2 * sqrt(w2))
beta.zx.2[is.na(beta.zx.2)] <- 0
XX.w1 <- t(as.matrix(X1)) %*% (as.matrix(X1) * w1)
XX.w2 <- t(as.matrix(X2)) %*% (as.matrix(X2) * w2)

n12=length(intersect(id1,id2))
gamma.p=(n1-n12)/(n1+n2-2*n12)  

# gamma.p <- n1/(n1 + n2)
XX.pool <- gamma.p * XX.w1 + (1 - gamma.p) * XX.w2
YZ.CIA <- t(beta.yx.1) %*% XX.pool %*% beta.zx.2


res_ren <- as.vector(YZ.CIA)/tapply(w1,list(Y1$ecostat),sum)
sum((res_ren - ZbyY)^2)



sum((res_opt - ZbyY)^2)*1e-6
sum((res_ran - ZbyY)^2)*1e-6
sum((res_ren - ZbyY)^2)*1e-6








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




