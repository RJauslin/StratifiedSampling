

rm(list = ls())
library(MASS)
library(BalancedSampling)
devtools::load_all(".")


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
b1 <- matrix(c(0.2,-0.3,1,1.2,0.4,-0.5),ncol = 2)
Y <- Xm%*%b1 + matrix(rnorm(2*N),ncol = 2)
# Y <- Xm^2%*%b1 + matrix(rnorm(2*N),ncol = 2)


# b2 <- matrix(c(-4,1,-3),ncol = 1)
b2 <- matrix(c(-0.4,1,-0.3,-1.4,0.3,-0.6),ncol = 2)
Z <- Xm%*%b2 + matrix(rnorm(2*N),ncol = 2)
# Z <- Xm^3%*%b2 + matrix(rnorm(2*N),ncol = 2)


ALL <- cbind(Xm,Y,Z)
S <- cov(ALL)
S

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



source("C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/matching/simu_c.R")

simu_c(Xm,Y,Z,id,n1,n2,totals = FALSE)
S[4:5,6:7]





#########################################################

y <- colnames(Y)
z <- colnames(Z)
N <- nrow(Xm)

s1 <- srswor(n1,N)
s2 <- srswor(n2,N)
id <- seq(1,N,1)

s1 <- srswor(n1,N)
comple <- which(s1 == 0)
s2 <- rep(0,N)
s2[comple[sample.int(length(comple),n2)]] <- 1

id1 <- id[s1 == 1]
id2 <- id[s2 == 1]
base::intersect(id1,id2)

X1 <- Xm[s1 == 1,]
X2 <- Xm[s2 == 1,]


# Y1 observed Y for S1
# Z2 observed Z for S2

Y1 <- data.frame(Y[s1 == 1,])
colnames(Y1) <- y
Z2 <- data.frame(Z[s2 == 1,])
colnames(Z2) <- z

id1 <- id[s1 == 1]
id2 <- id[s2 == 1]

rownames(Y1) <- id1
rownames(Z2) <- id2

d1 <- rep(N/n1,n1)
d2 <- rep(N/n2,n2)

# harmonization

re <- harmonize(X1,d1,id1,X2,d2,id2)

w1 = re$w1
w2 = re$w2

cat("Harmonization done \n\n")
if(abs(sum(w1) - sum(w2)) > 1e-7){
  w2 <- w2/sum(w2)*sum(w1)
}


object = otmatch(X1,id1,X2,id2,w1,w2,transport_method = "revsimplex")

w_opt <- object[,3]

y1_opt <- as.matrix(Y1[as.character(object$id1),])
z2_opt <-  as.matrix(Z2[as.character(object$id2),])

apply(y1_opt,MARGIN = 2, FUN = function(x){w_opt*x})[1:7,]


m1_opt <- colSums(w_opt*y1_opt)/sum(w_opt)
m2_opt <- colSums(w_opt*z2_opt)/sum(w_opt)

c_opt <- t(w_opt*(y1_opt - m1_opt))%*%(z2_opt - m2_opt)/sum(w_opt)

c_opt
S[4:5,6:7]





#########################################################

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
  l <- simu_c(Xm,Y,Z,id,n1,n2,totals = TRUE)
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



# saveRDS(l_sim, file = "C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/matching/l_sim.rds")
# l_sim <- readRDS(file = "C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/matching/l_sim_10000.rds")

c_opt <- lapply(l_sim,function(x){return(x$c_opt)})
c_ren1 <- lapply(l_sim,function(x){return(x$c_ren1)})
c_ren2 <- lapply(l_sim,function(x){return(x$c_ren2)})
c_ran <- lapply(l_sim,function(x){return(x$c_ran)})



cov_opt <- c_opt[[1]]
cov_ren1 <- c_ren1[[1]]
cov_ren2 <- c_ren2[[1]]
cov_ran <- c_ran[[1]]

for(i in 2:SIM){
  cov_opt <- cov_opt + c_opt[[i]]
  cov_ren1 <- cov_ren1 + c_ren1[[i]]
  cov_ren2 <- cov_ren2 + c_ren2[[i]]
  cov_ran <- cov_ran + c_ran[[i]]
  
}


cov_opt/SIM -S[4:5,6:7]
cov_ren1/SIM - S[4:5,6:7]
cov_ren2/SIM - S[4:5,6:7]
cov_ran/SIM - S[4:5,6:7]
# S[4:5,6:7]




b2_opt <- (cov_opt - S[4:5,6:7])^2
b2_ren1 <- (cov_ren1 - S[4:5,6:7])^2
b2_ren2 <- (cov_ren2 - S[4:5,6:7])^2
b2_ran <- (cov_ran - S[4:5,6:7])^2


v_opt <- (c_opt[[1]] -  cov_opt)^2/SIM
v_ren1 <- (c_ren1[[1]] - cov_ren1)^2/SIM
v_ren2 <- (c_ren2[[1]] - cov_ren2)^2/SIM
v_ran <- (c_ran[[1]] - cov_ran)^2/SIM

for(i in 2:SIM){
  v_opt <- v_opt + (c_opt[[i]] - cov_opt)^2/SIM
  v_ren1 <- v_ren1 + (c_ren1[[i]] - cov_ren1)^2/SIM
  v_ren2 <- v_ren2 + (c_ren2[[i]] - cov_ren2)^2/SIM
  v_ran <- v_ran + (c_ran[[i]] - cov_ran)^2/SIM
}

v_opt
v_ren1 
v_ren2
v_ran






mse_opt <- (c_opt[[1]] - S[4:5,6:7])^2
mse_ren1 <- (c_ren1[[1]] - S[4:5,6:7])^2
mse_ren2 <- (c_ren2[[1]] - S[4:5,6:7])^2
mse_ran <- (c_ran[[1]] - S[4:5,6:7])^2

for(i in 2:SIM){
  mse_opt <- mse_opt + (c_opt[[i]] - S[4:5,6:7])^2
  mse_ren1 <- mse_ren1 + (c_ren1[[i]] - S[4:5,6:7])^2
  mse_ren2 <- mse_ren2 + (c_ren2[[i]] - S[4:5,6:7])^2
  mse_ran <- mse_ran + (c_ran[[i]] - S[4:5,6:7])^2
}

mse_opt
# b2_opt + v_opt
mse_ren1 
# b2_ren1 + v_ren1
mse_ren2
# b2_ren2 + v_ren2
mse_ran
# b2_ran + v_ran



tab_gauss <- round(rbind(cbind(b2_opt,v_opt,mse_opt),
      cbind(b2_ren1,v_ren1,mse_ren1),
      cbind(b2_ren2,v_ren2,mse_ren2),
      cbind(b2_ran,v_ran,mse_ran)),3)

library(magrittr) 
library(kableExtra)
kable(tab_gauss, format = "latex",digits = 3, booktabs = T,linesep = "", caption = " lknfgldfng ",row.names = TRUE,escape = FALSE) %>%
  kable_styling(latex_options = c("hold_position")) %>% 
  add_header_above(c("Biais$^2$" = 2,"Variance" = 2, "MSE" = 2)) %>%
  group_rows("Optimal transport",1,2,escape = F,bold =  F,latex_gap_space = "1ex") %>%
  group_rows("Renssen $S_1$",3,4,escape = F,bold =  F,latex_gap_space = "1ex") %>%
  group_rows("Renssen $S_2$",5,6,escape = F,bold =  F,latex_gap_space = "1ex") %>%
  group_rows("Balanced imputation",7,8,escape = F,bold =  F,latex_gap_space = "1ex")




















