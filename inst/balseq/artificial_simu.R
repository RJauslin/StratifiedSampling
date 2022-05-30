

######### Package loading

library(sampling)
library(Spbsampling)
library(BalancedSampling)
library(WaveSampling)
library(parallel)
library(SDraw)
devtools::load_all(".")


######### set seed and rm environment

rm(list = ls())



findIndex <- function(x1,X){
  final <- c()
  for(i in 1:nrow(x1)){
    id <- seq(1,nrow(X),1)
    for(j in 1:ncol(x1)){
      out = which(x1[i,j] == X[id,j])
      if(length(out) > 1){
        id <- out
      }else if(length(out) == 1 & j == 1){
        final <- c(final,out)
        break;
      }else if(length(out) == 1){
        final <- c(final,id[out])
        break;
      }
    }
  }
  return(final)
}



source("C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/balseq/sim.R")
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


X_ns <- readRDS("C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/balseq/X_ns.rds")


######### Auxiliary variables and variable of interest


p <- 5 
n <- 50

# pik 

pik <- list(eq = rep(n/N,N),
            uneq = inclusionprobabilities(runif(N),n))

# auxiliary variables

Xaux_tmp <- matrix(rgamma(N*p,1,1/5),ncol = p)


Xaux <- list(eq = cbind(pik$eq,Xaux_tmp),
             uneq = cbind(pik$uneq,Xaux_tmp))

# Xaux <- list(eq = cbind(pik$eq),
             # uneq = cbind(pik$uneq))


# variable of interest

beta <- c(1,2,3,4,5)

y_tmp <- list(y1 = Xaux_tmp%*%beta + rnorm(N,0,20),
          y2 = Xaux_tmp%*%beta^2 + rnorm(N,0,2))
          

y <- rep(list(0),2)
names(y) <- c("eq","uneq")
y[[1]] <- list(y1 = y_tmp[[1]],
               y2 = y_tmp[[2]])
y[[2]] <- list(y1 = y_tmp[[1]],
               y2 = y_tmp[[2]])

pairs(cbind(Xaux_tmp,y[[1]][[1]]))
pairs(cbind(Xaux_tmp,y[[1]][[2]]))


### CHECK that hip works
# X_sp <- SpatialPoints(X_pp)
# tmp1 = try(SDraw::hip.point(x = X_sp,
#                             n = sum(pik$eq),
#                             plot.lattice = FALSE),
#            silent = TRUE)
# tmp1 <- coordinates(tmp1)
# s_01 <- rep(0,N)
# s_01[findIndex(tmp1,X_pp)] <- 1
# plot(X_pp)
# lines(X_pp[s_01 == 1,1],X_pp[s_01 == 1,2],pch = 16,type = "p")
# lines(tmp1,pch = 16,type = "p",col = "red")

######### Distance matrix


D_ns <- as.matrix(proxy::dist(X_ns))
D_pp <- as.matrix(proxy::dist(X_pp))

D <- list(pp = D_pp,
          ns = D_ns)



######### Xspread

Xspread <- list(pp = X_pp,
                ns = X_ns)


######### Sampling functions 



f <- list(balseq,
          lcube,
          lpm1,
          pwd,
          hip.point,
          wave,
          srswor)

names(f) <- c("stream",
              "lcube",
              "lpm1",
              "pwd",
              "hip",
              "wave",
              "srswor")

######### Simulation


test <- sim(1,f,Xaux,Xspread$pp,D$pp,pik,y)
test <- sim(1,f,Xaux,Xspread$ns,D$ns,pik,y)
test



SIM <- 1000

cl <- parallel::makeCluster(parallel::detectCores())
clusterEvalQ(cl,{
  library(WaveSampling)
  library(BalancedSampling)
  library(sampling)
  library(SDraw)
  library(Spbsampling)
  devtools::load_all("C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling")
})

start <- Sys.time()

res_pp = parallel::parLapply(cl = cl,
                          X = 1:SIM,
                          fun =  sim,
                          f = f,
                          Xaux = Xaux,
                          Xspread = Xspread$pp,
                          D = D$pp,
                          pik = pik,
                          y = y)

print(Sys.time() -  start)
start <- Sys.time()
res_ns = parallel::parLapply(cl = cl,
                             X = 1:SIM,
                             fun =  sim,
                             f = f,
                             Xaux = Xaux,
                             Xspread = Xspread$ns,
                             D = D$ns,
                             pik = pik,
                             y = y)

print(Sys.time() -  start)


res <- list(res_pp,res_ns)
names(res) <- c("pp","ns")
 

############# save RDS and read RDS

saveRDS(res,file = "C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/balseq/res_28032022.rds")
res <- readRDS(file = "C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/balseq/res_28032022.rds")


############# Creation table
l0 <- c("pp","ns")
l1 <- c("eq","uneq")
for(u in 1:length(l0)){
  for(i in 1:length(l1)){
    assign(paste(l0[u],"spread",l1[i],sep = "_"),lapply(res[[u]],function(x){return(x$spread[[i]])}))
    assign(paste(l0[u],"v",l1[i],sep = "_"),lapply(res[[u]],function(x){return(x$var[[i]])}))
  }
}

ns_spread_eq <- colMeans(do.call(rbind,lapply(ns_spread_eq,unlist)))
ns_spread_uneq <- colMeans(do.call(rbind,lapply(ns_spread_uneq,unlist)))

pp_spread_eq <- colMeans(do.call(rbind,lapply(pp_spread_eq,unlist)))
pp_spread_uneq <- colMeans(do.call(rbind,lapply(pp_spread_uneq,unlist)))

ns_v_eq <- colMeans(do.call(rbind,lapply(ns_v_eq,unlist)))
ns_v_uneq <- colMeans(do.call(rbind,lapply(ns_v_uneq,unlist)))

pp_v_eq <- colMeans(do.call(rbind,lapply(pp_v_eq,unlist)))
pp_v_uneq <- colMeans(do.call(rbind,lapply(pp_v_uneq,unlist)))


o_sb <- grepl(pattern= 'sb',names(ns_spread_eq))
o_IB <- grepl(pattern= 'IB',names(ns_spread_eq))
o_y1 <- grepl(pattern= 'y1',names(ns_v_eq))
o_y2 <- grepl(pattern= 'y1',names(ns_v_eq))

v_res_ns <- t(rbind(c(ns_v_eq[o_y1]/ns_v_eq[o_y1][7]*100,ns_v_uneq[o_y1]/ns_v_uneq[o_y1][7]*100),
                 c(ns_v_eq[o_y2]/ns_v_eq[o_y2][7]*100,ns_v_uneq[o_y2]/ns_v_uneq[o_y2][7]*100)))
v_res_pp <- t(rbind(c(pp_v_eq[o_y1]/pp_v_eq[o_y1][7]*100,pp_v_uneq[o_y1]/pp_v_uneq[o_y1][7]*100),
                 c(pp_v_eq[o_y2]/pp_v_eq[o_y2][7]*100,pp_v_uneq[o_y2]/pp_v_uneq[o_y2][7]*100)))

v_res_pp
v_res_ns

spread_res_ns <- t(rbind(c(ns_spread_eq[o_sb],ns_spread_uneq[o_sb])
,c(ns_spread_eq[o_IB],ns_spread_uneq[o_IB])))
spread_res_pp <- t(rbind(c(pp_spread_eq[o_sb],pp_spread_uneq[o_sb])
,c(pp_spread_eq[o_IB],pp_spread_uneq[o_IB])))
spread_res_pp
spread_res_ns


tab_ns <- cbind(v_res_ns,spread_res_ns)
tab_pp <- cbind(v_res_pp,spread_res_pp)
colnames(tab_ns) <- colnames(tab_pp) <- c("y1","y2","sb","IB")
rownames(tab_ns) <- rownames(tab_pp) <- rep(names(f),2)
tab <- rbind(tab_ns,tab_pp)
names_tab <- rep(c("Proposed Method","Local Cube","Local Pivotal","Proportional within distance","Halton Iterative","Wave","Max entropy"),times = 4)
tab = cbind(names_tab,tab)
tab 














library(kableExtra)
kable(tab, format = "latex",digits = 3, booktabs = T, caption = "Ratio of simulated variance of the Horvitz-Thompson estimator for the two variables of interests. Two measures of spatial balance are displayed for each design.",row.names = FALSE,escape = FALSE) %>%
  add_header_above(c(" " = 1,"$v_{sim}/v_{sim}^{CP}$ " = 2, "Spread measures" = 2),escape = F) %>%
  pack_rows("Neyman-Scott process", 1, 12,latex_gap_space = "2em") %>%
  group_rows("Equal",1,6,escape = F,bold = F,latex_gap_space = "1ex")%>%
  group_rows("Unequal",7,12,escape = F,bold = F,latex_gap_space = "1ex")%>%
  pack_rows("Complete Spatial Randomness", 13, 24,latex_gap_space = "2em") %>%
  group_rows("Equal",13,18,escape = F,bold = F,latex_gap_space = "1ex")%>%
  group_rows("Unequal",19,24,escape = F,bold = F,latex_gap_space = "1ex")


rep(c("Proposed Method","Local Cube","Local Pivotal","Proportional within distance","Halton Iterative","Max entropy"),times = 4)
