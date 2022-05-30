rm(list = ls())

######### Package loading
library(sampling)
library(Spbsampling)
library(BalancedSampling)
library(WaveSampling)
library(parallel)
library(sp)
library(SDraw)
devtools::load_all(".")

######### set seed and rm environment


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
  
data("meuse")
meuse <- na.omit(meuse)
  
N <- nrow(meuse)
n <- 50

X <- as.matrix(meuse[,1:2])

######### Auxiliary variables and variable of interest

# pik 
pik_uneq <- inclusionprobabilities(meuse$copper,n)
index_uneq <- which(pik_uneq < (1 - 1e-7) & pik_uneq > (1e-7))
pik_eq <- rep(n/N,N)
pik <- list(eq = pik_eq, 
            uneq = pik_uneq[index_uneq])
  
# auxiliary variables
# Xaux_tmp <- as.matrix(meuse[,c(5,6,7,9)])
Xaux_tmp <- as.matrix(meuse[,which(colnames(meuse) == "copper" | colnames(meuse) == "elev" | colnames(meuse) == "om")])

Xaux <- list(eq = cbind(pik$eq,Xaux_tmp),
             uneq = cbind(pik_uneq[index_uneq],Xaux_tmp[index_uneq,]))


# variable of interest
  
y_tmp <- list(y1 = meuse$lead,
              y2 = meuse$zinc,
              y3 = meuse$cadmium)


y_eq <- list(y1 = meuse$lead,
             y2 = meuse$zinc,
             y3 = meuse$cadmium)
y_uneq <- list(y1 = meuse$lead[index_uneq],
               y2 = meuse$zinc[index_uneq],
               y3 = meuse$cadmium[index_uneq])

# y <- rep(list(0),2)
# names(y) <- c("eq","uneq")
# y$eq <- y_eq
# y$uneq <- y_uneq

pairs(cbind(Xaux_tmp,y_eq$y3))
  

######### Distance matrix

X_eq <- as.matrix(scale(as.matrix(meuse[,1:2])))
X_uneq <- as.matrix(scale(as.matrix(meuse[index_uneq,1:2])))

D_uneq <- as.matrix(proxy::dist(X_uneq))
D_eq <- as.matrix(proxy::dist(X_eq))

D <- list(eq = D_eq,
          uneq = D_uneq)
  
######### Xspread

Xspread <- list(eq = X_eq,
                uneq = X_uneq)



######### Sampling functions 



# system.time(wave(Xspread$amphib,pik = pik$eq,comment = TRUE))

f_eq <- list(balseq,
             lcube,
             cube,
             lpm1,
             pwd,
             hip.point,
             wave,
             srswor)

names(f_eq) <- c("stream",
                 "lcube",
                 "cube",
                 "lpm1",
                 "pwd",
                 "hip",
                 "wave",
                 "srswor")

f_eq_name <- c("Proposed Method",
               "Local Cube",
               "Cube Method",
               "Local Pivotal",
               "Proportional within distance",
               "Halton Iterative Partitioning",
               "Wave",
               "Simple random sampling")


f_uneq <- list(balseq,
               lcube,
               cube,
               lpm1,
               wave,
               cps)

names(f_uneq) <- c("stream",
                   "lcube",
                   "cube",
                   "lpm1",
                   "wave",
                   "cps")

f_uneq_name <- c("Proposed Method",
                 "Local Cube",
                 "Cube Method",
                 "Local Pivotal",
                 "Wave",
                 "Conditional Poisson sampling")


# sim_design(1,f = f_uneq,
#            Xaux = Xaux$uneq,
#            Xspread = Xspread$uneq,
#            D = D$uneq,
#            pik = pik$uneq,
#            y = y_uneq)
# 
# sim_design(1,f = f_eq,
#            Xaux = Xaux$eq,
#            Xspread = Xspread$eq,
#            D = D$eq,
#            pik = pik$eq,
#            y = y$eq)

######### Simulation
  
  
  
  
SIM <- 1000
  
cl <- parallel::makeCluster(parallel::detectCores())
clusterEvalQ(cl,{
  library(WaveSampling)
  library(BalancedSampling)
  library(sampling)
  library(Spbsampling)
  devtools::load_all("C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling")
})



start <- Sys.time()
res_meuse_eq = parallel::parLapply(cl = cl,
                                    X = 1:SIM,
                                    fun =  sim_design,
                                    f = f_eq,
                                    Xaux = Xaux$eq,
                                    Xspread = Xspread$eq,
                                    D = D$eq,
                                    pik = pik$eq,
                                    y = y_eq)
print(Sys.time() -  start)

start <- Sys.time()
res_meuse_uneq = parallel::parLapply(cl = cl,
                                      X = 1:SIM,
                                      fun =  sim_design,
                                      f = f_uneq,
                                      Xaux = Xaux$uneq,
                                      Xspread = Xspread$uneq,
                                      D = D$uneq,
                                      pik = pik$uneq,
                                      y = y_uneq)
print(Sys.time() -  start)

stopCluster(cl)

# set.seed(41) # fail for hip.point
# sim(1,f,Xaux,Xspread$meuse,D$meuse,pik,y)


# res_meuse <- list()
# for(i in 1:SIM){
#   print(i)
#   set.seed(i)
#   res_meuse[[i]] <- sim(1,f,Xaux,Xspread$meuse,D$meuse,pik,y)
# }

res_meuse <- list(eq = res_meuse_eq,
                   uneq = res_meuse_uneq)

res <- list(res_meuse)
names(res) <- c("meuse")
 

############# save RDS and read RDS

saveRDS(res,file = "C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/balseq/res_meuse_09052022.rds")
res <- readRDS(file = "C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/balseq/res_meuse_09052022.rds")


############# Creation table
l0 <- c("meuse")
l1 <- c("eq","uneq")
for(u in 1:length(l0)){
  for(i in 1:length(l1)){

    assign(paste(l0[u],l1[i],"vApp",sep = "_"),lapply(res[[u]][[i]],function(x){return(x$vApp)}))
    assign(paste(l0[u],l1[i],"spread",sep = "_"),lapply(res[[u]][[i]],function(x){return(x$spread)}))
    assign(paste(l0[u],l1[i],"v",sep = "_"),lapply(res[[u]][[i]],function(x){return(x$v)}))
    assign(paste(l0[u],l1[i],"vEst",sep = "_"),lapply(res[[u]][[i]],function(x){return(x$vEst)}))
    assign(paste(l0[u],l1[i],"HT",sep = "_"),lapply(res[[u]][[i]],function(x){return(x$HT)}))

  }
}


meuse_eq_spread <- colMeans(do.call(rbind,lapply(meuse_eq_spread,unlist)))
meuse_uneq_spread <- colMeans(do.call(rbind,lapply(meuse_uneq_spread,unlist)))

meuse_eq_vEst <- colMeans(do.call(rbind,lapply(meuse_eq_vEst,unlist)))
meuse_uneq_vEst <- colMeans(do.call(rbind,lapply(meuse_uneq_vEst,unlist)))

meuse_eq_v <- colMeans(do.call(rbind,lapply(meuse_eq_v,unlist)))
meuse_uneq_v <- colMeans(do.call(rbind,lapply(meuse_uneq_v,unlist)))

o_eq_sb <- grepl(pattern= 'sb',names(meuse_eq_spread))
o_eq_IB <- grepl(pattern= 'IB',names(meuse_eq_spread))
o_uneq_sb <- grepl(pattern= 'sb',names(meuse_uneq_spread))
o_uneq_IB <- grepl(pattern= 'IB',names(meuse_uneq_spread))


o_eq_y1 <- grepl(pattern= 'y1',names(meuse_eq_v))
o_uneq_y1 <- grepl(pattern= 'y1',names(meuse_uneq_v))

o_eq_y2 <- grepl(pattern= 'y2',names(meuse_eq_v))
o_uneq_y2 <- grepl(pattern= 'y2',names(meuse_uneq_v))

o_eq_y3 <- grepl(pattern= 'y3',names(meuse_eq_v))
o_uneq_y3 <- grepl(pattern= 'y3',names(meuse_uneq_v))




meuse_v <- cbind(c(meuse_eq_v[o_eq_y1]/meuse_eq_v[o_eq_y1][which(grepl('srswor',names(meuse_eq_v[o_eq_y1])) == TRUE)]*100,
                   meuse_uneq_v[o_uneq_y1]/meuse_uneq_v[o_uneq_y1][which(grepl('cps',names(meuse_uneq_v[o_uneq_y1])) == TRUE)]*100),
                 c(meuse_eq_v[o_eq_y2]/meuse_eq_v[o_eq_y2][which(grepl('srswor',names(meuse_eq_v[o_eq_y2])) == TRUE)]*100,
                   meuse_uneq_v[o_uneq_y2]/meuse_uneq_v[o_uneq_y2][which(grepl('cps',names(meuse_uneq_v[o_uneq_y2])) == TRUE)]*100),
                 c(meuse_eq_v[o_eq_y3]/meuse_eq_v[o_eq_y3][which(grepl('srswor',names(meuse_eq_v[o_eq_y3])) == TRUE)]*100,
                   meuse_uneq_v[o_uneq_y3]/meuse_uneq_v[o_uneq_y3][which(grepl('cps',names(meuse_uneq_v[o_uneq_y3])) == TRUE)]*100))

meuse_vEst <- cbind(c(meuse_eq_vEst[o_eq_y1]/meuse_eq_v[o_eq_y1]*100, meuse_uneq_vEst[o_uneq_y1]/meuse_uneq_v[o_uneq_y1]*100),
                 c(meuse_eq_vEst[o_eq_y2]/meuse_eq_v[o_eq_y2]*100, meuse_uneq_vEst[o_uneq_y2]/meuse_uneq_v[o_uneq_y2]*100),
                 c(meuse_eq_vEst[o_eq_y3]/meuse_eq_v[o_eq_y3]*100, meuse_uneq_vEst[o_uneq_y3]/meuse_uneq_v[o_uneq_y3]*100))
meuse_vEst



meuse_spread <- t(rbind(c(meuse_eq_spread[o_eq_sb],meuse_uneq_spread[o_uneq_sb])
                         ,c(meuse_eq_spread[o_eq_IB],meuse_uneq_spread[o_uneq_IB])))
meuse_spread



tab_meuse <- round(cbind(meuse_v,meuse_vEst,meuse_spread),3)
# colnames(tab_meuse) <- c("y3","v1","sb","IB")
rownames(tab_meuse) <- c(names(f_eq),names(f_uneq))
# rownames(tab_meuse) <-  rep(names(f),2)
# names_tab <- rep(c("Proposed Method","Local Cube","Local Pivotal","Proportional within distance","Max entropy"),times = 2)
# tab_meuse = cbind(rownames(tab_meuse),tab_meuse)
tab_meuse

kable(tab_meuse, format = "latex",digits = 3, booktabs = T, caption = "blab bla",row.names = FALSE,escape = FALSE) %>%
  add_header_above(c(" " = 1,"$v_{sim}/v_{sim}^{CP}$ " = 2, "Spread measures" = 2),escape = F) %>%
  pack_rows("Neyman-Scott process", 1, 5,latex_gap_space = "2em") %>%
  group_rows("Equal",1,5,escape = F,bold = F,latex_gap_space = "1ex")%>%
  group_rows("Unequal",6,10,escape = F,bold = F,latex_gap_space = "1ex")