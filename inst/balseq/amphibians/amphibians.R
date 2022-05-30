# source load data



source("./inst/balseq/amphibians/data_loading.R")

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

######### Auxiliary variables and variable of interest

N <- nrow(amphib)
n <- 100

# pik 

pik_eq <- rep(n/N,N)
pik_uneq <- inclusionprobabilities(amphib$AREA,n)
index_uneq <- which(pik_uneq < (1 - 1e-7) & pik_uneq > (1e-7))

pik <- list(eq = pik_eq,
            uneq = pik_uneq[index_uneq])

# auxiliary variables

bioreg <- disj(amphib$NO_BIOREG)
colnames(bioreg) <- seq(1,6,1)

# Xaux_tmp <- data.frame(altitude = amphib$ALTITUDE,
#                        area = amphib$AREA,
#                        y1 = amphib$y1,
#                        y2 = amphib$y2,
#                        y3 = amphib$y3,
#                        y4 = amphib$y4)

# Xaux_tmp <- data.frame(altitude = amphib$ALTITUDE,
                       # area = amphib$AREA)
Xaux_tmp <- data.frame(area = amphib$ALTITUDE)


Xaux_tmp <- cbind(Xaux_tmp,
                  bioreg)

# Xaux_eq <- list(eq = as.matrix(cbind(pik_eq,Xaux_tmp)))
# Xaux_uneq <- list(uneq = as.matrix(cbind(pik_uneq,Xaux_tmp)))

Xaux <- list(eq = as.matrix(cbind(pik$eq,Xaux_tmp)),
             uneq = as.matrix(cbind(pik$uneq,Xaux_tmp[index_uneq,])))

# variable of interest

# y_tmp <- list(y1 = amphib$y1 + amphib$y2 + amphib$y3 +amphib$y4)
y_tmp <- list(y1 = amphib$diversity)

y_eq <- list(y1 = y_tmp[[1]])
y_uneq <- list(y1 = y_tmp[[1]][index_uneq])

y <- rep(list(0),2)
names(y) <- c("eq","uneq")
y$eq <- y_eq
y$uneq <- y_uneq
# y[[1]] <- list(y1 = y_tmp[[1]])
# y[[2]] <- list(y1 = y_tmp[[1]])


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

X_eq <- as.matrix(scale(data.frame(x = amphib$COORD_X,y = amphib$COORD_Y)))
X_uneq <- as.matrix(scale(data.frame(x = amphib$COORD_X[index_uneq],y = amphib$COORD_Y[index_uneq])))

D_uneq <- as.matrix(proxy::dist(X_uneq))
D_eq <- as.matrix(proxy::dist(X_eq))

D <- list(eq = D_eq,
          uneq = D_uneq)

######### Xspread

# Xspread <- list(amphib = X)

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
          # wave,
          srswor)

names(f_eq) <- c("stream",
              "lcube",
              "cube",
              "lpm1",
              "pwd",
              "hip",
              # "wave",
              "srswor")

f_eq_name <- c("Proposed Method",
               "Local Cube",
               "Cube Method",
               "Local Pivotal",
               "Proportional within distance",
               "Halton Iterative Partitioning",
               # "Wave",
               "Simple random sampling")


f_uneq <- list(balseq,
             lcube,
             cube,
             lpm1,
             # wave,
             cps)

names(f_uneq) <- c("stream",
                 "lcube",
                 "cube",
                 "lpm1",
                 # "wave",
                 "cps")

f_uneq_name <- c("Proposed Method",
                 "Local Cube",
                 "Cube Method",
                 "Local Pivotal",
                 # "Wave",
                 "Conditional Poisson sampling")


sim_design(1,f = f_uneq,
           Xaux = Xaux$uneq,
           Xspread = Xspread$uneq,
           D = D$uneq,
           pik = pik$uneq,
           y = y$uneq)

sim_design(1,f = f_eq,
           Xaux = Xaux$eq,
           Xspread = Xspread$eq,
           D = D$eq,
           pik = pik$eq,
           y = y$eq)


######### Simulation
# 
# s <- BalancedSampling::cube(pik$eq,Xaux$eq)
# s_01 <- rep(0,nrow(Xaux$eq))
# s_01[s] <- 1
# vEst(Xaux$uneq,pik$uneq,y$uneq$y1,s_01)

# sim_design(1,f_eq,Xaux$eq,Xspread$amphib,D$amphib,pik$eq,y = y$eq)
# sim_design(1,f_uneq,Xaux$uneq,Xspread$amphib,D$amphib,pik$uneq,y = y$uneq)

# test <- sim(1,f,Xaux,Xspread$amphib,D$amphib,pik,y)
# test

SIM <- 10000

cl <- parallel::makeCluster(parallel::detectCores())
clusterEvalQ(cl,{
  library(WaveSampling)
  library(BalancedSampling)
  library(sampling)
  # library(SDraw)
  library(Spbsampling)
  devtools::load_all("C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling")
})


# res_amphib_eq <- list()
# for(i in 1:100){
#   # set.seed(i)
#   print(i)
#   res_amphib_eq[[i]] <- sim_design(1,f_eq,Xaux$eq,Xspread$amphib,D$amphib,pik$eq,y = y$eq)
# }


start <- Sys.time()
res_amphib_eq = parallel::parLapply(cl = cl,
                                    X = 1:SIM,
                                    fun =  sim_design,
                                    f = f_eq,
                                    Xaux = Xaux$eq,
                                    Xspread = Xspread$eq,
                                    D = D$eq,
                                    pik = pik$eq,
                                    y = y$eq)
print(Sys.time() -  start)

start <- Sys.time()
res_amphib_uneq = parallel::parLapply(cl = cl,
                                    X = 1:SIM,
                                    fun =  sim_design,
                                    f = f_uneq,
                                    Xaux = Xaux$uneq,
                                    Xspread = Xspread$uneq,
                                    D = D$uneq,
                                    pik = pik$uneq,
                                    y = y$uneq)
print(Sys.time() -  start)

stopCluster(cl)

res_amphib <- list(eq = res_amphib_eq,
                   uneq = res_amphib_uneq)

res <- list(res_amphib)
names(res) <- c("amphib")


############# save RDS and read RDS

saveRDS(res,file = "C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/balseq/res_amphib_25052022_2.rds")
# res <- readRDS(file = "C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/balseq/res_amphib_11052022_2.rds")


############# Creation table
# l0 <- c("amphib")
# l1 <- c("eq","uneq")
# for(u in 1:length(l0)){
#   for(i in 1:length(l1)){
# 
#     assign(paste(l0[u],l1[i],"vApp",sep = "_"),lapply(res[[u]][[i]],function(x){return(x$vApp)}))
#     assign(paste(l0[u],l1[i],"spread",sep = "_"),lapply(res[[u]][[i]],function(x){return(x$spread)}))
#     assign(paste(l0[u],l1[i],"v",sep = "_"),lapply(res[[u]][[i]],function(x){return(x$v)}))
#     assign(paste(l0[u],l1[i],"vEst",sep = "_"),lapply(res[[u]][[i]],function(x){return(x$vEst)}))
#     assign(paste(l0[u],l1[i],"HT",sep = "_"),lapply(res[[u]][[i]],function(x){return(x$HT)}))
#     assign(paste(l0[u],l1[i],"dev",sep = "_"),lapply(res[[u]][[i]],function(x){return(x$dev)}))
# 
#   }
# }
# 
# ##### RELATIVE DEVIATION
# 
# amphib_dev <- rep(list(0),2)
# names(amphib_dev) <- c("eq","uneq")
# 
# amphib_dev$eq <- amphib_eq_dev[[1]]
# amphib_dev$uneq <- amphib_uneq_dev[[1]]
# 
# 
# for(i in 2:SIM){
#   amphib_dev$eq  <- amphib_dev$eq + amphib_eq_dev[[i]]/SIM
#   amphib_dev$uneq  <- amphib_dev$uneq + amphib_uneq_dev[[i]]/SIM
# }
# amphib_dev$eq
# amphib_dev$uneq
# 
# dev <- round(as.matrix(rbind(amphib_dev$eq,
#                              amphib_dev$uneq)),3)
# rownames(dev) <- c(f_eq_name,f_uneq_name)
# 
# 
# amphib_eq_spread <- colMeans(do.call(rbind,lapply(amphib_eq_spread,unlist)))
# amphib_uneq_spread <- colMeans(do.call(rbind,lapply(amphib_uneq_spread,unlist)))
# 
# amphib_eq_vEst <- colMeans(do.call(rbind,lapply(amphib_eq_vEst,unlist)))
# amphib_uneq_vEst <- colMeans(do.call(rbind,lapply(amphib_uneq_vEst,unlist)))
# 
# amphib_eq_v <- colMeans(do.call(rbind,lapply(amphib_eq_v,unlist)))
# amphib_uneq_v <- colMeans(do.call(rbind,lapply(amphib_uneq_v,unlist)))
# 
# o_eq_sb <- grepl(pattern= 'sb',names(amphib_eq_spread))
# o_eq_IB <- grepl(pattern= 'IB',names(amphib_eq_spread))
# o_uneq_sb <- grepl(pattern= 'sb',names(amphib_uneq_spread))
# o_uneq_IB <- grepl(pattern= 'IB',names(amphib_uneq_spread))
# 
# o_eq_y1 <- grepl(pattern= 'y1',names(amphib_eq_v))
# o_uneq_y1 <- grepl(pattern= 'y1',names(amphib_uneq_v))
# 
# amphib_v <- cbind(c(amphib_eq_v[o_eq_y1]/amphib_eq_v[o_eq_y1][which(grepl('srswor',names(amphib_eq_v[o_eq_y1])) == TRUE)]*100,
#                           amphib_uneq_v[o_uneq_y1]/amphib_uneq_v[o_uneq_y1][which(grepl('cps',names(amphib_uneq_v[o_uneq_y1])) == TRUE)]*100))
# 
# amphib_vEst <- cbind(c(amphib_eq_vEst/amphib_eq_v,amphib_uneq_vEst/amphib_uneq_v))
# 
# amphib_spread <- t(rbind(c(amphib_eq_spread[o_eq_sb],amphib_uneq_spread[o_uneq_sb])
#                          ,c(amphib_eq_spread[o_eq_IB],amphib_uneq_spread[o_uneq_IB])))
# amphib_spread
# 
# 
# 
# tab_amphib <- round(cbind(amphib_v,amphib_vEst,amphib_spread),3)
# colnames(tab_amphib) <- c("y1","v1","sb","IB")
# rownames(tab_amphib) <- c(names(f_eq),names(f_uneq))
# # rownames(tab_amphib) <-  rep(names(f),2)
# # names_tab <- rep(c("Proposed Method","Local Cube","Local Pivotal","Proportional within distance","Max entropy"),times = 2)
# # tab_amphib = cbind(rownames(tab_amphib),tab_amphib)
# tab_amphib
# 
# kable(tab_amphib, format = "latex",digits = 3, booktabs = T, caption = "blab bla",row.names = FALSE,escape = FALSE) %>%
#   add_header_above(c(" " = 1,"$v_{sim}/v_{sim}^{CP}$ " = 2, "Spread measures" = 2),escape = F) %>%
#   pack_rows("Neyman-Scott process", 1, 5,latex_gap_space = "2em") %>%
#   group_rows("Equal",1,5,escape = F,bold = F,latex_gap_space = "1ex")%>%
#   group_rows("Unequal",6,10,escape = F,bold = F,latex_gap_space = "1ex")
# 
# 
# 
# 
