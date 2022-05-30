

######### Package loading

library(sampling)
library(Spbsampling)
library(BalancedSampling)
library(WaveSampling)
library(parallel)
library(ggplot2)
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

gauss2d <- function(x,y,x0,y0,a,b,c,A){
  return(A*exp(-(a*(x-x0)^2 + 2*b*(x-x0)*(y-y0) + c*(y-y0)^2)))
}

# x <- seq(0,1,0.005)
# y <- seq(0,1,0.005)
# dat <- expand.grid(x,y)
# 
# gauss2d(x,y)
# library(ggplot2)
# ggplot() +
#   geom_point(data = data.frame(x = dat$Var1,y = dat$Var2,z = gauss2d(dat$Var1,dat$Var2,0.5,0.5,1,-0.3,1,1,10) ),aes(x = x, y = y,color = z))
# 

source("C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/balseq/sim.R")
# set.seed(1234)
set.seed(3)

######### Data generation loading

N <- 200

# CSR


pp <- spatstat.random::rpoispp(N)
X_pp <- matrix(cbind(pp$x,pp$y),nrow = pp$n,ncol = 2)
N_pp <- nrow(X_pp)
while(N_pp != N){
  pp <- spatstat.random::rpoispp(N)
  X_pp <- matrix(cbind(pp$x,pp$y),nrow = pp$n,ncol = 2)
  N_pp <- nrow(X_pp)
}

######### Auxiliary variables and variable of interest


grad_pp <- gauss2d(X_pp[,1],X_pp[,2],0.7,0.5,3,-0.5,3,15) 

X1 <- rnorm(N)
X2 <- rexp(N,rate = 1)
X3 <- rgamma(N,3,1)
X4 <- rbeta(N,2,5)
X5 <- runif(N,0,3)

n <- 20


pik_tmp_uneq_pp <- inclusionprobabilities(grad_pp,n)

pik_pp <- list(eq = rep(n/N,N),
               uneq = pik_tmp_uneq_pp)


# auxiliary variables

Xaux_tmp <- cbind(X1,X2,X3,X4,X5)


Xaux_pp <- list(eq = cbind(pik_pp$eq,Xaux_tmp),
                uneq = cbind(pik_pp$uneq,Xaux_tmp))

# variable of interest


beta <- c(1,1,1,1,1)
y_pp <- list(y1 = grad_pp + Xaux_tmp%*%beta + rnorm(N,0,0.1))


x <- seq(0,1,length.out = 500)
sha <- expand.grid(x,x)

gauss2d(X_pp[,1],X_pp[,2],0.7,0.5,3,-0.5,3,15) 

ggplot() +
  geom_raster(data = data.frame(x = sha[,1],y = sha[,2],z = gauss2d(sha[,1],sha[,2],0.7,0.5,3,-0.5,3,15),shape = 23),aes(x = x,y =y ,fill = z))+
  scale_fill_gradient("Spatial Correlation",high = "#1c658c",low = "#d8d2cb")+
  geom_point(data = data.frame(x = X_pp[,1],y = X_pp[,2],z = pik_pp$uneq ),aes(x = x, y = y,size = z),shape = 16)+
  scale_size(range = c(0.5, 3.5)) +  
  theme_minimal() +
  theme(
    text = element_text(family="sans",color = "black",size = 9),
    panel.spacing = unit(2, "lines"),
    # title
    plot.title = element_text(hjust = 0.5,size = 9),
    # axes
    axis.line=element_blank(),
    axis.text = element_blank(),
    axis.ticks=element_blank(),
    # legend
    legend.position="bottom",
    legend.title = element_text(size = 9,vjust = +1.0),
    legend.key.size = unit(0.3, "cm"),
    legend.key.width = unit(1.1,"cm") ,
    # background colors
    panel.background=element_blank(),
    panel.border=element_rect(colour = "black",fill = "transparent"),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    # keep edge black facet_wrap
    # strip.background = element_rect(fill="white"),
    strip.text =element_text(color = "black",size = 8)
  )


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

D_pp <- as.matrix(proxy::dist(X_pp))
D <- list(pp = D_pp)



######### Xspread

Xspread <- list(pp = X_pp)


######### Sampling functions 


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

######### Simulation




# sim_design(1,f_eq,Xaux_pp$eq,Xspread$pp,D$pp,pik_pp$eq,y_pp)
# sim_design(1,f_uneq,Xaux_pp$uneq,Xspread$pp,D$pp,pik_pp$uneq,y_pp)



# test <- sim(1,f,Xaux,Xspread$pp,D$pp,pik,y)
# test <- sim(1,f,Xaux,Xspread$ns,D$ns,pik,y)
# test



########################### vEst pas normal

# SIM <- 100
# v_estimated_srs <- v_estimated <- v_estimated2 <- v_simulated <- 0
# 
# HT <- c()
# for(i in 1:SIM){
#   print(i)
#   s_srs <- sampling::srswor(sum(pik$eq),nrow(Xaux$eq))
#   # s <- BalancedSampling::cube(pik$eq,Xaux$eq,Xspread$pp)
#   s <- BalancedSampling::cube(pik$eq,Xaux$eq)
#   # print(s)
#   s_01 <- rep(0,nrow(Xaux$eq))
#   s_01[s] <- 1
#   # HT_srs[i] <- sum(y$eq[s]/pik$eq[s])
#   # HT_balseq[i] <- sum(y$eq[s]/pik$eq[s])
#   
#   v_estimated <- v_estimated + vEst(Xaux$eq,pik$eq,y$eq$y1,s_01)
#   # v_estimated2 <- v_estimated2 + varEst2(Xaux$eq,pik$eq,y$eq$y1,s_01)
#   v_estimated2 <- v_estimated2 + vEst2(Xaux$eq[s,],pik$eq[s],y$eq$y1[s])
#   v_estimated_srs <- v_estimated_srs + vEst(Xaux$eq,pik$eq,y$eq$y1,s_srs)
#   v_simulated <- v_simulated + (sum(y$eq$y1[s]/pik$eq[s]) - sum(y$eq$y1))^2
#   v_simulated_srs <- v_simulated_srs + (sum(y$eq$y1[s_srs]/pik$eq[s_srs]) - sum(y$eq$y1))^2
# }
# 
# 
# v_estimated/SIM
# v_estimated2/SIM
# v_estimated_srs/SIM
# v_simulated/SIM
# v_simulated_srs/SIM





###########################

SIM <- 10000


########################### SIMU NON PARALLE


# res_pp_eq <- res_pp_uneq <- res_ns_eq <- res_ns_uneq <- list()
# 
# for(i in 1:SIM){
#   print(i)
#   res_pp_eq[[i]] <- sim_design(1,f_eq,Xaux_pp$eq,Xspread$pp,D$pp,pik_pp$eq,y_pp) 
#   res_pp_uneq[[i]] <- sim_design(1,f_uneq,Xaux_pp$uneq,Xspread$pp,D$pp,pik_pp$eq,y_pp) 
#   res_ns_eq[[i]] <- sim_design(1,f_eq,Xaux_ns$eq,Xspread$ns,D$ns,pik_ns$eq,y_ns) 
#   res_ns_uneq[[i]] <- sim_design(1,f_uneq,Xaux_ns$uneq,Xspread$ns,D$ns,pik_ns$eq,y_ns) 
# }


###########################


cl <- parallel::makeCluster(parallel::detectCores())
# cl <- parallel::makeCluster(4)
clusterEvalQ(cl,{
  library(WaveSampling)
  library(BalancedSampling)
  library(sampling)
  # library(SDraw)
  library(Spbsampling)
  devtools::load_all("C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling")
})

start <- Sys.time()
res_pp_eq = parallel::parLapply(cl = cl,
                                X = 1:SIM,
                                fun =  sim_design,
                                f = f_eq,
                                Xaux = Xaux_pp$eq,
                                Xspread = Xspread$pp,
                                D = D$pp,
                                pik = pik_pp$eq,
                                y = y_pp)
print(Sys.time() -  start)
start <- Sys.time()
res_pp_uneq = parallel::parLapply(cl = cl,
                                  X = 1:SIM,
                                  fun =  sim_design,
                                  f = f_uneq,
                                  Xaux = Xaux_pp$uneq,
                                  Xspread = Xspread$pp,
                                  D = D$pp,
                                  pik = pik_pp$uneq,
                                  y = y_pp)
print(Sys.time() -  start)

stopCluster(cl)

res_pp <- list(eq = res_pp_eq,
               uneq = res_pp_uneq)
res <- list(res_pp)
names(res) <- c("pp")


############# save RDS and read RDS

saveRDS(res,file = "C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/balseq/res_perugia_25052022.rds")
# res <- readRDS(file = "C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/inst/balseq/res_artificial_11052022.rds")


########### Creation table
# 
# l0 <- c("pp")
# l1 <- c("eq","uneq")
# for(u in 1:length(l0)){
#   for(i in 1:length(l1)){
#     assign(paste(l0[u],l1[i],"spread",sep = "_"),lapply(res[[u]][[i]],function(x){return(x$spread)}))
#     assign(paste(l0[u],l1[i],"v",sep = "_"),lapply(res[[u]][[i]],function(x){return(x$v)}))
#     assign(paste(l0[u],l1[i],"vEst",sep = "_"),lapply(res[[u]][[i]],function(x){return(x$vEst)}))
#     assign(paste(l0[u],l1[i],"HT",sep = "_"),lapply(res[[u]][[i]],function(x){return(x$HT)}))
#     assign(paste(l0[u],l1[i],"dev",sep = "_"),lapply(res[[u]][[i]],function(x){return(x$dev)}))
#   }
# }
# 
# ##### RELATIVE DEVIATION
# 
# pp_dev <- rep(list(0),2)
# names(pp_dev) <- c("eq","uneq")
# 
# pp_dev$eq <- abs(pp_eq_dev[[1]])/SIM
# pp_dev$uneq <- abs(pp_uneq_dev[[1]])/SIM
# 
# 
# for(i in 2:SIM){
#   pp_dev$eq  <- pp_dev$eq + abs(pp_eq_dev[[i]])/SIM
#   pp_dev$uneq  <- pp_dev$uneq + abs(pp_uneq_dev[[i]])/SIM
# }
# 
#  do.call(rbind,apply(pp_dev$eq,MARGIN = 1,FUN = function(x){x/pp_dev$eq[nrow(pp_dev$eq),]}))*100
# pp_dev$eq/
# pp_dev$uneq
# 
# dev <- round(as.matrix(rbind(pp_dev$eq,
#              pp_dev$uneq)),3)
# 
# dev
# 
# rownames(dev) <- c(f_eq_name,f_uneq_name)
# colnames(dev) <- c("Inclusion Probabilities","$X_1$","$X_2$","$X_3$","$X_4$","$X_5$")
# 
# 
# ############# SPREAD AND VARIANCE
# 
# pp_eq_spread <- colMeans(do.call(rbind,lapply(pp_eq_spread,unlist)))
# pp_uneq_spread <- colMeans(do.call(rbind,lapply(pp_uneq_spread,unlist)))
# 
# 
# pp_eq_v <- colMeans(do.call(rbind,lapply(pp_eq_v,unlist)))
# pp_uneq_v <- colMeans(do.call(rbind,lapply(pp_uneq_v,unlist)))
# 
# pp_eq_vEst <- colMeans(do.call(rbind,lapply(pp_eq_vEst,unlist)))
# pp_uneq_vEst <- colMeans(do.call(rbind,lapply(pp_uneq_vEst,unlist)))
# 
# 
# o_eq_sb <- grepl(pattern= 'sb',names(pp_eq_spread))
# o_uneq_sb <- grepl(pattern= 'sb',names(pp_uneq_spread))
# o_eq_IB <- grepl(pattern= 'IB',names(pp_eq_spread))
# o_uneq_IB <- grepl(pattern= 'IB',names(pp_uneq_spread))
# 
# o_eq_y1 <- grepl(pattern= 'y1',names(pp_eq_v))
# o_uneq_y1 <- grepl(pattern= 'y1',names(pp_uneq_v))
# # o_eq_y2 <- grepl(pattern= 'y2',names(ns_eq_v))
# # o_uneq_y2 <- grepl(pattern= 'y2',names(ns_uneq_v))
# 
# # which(grepl('srswor',names(ns_eq_v[o_eq_y1])))
# # which(grepl('srswor',names(ns_uneq_v[o_uneq_y1])))
# 
# 
# pp_eq_vEst <- pp_eq_vEst/pp_eq_v*100
# pp_uneq_vEst <- pp_uneq_vEst/pp_uneq_v*100
# 
# pp_v_y1 <- c(pp_eq_v[o_eq_y1]/pp_eq_v[o_eq_y1][which(grepl('srswor',names(pp_eq_v[o_eq_y1])))]*100,
#                  pp_uneq_v[o_uneq_y1]/pp_uneq_v[o_eq_y1][which(grepl('cps',names(pp_uneq_v[o_eq_y1])))]*100)
# 
# 
# pp_v <- cbind(pp_v_y1)
# 
# 
# pp_vEst_y1 <- c(pp_eq_vEst[o_eq_y1],pp_uneq_vEst[o_uneq_y1])
# 
# 
# 
# # ns_vEst_y1 <- c(ns_eq_vEst[o_eq_y1]/ns_eq_vEst[o_eq_y1][which(grepl('srswor',names(ns_eq_vEst[o_eq_y1])))]*100,
# #              ns_uneq_vEst[o_uneq_y1]/ns_uneq_vEst[o_eq_y1][which(grepl('cps',names(ns_uneq_vEst[o_eq_y1])))]*100)
# # ns_vEst_y2 <- c(ns_eq_vEst[o_eq_y2]/ns_eq_vEst[o_eq_y2][which(grepl('srswor',names(ns_eq_vEst[o_eq_y2])))]*100,
# #              ns_uneq_vEst[o_uneq_y2]/ns_uneq_vEst[o_eq_y2][which(grepl('cps',names(ns_uneq_vEst[o_eq_y2])))]*100)
# # pp_vEst_y1 <- c(pp_eq_vEst[o_eq_y1]/pp_eq_vEst[o_eq_y1][which(grepl('srswor',names(pp_eq_vEst[o_eq_y1])))]*100,
# #              pp_uneq_vEst[o_uneq_y1]/pp_uneq_vEst[o_eq_y1][which(grepl('cps',names(pp_uneq_vEst[o_eq_y1])))]*100)
# # pp_vEst_y2 <- c(pp_eq_vEst[o_eq_y2]/pp_eq_vEst[o_eq_y2][which(grepl('srswor',names(pp_eq_vEst[o_eq_y2])))]*100,
# #              pp_uneq_vEst[o_uneq_y2]/pp_uneq_vEst[o_eq_y2][which(grepl('cps',names(pp_uneq_vEst[o_eq_y2])))]*100)
# 
# # ns_vEst <- cbind(ns_vEst_y1,ns_vEst_y2)
# 
# pp_vEst <- cbind(pp_vEst_y1)
# # pp_vEst <- cbind(pp_vEst_y1,pp_vEst_y2)
# 
# 
# pp_spread_sb <- c(pp_eq_spread[o_eq_sb],pp_uneq_spread[o_uneq_sb])
# pp_spread_IB <- c(pp_eq_spread[o_eq_IB],pp_uneq_spread[o_uneq_IB])
# 
# 
# 
# pp_spread <- cbind(pp_spread_sb,pp_spread_IB)
# 
# 
# 
# pp <- cbind(pp_v,pp_vEst,pp_spread)
# 
# tab <- rbind(pp)
# tab
# 
# # rownames(tab) <- c(f_eq_name,f_uneq_name)
# # colnames(tab) <- c("$v_1/v_{srs,1}$",
# #                    # "$v_2/v_{srs,2}$",
# #                    "$\\widehat{\\text{var}}(\\widehat{Y}_1)/v_1$",
# #                    # "$\\widehat{\\text{var}}(\\widehat{Y}_2)/v_2$",
# #                    "sb","IB")
# # # colnames(tab) <- c("$v_1/v_{srs,1}$","$v_2/v_{srs,2}$","$\\widehat{\\text{var}}(\\widehat{Y}_1)/v_1$","$\\widehat{\\text{var}}(\\widehat{Y}_2)/v_2$","sb","IB")
# # tab
# # tab <- round(tab,3)
# # 
# # 
# # nfeq <- length(f_eq_name)
# # nfuneq <- length(f_uneq_name)
# # 
# # library(kableExtra)
# # kable(tab, format = "latex",digits = 3, booktabs = T, caption = "Ratio of simulated variance of the Horvitz-Thompson estimator for the two variables of interests. Two measures of spatial balance are displayed for each design.",row.names = TRUE,escape = FALSE) %>%
# #   add_header_above(c(" " = 1,"Simulated Variances " = 2,"Variance Estimators" = 2, "Spread measures" = 2),escape = F) %>%
# #   pack_rows("Neyman-Scott process", 1, nfeq + nfuneq,latex_gap_space = "2em") %>%
# #   group_rows("Equal",1,nfeq,escape = F,bold = F,latex_gap_space = "1ex")%>%
# #   group_rows("Unequal",nfeq + 1, nfeq + nfuneq,escape = F,bold = F,latex_gap_space = "1ex")%>%
# #   pack_rows("Complete Spatial Randomness", nfeq + nfuneq + 1, 2*(nfeq + nfuneq),latex_gap_space = "2em") %>%
# #   group_rows("Equal",nfeq + nfuneq + 1,2*nfeq + nfuneq,escape = F,bold = F,latex_gap_space = "1ex")%>%
# #   group_rows("Unequal",2*nfeq + nfuneq +1,2*(nfeq + nfuneq),escape = F,bold = F,latex_gap_space = "1ex")
# # 
# # 
# # 
