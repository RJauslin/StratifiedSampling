# install.packages("laeken")
rm(list = ls())
library(laeken)
library(sampling)
library(parallel)
library(BalancedSampling)
devtools::load_all(".")

# path <- "C:/Users/Raph/switchdrive/StratifiedSampling/StratifiedSampling"
path <- "C:/Users/Raphael/switchdrive/StratifiedSampling/StratifiedSampling"



devtools::load_all(path)
source(file.path(path,"./inst/matching/simu2.R"))

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
res_true <- tapply(YZ$income,list(YZ$ecostat),mean)


# library(ggplot2)
# library(ggridges)
# library(viridis)
# ggplot(YZ,aes(x = ecostat,y = income, fill = ecostat))+
  # scale_fill_viridis( option = "G",begin = 0.1,end = 0.7) +
  # scale_fill_grey( start = 0.4, end = 1)+
  # geom_violin(width = 0.9) +
  # geom_jitter(position=position_jitter(0.2))
  # stat_boxplot(geom ='errorbar',width = 0.5)+
  # geom_boxplot(color = "black",outlier.alpha = 0.5)
  
# ggplot(YZ,aes(x = income, y = ecostat,fill = ..x..))+
#   geom_density_ridges_gradient(scale = 1.5) +
#   theme_ridges() +
#   theme(legend.position = "none")+
#   scale_fill_viridis( option = "H",begin = 0.1,end = 0.7)
# 


N <- nrow(eusilc)
n1 <- 1000
n2 <- 4000
strata <- as.numeric(as.vector(eusilc$pl030))

sd_strata <- aggregate(x = YZ$income,by= list(YZ$ecostat),FUN = "sd")
size_strata <- aggregate(x = YZ$income,by= list(YZ$ecostat),FUN = "sum")

p <- sd_strata$x*size_strata$x

nh1 <- p/sum(p)*n1
nh2 <- p/sum(p)*n2
table(eusilc$pl030)

pik1 <- inclusionprobastrata(strata,nh1)
pik2 <- inclusionprobastrata(strata,nh2)

Xcat <- disj(strata)
t(Xcat)%*%pik1
t(Xcat)%*%pik2


s <- stratifiedcube(X = as.matrix(Xm),
                    strata = strata,
                    pik = pik1)
sum(s)

colSums(Xm)
t(Xm/pik1)%*%s

Xcat <- disj(strata)
  
t(Xcat)%*%s
t(Xcat)%*%pik1

SIM <- 1000

# set.seed(156)
# simu_strata(Xm,Y,Z,id,pik1,pik2,strata,totals = TRUE)

l_sim <- list()
for(i in 1:SIM){
  print(i)
  set.seed(i)
  # l_sim[[i]] <- simu2(Xm,Y,Z,id,n1,n2,totals = TRUE)
  l_sim[[i]] <- simu_strata(Xm,Y,Z,id,pik1,pik2,strata,totals = TRUE)
}


##------------- simu parrallel
# cl <- makeCluster(detectCores())
# clusterEvalQ(cl,{
#   path <- "C:/Users/Raph/switchdrive/StratifiedSampling/StratifiedSampling"
#   library(BalancedSampling)
#   devtools::load_all(path)
#   source(file.path(path,"./inst/matching/simu2.R"))
# })
# 
# f1 <- function(n,Xm,Y,Z,id,n1,n2){
#   l <- simu2(Xm,Y,Z,id,n1,n2,totals = FALSE)
#   return(list(res_ren =  l$res_ren,res_opt = l$res_opt, res_ran = l$res_ran))
# }
# 
# 
# # f1(1,Xm,Y,Z,id,n1,n2)
# 
# # Sys.sleep(1800)
# 
# 
# start <- Sys.time ()
# 
# l_sim <- parLapply(cl = cl,
#                    X = 1:SIM,
#                    fun = f1,
#                    Xm = Xm,
#                    Y = Y,
#                    Z = Z,
#                    id = id,
#                    n1 = n1,
#                    n2 = n2)
# print(Sys.time () - start)
# 
# stopCluster(cl)

#------ save

saveRDS(l_sim, file = "C:/Users/Raphael/switchdrive/StratifiedSampling/StratifiedSampling/inst/matching/l_sim3.rds")


#------ load


l_sim <- readRDS("C:/Users/Raphael/switchdrive/StratifiedSampling/StratifiedSampling/inst/matching/l_sim3.rds")
length(l_sim)


#----- processing


res_ren <- lapply(l_sim,function(x){x$res_ren})
res_opt <- lapply(l_sim,function(x){x$res_opt})
res_ran <- lapply(l_sim,function(x){x$res_ran})


m_opt <- colMeans(do.call(rbind,res_opt),na.rm =TRUE)
m_ran <- colMeans(do.call(rbind,res_ran),na.rm= TRUE)
m_ren <- colMeans(do.call(rbind,res_ren),na.rm = TRUE)


mse_opt <- colMeans(do.call(rbind,lapply(res_opt, function(x){(x - res_true)^2})),na.rm = TRUE)
mse_ren <- colMeans(do.call(rbind,lapply(res_ren, function(x){(x - res_true)^2})),na.rm = TRUE)
mse_ran <- colMeans(do.call(rbind,lapply(res_ran, function(x){(x - res_true)^2})),na.rm = TRUE)


b2_opt <- (m_opt - res_true)^2
b2_ren <- (m_ren - res_true)^2
b2_ran <- (m_ran - res_true)^2


v_opt <- colMeans(do.call(rbind,lapply(res_opt, function(x){(x - m_opt)^2})),na.rm = TRUE)
v_ren <- colMeans(do.call(rbind,lapply(res_ren, function(x){(x - m_ren)^2})),na.rm = TRUE)
v_ran <- colMeans(do.call(rbind,lapply(res_ran, function(x){(x - m_ran)^2})),na.rm = TRUE)




mse_opt
b2_opt + v_opt





tab <- rbind(mse_opt,
             mse_ran,
             mse_ren,
             b2_opt,
             b2_ran,
             b2_ren,
             v_opt,
             v_ran,
             v_ren) 

tab/1e6

sum(mse_opt)/1e6
sum(mse_ran)/1e6
sum(mse_ren)/1e6


