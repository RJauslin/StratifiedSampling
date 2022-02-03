
library(sf)
library(ggplot2)
library(sampling)
library(viridis)
library("RColorBrewer")  
library(BalancedSampling)

rm(list = ls())
data("belgianmunicipalities")
pik <- inclusionprobabilities(belgianmunicipalities$Tot04,200)
# index <- order(belgianmunicipalities$Tot04) 

# pik <- pik[index]

# index <- which(pik < (1-1e-8))
# pik <- pik[index]

N <- length(pik)
n <- sum(pik)


y=belgianmunicipalities$TaxableIncome
# y=belgianmunicipalities$averageincome
sim = 10000
nom <- c("Osod","Random systematic","Random pivotal","Tille","Midzuno","Max entropy")
ss=array(0,c(sim,length(nom)))
colnames(ss) <- nom
ht=numeric(length(nom))



l <- list( osod,
           UPrandomsystematic,
           rpm,
           UPtille,
           UPmidzuno,
           maxent)

for(i in 1:sim){
  cat("Step ",i,"\n")
  for(k in 1:length(l)){
    # if(k == 1){
      o <- sample(1:length(y))
      y_tmp <- y[o]
      pik_tmp <- pik[o]
      # s <- l[[k]](pik_tmp)
      # ht[k] <- HTestimator(y_tmp[s==1],pik_tmp[s==1])
    # }else{
      s <- l[[k]](pik_tmp)
      if(nom[k] == "Random pivotal"){
        s_i <- s
        s <- rep(0,length(pik))
        s[s_i] <- 1
      }
      ht[k] <- HTestimator(y_tmp[s==1],pik_tmp[s==1])
    }
  # }
  ss[i,]=ht
}


ss2 <- melt(ss)
saveRDS(ss2,file = "C:/Users/jauslinr/switchdrive/flow_project/1- One Step One Decision Stream Sampling/Rnw/ss2.rds")
ss2 <- readRDS("C:/Users/jauslinr/switchdrive/flow_project/1- One Step One Decision Stream Sampling/Rnw/ss2.rds")

bp <- ggplot(data=ss2, aes(x = Var2, y=value)) +
  geom_hline(yintercept = sum(y),lty = 2)+
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(width = 0.6, fill = "lightgrey") 

# geom_boxplot()+
# theme_perso()
bp
