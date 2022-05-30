
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

plot(X_ns)


dat_ggplot <- rbind(data.frame(x = X_ns[,1],y = X_ns[,2],type = rep("Neyman-Scott process",nrow(X_ns))),
                    data.frame(x = X_pp[,1],y = X_pp[,2],type = rep("Complete Spatial Randomness",nrow(X_pp))))

library(ggplot2)
p_artificial <- ggplot() + 
  geom_point(data = dat_ggplot, aes(x = x,y = y)) +
  facet_wrap(~type)+
  theme_minimal() + 
  theme(
      text = element_text(family="sans",color = "black",size = 9),
      panel.spacing = unit(2, "lines"),
      # title
      plot.title = element_text(hjust = 0.5,size = 9),
      # axes
      axis.line=element_blank(),
      axis.ticks=element_blank(),
      # legend
      legend.position="bottom",
      legend.title = element_text(size = 9,vjust = +1.0),
      legend.key.size = unit(0.3, "cm"),
      legend.key.width = unit(0.7,"cm") ,
      # background colors
      panel.background=element_blank(),
      panel.border=element_rect(colour = "black",fill = "transparent"),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      # keep edge black facet_wrap
      # strip.background = element_rect(fill="white"),
      strip.text =element_text(color = "black",size = 8)
    )
print(p_artificial)
