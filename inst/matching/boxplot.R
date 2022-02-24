rm(list = ls())
library(laeken)
library(sampling)
library(parallel)

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

theme_wave <- function(...) {
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
      legend.position= "none",
      # legend.title = element_text(size = 9,vjust = +1.0),
      # legend.key.size = unit(0.3, "cm"),
      # legend.key.width = unit(0.7,"cm") ,
      # background colors
      panel.background=element_blank(),
      panel.border=element_rect(colour = "black",fill = "transparent"),
      # panel.grid.major=element_blank(),
      # panel.grid.minor=element_blank(),
      # keep edge black facet_wrap
      # strip.background = element_rect(fill="white"),
      strip.text =element_text(color = "black",size = 8)
    )
}

library(ggplot2)
# library(ggridges)
library(viridis)
ggplot(YZ,aes(x = ecostat,y = income, fill = ecostat))+
  # scale_fill_viridis( option = "G",begin = 0.1,end = 0.7) +
  scale_fill_grey("Economical Status", start = 0.4, end = 1)+
  # geom_violin(width = 0.9) +
  # geom_jitter(position=position_jitter(0.2))
  stat_boxplot(geom ='errorbar',width = 0.3)+
  geom_boxplot(color = "black",outlier.alpha = 0.5,outlier.size = 0.8) +  
  labs(x=" Economical Status",y = " Income")+
  theme_wave()



 