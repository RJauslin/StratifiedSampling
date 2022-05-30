
library(Rcpp)
library(WaveSampling)
library(sampling)
library(BalancedSampling)
library(spatstat)
library(spsurvey)
library(sp)
library(SDraw)

library(knitr)
library(readr)
library(tikzDevice)

library(raster)
library(rgdal)
library(sf)
library(rgeos)
library(gridExtra)
library(ggvoronoi)
library(grid)
library(lattice)
library(ggplot2)
library(ggrepel)
library(scico)

library(kableExtra)
library(magrittr)
library(plyr)

library(gganimate)
library(gifski)

rm(list = ls())

set.seed(3)

gauss2d <- function(x,y,x0,y0,a,b,c,A){
  return(A*exp(-(a*(x-x0)^2 + 2*b*(x-x0)*(y-y0) + c*(y-y0)^2)))
}


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

n <- 10

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


# s <- balseq(pik_pp$uneq,Xaux_pp$uneq,Xspread = X_pp)

pik <- pik_pp$eq
Xaux <- Xaux_pp$eq
Xspread <- X_pp

############################################################################################

deg = 1
N <- length(pik)
eps <- 1e-6

pikInit <- pik
index <- which(pik > eps & pik < (1-eps))

n <- 0
counter <- 1


dat_index2n <- data.frame(x = NULL,
           y = NULL,
           type = NULL)
dat_center <- data.frame(x = NULL,
                         y = NULL,
                         type = NULL)
dat_unit1 <- data.frame(x = NULL,
                        y = NULL,
                        type = NULL)
dat_unit0 <- data.frame(x = NULL,
                        y = NULL,
                        type = NULL)


#----------- MAIN LOOP
while(length(index) > 0){
  # cat("Step :",counter,"\n")
  # print(length(index))
  # print(n)
  
  i <- which.max(pik[index])
  
  if(!is.null(Xspread)){
    #take distance of the considered unit 
    d <- distUnitk(Xspread,index[i],F,F)
    # modify index respect to distance
    index <- index[order(d[index])]
  }else{
    
    # tmp <- index[1]
    # index[1] <- index[i]
    # index[i] <- tmp
    
    index <- index[order(pik[index],decreasing = TRUE)]
  }
  
  l <- balseq_onestep(Xaux,pik,pikInit,index,deg)
  status <- l$status
  v = l$v
  n <- l$n
  
  unit0 <- which(pik < eps)
  unit1 <- which(pik > 1- eps)
  
  
  
  # tmp <- Xspread[index[1:n],]
  # tmp <- tmp[chull(tmp),]
  
  # plot(Xspread)
  # lines(Xspread[index[1:n],1],Xspread[index[1:n],2],type ="p",pch = 16,col = "cyan")
  # # lines(tmp[,1],tmp[,2],type ="p",pch = 16,col = "blue")
  # 
  # lines(Xspread[index[1],1],Xspread[index[1],2],type ="p",pch = 16,col = "red")
  # lines(Xspread[unit0,1],Xspread[unit0,2],type ="p",pch = 16)
  # lines(Xspread[unit1,1],Xspread[unit1,2],type ="p",pch = 16,col = "orange")
  # 
  
  dat_index2n <- rbind(dat_index2n,data.frame(x = Xspread[index[2:n],1],
                            y = Xspread[index[2:n],2],
                            type = rep(counter,length(index[2:n]))))
  dat_center <- rbind(dat_center,data.frame(x = Xspread[index[1],1],
                           y = Xspread[index[1],2],
                           type = rep(counter,length(index[1]))))
  dat_unit1 <- rbind(dat_unit1,data.frame(x = Xspread[unit1,1],
                          y = Xspread[unit1,2],
                          type = rep(counter,length(unit1))))
  dat_unit0 <- rbind(dat_unit0,data.frame(x = Xspread[unit0,1],
                          y = Xspread[unit0,2],
                          type = rep(counter,length(unit0))))
  
  # p <- ggplot()+
  #   geom_point(data = data.frame(x = Xspread[,1],y = Xspread[,2]),aes(x = x,y = y),shape = 1,size = 2)+
  #   geom_point(data = data.frame(x = Xspread[index[2:n],1],y = Xspread[index[2:n],2]),aes(x = x,y = y),shape = 16,size = 2,color = "cyan") +
  #   geom_point(data = data.frame(x = Xspread[index[1],1],y = Xspread[index[1],2]),aes(x = x,y = y),shape = 16,size = 2,color = "red")+
  #   geom_point(data = data.frame(x = Xspread[unit0,1],y = Xspread[unit0,2]),aes(x = x,y = y),shape = 1,size = 2,color = "white") +
  #   geom_point(data = data.frame(x = Xspread[unit1,1],y = Xspread[unit1,2]),aes(x = x,y = y),shape = 16,size = 2) + 
  #   theme_minimal() +
  #   theme(
  #     text = element_text(family="sans",color = "black",size = 9),
  #     panel.spacing = unit(2, "lines"),
  #     # title
  #     plot.title = element_text(hjust = 0.5,size = 9),
  #     # axes
  #     axis.line=element_blank(),
  #     axis.ticks=element_blank(),
  #     # legend
  #     legend.position="bottom",
  #     legend.title = element_text(size = 9,vjust = +1.0),
  #     legend.key.size = unit(0.4, "cm"),
  #     legend.key.width = unit(1,"cm") ,
  #     # background colors
  #     panel.background = element_rect(fill = "grey70",
  #                                     colour = "grey70"),
  #     # panel.grid.major = element_line(size = 0.5, linetype = 'solid',
  #                                     # colour = "white"), 
  #     # panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
  #                                     # colour = "white"),
  #     # panel.border=element_rect(colour = "black",fill = "transparent"),
  #     panel.border = element_blank(),
  #     panel.grid.major=element_blank(),
  #     panel.grid.minor=element_blank(),
  #     # keep edge black facet_wrap
  #     # strip.background = element_rect(fill="white"),
  #     strip.text =element_text(color = "black",size = 8)
  #   )
  # 
  #   # geom_polygon(data = data.frame(x = tmp[,1],y = tmp[,2]),aes(x = x,y = y),alpha = 0.4)+
  #   # geom_point(data = data.frame(x = Xspread[-unit0,1],y = Xspread[-unit0,2]),aes(x = x,y = y),shape = 1,size = 2) +
  #   # geom_point(data = data.frame(x = Xspread[unit1,1],y = Xspread[unit1,2]),aes(x = x,y = y),shape = 16,size = 2) +
  #   # geom_point(data = data.frame(x = Xspread[index[1],1], y = Xspread[index[1],2]),aes(x = x,y = y),color = "red")
  # print(p)
  # #   
  # # 
  # #   ggplot(data = data.frame(x = tmp[,1],y = tmp[,2]),aes(x = x,y = y))+
  # #     geom_polygon(alpha = 0.4)+
  # #     geom_shape(radius = unit(3, 'cm'))
  # #   
  # 
  # Sys.sleep(1)
  
  # if we can no longer find solution and index is at the end of the vector then exit and return pikstar
  if(status == 1 & is.na(index[n])){
    # return(pik)
    break;
  }else{
    
    v <-  v - pmin(pik[index[2:n]],(1-pik[index[2:n]])*pik[index[1]]/(1-pik[index[1]]))
    
    if(stats::runif(1) < pik[index[1]]){
      pik[index[2:n]] <- pik[index[2:n]] - v*(1-pik[index[1]])/pik[index[1]]
      pik[index[1]] <- 1
    }else{
      pik[index[2:n]] <- pik[index[2:n]] + v
      pik[index[1]] <- 0
    }
    
    index <- which(pik > eps & pik < (1-eps))
    
  }
  counter <- counter + 1
}


# dat_index2n$type



# ntype = 40

dat_index2n <- dat_index2n[-nrow(dat_index2n),]

any(is.na(dat_index2n))
any(is.na(dat_center))
any(is.na(dat_unit0))
any(is.na(dat_unit1))

# dat_index2n <- dat_index2n[which(dat_index2n$type < ntype),]
# dat_center <- dat_center[which(dat_center$type < ntype),]
# dat_unit0 <- dat_unit0[which(dat_unit0$type < ntype),]
# dat_unit1 <- dat_unit1[which(dat_unit1$type < ntype),]

col1 <- "#1c658c"
col2 <- "#d8d2cb"

p <- ggplot()+
    geom_point(data = data.frame(x = Xspread[,1],y = Xspread[,2]),aes(x = x,y = y),shape = 1,size = 2,color = "grey30")+
    geom_point(data = dat_index2n,aes(x = x,y = y,group = type),shape = 16,size = 2,color = col1)+
    geom_point(data = dat_center,aes(x = x,y = y,group = type),shape = 16,size = 2,color = "red") +
    geom_point(data = dat_unit0,aes(x = x,y = y,group = type),shape = 1,size = 2,color = "white") +
    geom_point(data = dat_unit1,aes(x = x,y = y,group = type),shape = 16,size = 2)+ 
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
    legend.key.size = unit(0.4, "cm"),
    legend.key.width = unit(1,"cm") ,
    # background colors
    panel.background = element_rect(fill = col2,
                                    colour = col2),
    # panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    # colour = "white"),
    # panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    # colour = "white"),
    # panel.border=element_rect(colour = "black",fill = "transparent"),
    panel.border = element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    # keep edge black facet_wrap
    # strip.background = element_rect(fill="white"),
    strip.text =element_text(color = "black",size = 8)
  )
p

anim1 <- p + transition_states(states = type)
anim1
# gganimate::animate(anim1, nframes =  6*max(dat_center$type))

a <- gganimate::animate(anim1, nframes = 6*max(dat_center$type),
                        renderer = file_renderer("C:/Users/jauslinr/switchdrive/Perrugia"))


# 
# gganimate::animate(anim1, nframes = 3*max(dat_center$type),width = 1200, height = 1000,
#                    renderer = gifski_renderer("meuse.gif"))
