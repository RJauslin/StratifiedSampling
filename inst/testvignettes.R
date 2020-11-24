
rm(list = ls())

# findIndex <- function(x1,X){
#   final <- c()
#   for(i in 1:nrow(x1)){
#     id <- seq(1,nrow(X),1)
#     for(j in 1:ncol(x1)){
#       out = which(x1[i,j] == X[id,j])
#       if(length(out) > 1){
#         id <- out
#       }else if(length(out) == 1 & j == 1){
#         final <- c(final,out)
#         break;
#       }else if(length(out) == 1){
#         final <- c(final,id[out])
#         break;
#       }
#     }
#   }
#   return(final)
# }

library(sf)
library(ggplot2)

library(sampling)

data("SwissCanton_LV95")
# 
# ggplot()+
#   # geom_raster(data = relief,aes(x = x,y = y,alpha = SwissRelief))+
#   geom_sf(data = SwissCommune,fill = "transparent",color = "black",size = 0.5)+
#   # geom_sf(data = SwissRegion,fill = "transparent",color = "grey",size = 0.5)+
#   # geom_sf(data = SwissLake,fill = "grey50")+
#   scale_alpha(name = "",
#               range = c(0.6, 0),
#               guide = F)
# 
# 
# 
# 
# ent1 <- readRDS("./inst/article/ent3.rds")
# 
# 
# 
# dat <- ent1
# 
# cat <- readRDS("./inst/article/cat3.rds")
# dat <- dat[-which(cat == 0),]
# cat <- cat[-which(cat == 0)]
# cat <- cleanstrata(cat)
# 
# dat$cat <- cat
# dat <- dat[order(dat$cat),]
# dat$cat <- cleanstrata(dat$cat)
# 
# 
# # Extract data
# 
# N <- nrow(dat)
# n <- 1000
# Nh <- table(dat$cat)
# 
# nh <- ceiling(n*Nh/N)
# length(which(nh < 1))
# 
# pik <- as.vector(inclusionprobastrata(dat$cat,nh))
# X <- as.matrix(as.vector(pik))
# nrow(X)
# length(pik)
# strata <- dat$cat
# s <- stratifiedcube(as.matrix(pik),dat$cat,pik)
# 
# 
# 
# for(i in 1:length(unique(strata))){
#   if(round(sum(pik[strata == i]),7) != round(sum(s[strata == i]),7)){
#     print(sum(pik[strata == i]))
#     print(sum(s[strata == i]))
#     print("ERROR not ok")
#   }
# }
# 
# 
# establishement <- dat[s== 1,]
# 
# 
# 
# 
# 
# 
# 
# 
# Kanton <- rep(0,nrow(establishement))
# dat <- data.frame(x = establishement$E_KOORD,
#                   y = establishement$N_KOORD,
#                   z = establishement$B08EMPT)
# for( i in 1:nrow(SwissCanton_LV95)){
#   print(i)
# 
#   dat1 <- rasterFromXYZ(dat)
#                         # crs= CRS("+proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333 +k_0=1 +x_0=2600000 +y_0=1200000 +ellps=bessel +towgs84=674.374,15.056,405.346,0,0,0,0 +units=m +no_defs"))
#                         # crs = CRS("+proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333 +k_0=1 +x_0=600000 +y_0=200000 +ellps=bessel +units=m +no_defs"))
#   dat1 <- mask(dat1,SwissCanton_LV95[i,])
#   # dat1 <- crop(dat1,SwissCommune[i,])
#   # print(length(dat1$z@data@values))
#   test <- coordinates(dat1)[which(!is.na(dat1@data@values)),]
#   if(length(test) != 0){
#     if(length(test) == 2){
#       test <- matrix(test,ncol= 2)
#     }else{
#       test <- as.matrix(test)
#     }
#     index <- findIndex(test,dat[,1:2])
#     if(Kanton[index] != 0){
#       break;
#     }
#     Kanton[findIndex(as.matrix(test),dat[,1:2])] <- i
#     as(dat1,"SpatialPixels")
#   }
# 
# }
# 
# establishement$Kanton <- Kanton
# 
# saveRDS(establishement,"./inst/article/establishement.rds")
# saveRDS(establishement,"establishement.rds")





N <- nrow(establishement)
n <- 260
Nh <- table(establishement$Kanton)

nh <- ceiling(n*Nh/N)
nh <- n*Nh/N
# length(which(nh < 1))

pik <- as.vector(inclusionprobastrata(establishement$Kanton,nh))
X <- as.matrix(as.vector(pik))
nrow(X)
length(pik)
strata <- establishement$Kanton
s <- stratifiedcube(as.matrix(pik),strata,pik)

t(X/pik)%*%pik
t(X/pik)%*%s

Xcat <- disj(as.matrix(strata))


t(Xcat)%*%pik
t(Xcat)%*%s

library(viridis)
p_neuch <- ggplot()+
  geom_sf(data = SwissLake_LV95,fill = "lightskyblue",color = "grey40",size = 0.5)+
  geom_sf(data = SwissCanton_LV95,aes(fill = KtName),color = "grey40",size = 0.5)+
  geom_point(data = establishement,
             aes(x=E_KOORD,y = N_KOORD),
             shape = 1,
             size = 1,
             colour = "black")+
  geom_point(data = establishement[s == 1,],
             aes(x=E_KOORD,y = N_KOORD),
             shape = 16,
             colour = "black",
             size = 1)+
  # scale_alpha(name = "",
  #             range = c(0.6, 0),
  #             guide = F)+
  scale_fill_viridis_d()+
  labs(x = NULL,
       y = NULL,
       title = "Neuch\\^atel small medium-sized enterprise",
       # title = "",
       # subtitle = "STATENT 2016",
       size = "Employee",
       caption = NULL) 
  # scale_size(range = c(0.5, 3))
p_neuch

