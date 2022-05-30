rm(list = ls())
library(readxl)
library(reshape2)
library(ggplot2)
library(viridis)
library(sf)

############ LOAD DATA AND SHAPE SWITZERLAND

SwissCommune <- read_sf("C:/Users/jauslinr/switchdrive/Swiss2021/ggg_2021-LV95/shp/k4g21_01012021.shp")
SwissLake <- read_sf("C:/Users/jauslinr/switchdrive/Swiss2021/ggg_2021-LV95/shp/k4s21.shp")

SwissCommune <- read_sf("C:/Users/jauslinr/switchdrive/Swiss2021/ggg_2021-LV03/shp/k4g21_01012021.shp")
SwissCanton <- read_sf("C:/Users/jauslinr/switchdrive/Swiss2021/ggg_2021-LV03/shp/k4k21.shp")
SwissLake <- read_sf("C:/Users/jauslinr/switchdrive/Swiss2021/ggg_2021-LV03/shp/k4s21.shp")

amphib <- read_excel("C:/Users/jauslinr/switchdrive/Amphibians/amphib.xlsx")

########### REMOVE DUPLICATED SITE

# coord <- data.frame(x = amphib$X,
# y = amphib$Y)
# test <- coord[duplicated(coord),]
# amphib[which(amphib$X == test[1,1] & amphib$Y == test[1,2]),]
# amphib[which(amphib$X == test[2,1] & amphib$Y == test[2,2]),]
# amphib[which(amphib$X == test[3,1] & amphib$Y == test[3,2]),]
# amphib[which(amphib$X == test[4,1] & amphib$Y == test[4,2]),]
# amphib[which(amphib$X == test[5,1] & amphib$Y == test[5,2]),]
# amphib[which(amphib$X == test[6,1] & amphib$Y == test[6,2]),]
# amphib[which(amphib$X == test[7,1] & amphib$Y == test[7,2]),]
# amphib[which(amphib$X == test[8,1] & amphib$Y == test[8,2]),]
# amphib[which(amphib$X == test[9,1] & amphib$Y == test[9,2]),]
rm_site <- c(11149,28922,28742,6123,7224,11183,11184,11058,1772)
for(i in rm_site){
  amphib <- amphib[-which(amphib$ID1 == i),]  
}
# coord <- data.frame(x = amphib$X,
# y = amphib$Y)
# test <- coord[duplicated(coord),]

################ LAND SQUARES OF 1 KM2

amphib$xmin <- amphib$X - 500
amphib$xmax <- amphib$X + 500
amphib$ymin <- amphib$Y - 500
amphib$ymax <- amphib$Y + 500


species_name <- c("ALOB","BOVA","BUBU","BUCA","HYAR","HYIN","RATE","RALA","RADA","PERI","PEAG","TRAL","TRHE","TRVU","TRCR","TRCA","SASA")
M <- amphib[,species_name]
M[!is.na(M)] <- 1
rowSums(M,na.rm = TRUE)



amphib_long <- reshape2::melt(amphib,id.vars = c("CANTON","COMMUNE","OBNR","ID1","BED","NOM","X","Y","xmin","xmax","ymin","ymax"),
                              variable.name = "species",
                              value.name = "occ",
                              na.rm = TRUE)
amphib_long_2019 <- amphib_long[amphib_long$occ == '2019',]







ggplot() + 
  geom_sf(data = SwissCanton,color = "black",size = 0.2)+
  geom_sf(data = SwissLake,fill = "skyblue",color = "black",size = 0.2)+
  geom_rect(data = amphib_long_2019,
            aes(xmin = xmin,xmax = xmax,ymin = ymin,ymax = ymax,fill = species,colour = species),
            size = 0.4)+
  # geom_point(data = amphib_long,aes(x = X,y = Y, colour = species),size = 1)+
  scale_colour_viridis_d()




















