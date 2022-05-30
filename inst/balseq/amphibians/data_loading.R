

########## REMOVE PACKAGE

rm(list = ls())

########## LOAD LIBRARY

library(readxl)
library(reshape2)
library(ggplot2)
library(viridis)
library(sf)
library(sampling)
library(Spbsampling)
library(BalancedSampling)
library(WaveSampling)
library(parallel)
library(SDraw)
devtools::load_all("C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling")

############ LOAD DATA AND SHAPE SWITZERLAND


Swiss <- read_sf("C:/Users/jauslinr/switchdrive/Swiss2021/ggg_2021-LV95/shp/k4l21.shp")
SwissCommune <- read_sf("C:/Users/jauslinr/switchdrive/Swiss2021/ggg_2021-LV95/shp/k4g21_01012021.shp")
SwissCanton <- read_sf("C:/Users/jauslinr/switchdrive/Swiss2021/ggg_2021-LV95/shp/k4k21.shp")
SwissLake <- read_sf("C:/Users/jauslinr/switchdrive/Swiss2021/ggg_2021-LV95/shp/k4s21.shp")
SwissBioRegion <- read_sf("C:/Users/jauslinr/switchdrive/Swiss2021/BiogeographischeRegionen/N2020_Revision_BiogeoRegion.shp")

# SwissCommune <- read_sf("C:/Users/jauslinr/switchdrive/Swiss2021/ggg_2021-LV03/shp/k4g21_01012021.shp")
# SwissCanton <- read_sf("C:/Users/jauslinr/switchdrive/Swiss2021/ggg_2021-LV03/shp/k4k21.shp")
# SwissLake <- read_sf("C:/Users/jauslinr/switchdrive/Swiss2021/ggg_2021-LV03/shp/k4s21.shp")

amphib <- read_excel("C:/Users/jauslinr/switchdrive/Amphibians/AM_OBJ_DATA_171024.xlsx")
amphib <- amphib[which(amphib$AREA > 1e-7),]



SwissBioRegion$FRRegionNa[SwissBioRegion$FRRegionNa == "Jura"] = "Jura"
SwissBioRegion$FRRegionNa[SwissBioRegion$FRRegionNa == "Plateau"] = "Central Plateau"
SwissBioRegion$FRRegionNa[SwissBioRegion$FRRegionNa == "Versant nord des Alpes"] = "North side of the Alps"
SwissBioRegion$FRRegionNa[SwissBioRegion$FRRegionNa == "Versant sud des Alpes"] = "South side of the Alps"
SwissBioRegion$FRRegionNa[SwissBioRegion$FRRegionNa == "Alpes centrales occidentales"] = "Western Central Alps"
SwissBioRegion$FRRegionNa[SwissBioRegion$FRRegionNa == "Alpes centrales orientales"] = "Central Eastern Alps"


species_name <- as.vector(do.call(rbind,strsplit(colnames(amphib)[grepl('ESP',colnames(amphib))],'_'))[,2])

colnames(amphib)[grepl('ESP',colnames(amphib))]  <- species_name
colnames(amphib)

M <- amphib[,species_name]
# any(apply(M,MARGIN = 1,FUN = function(x){all(is.na(M))})) # FALSE no empty sites
M1 <- 1/M
M1[is.na(M1)] <- 0


y1  <- apply(M,MARGIN = 1, FUN = function(x){
  tmp <- x[!is.na(x)]
  return(length(which(tmp == 1)))
})
y2  <- apply(M,MARGIN = 1, FUN = function(x){
  tmp <- x[!is.na(x)]
  return(length(which(tmp == 2)))
})
y3  <- apply(M,MARGIN = 1, FUN = function(x){
  tmp <- x[!is.na(x)]
  return(length(which(tmp == 3)))
})
y4  <- apply(M,MARGIN = 1, FUN = function(x){
  tmp <- x[!is.na(x)]
  return(length(which(tmp == 4)))
})

amphib$y1 <- y1 # nbr species of type 1
amphib$y2 <- y2 # nbr species of type 2
amphib$y3 <- y3 # nbr species of type 3
amphib$y4 <- y4 # nbr species of type 4


amphib$diversity <- rowSums(M1)
# amphib$diversity <- y1+y2+y3+y4 # sum of occurence of all species, it means that the more you have some species at a certain site, the more you have richness... ?


# amphib$xmin <- amphib$COORD_X - 500
# amphib$xmax <- amphib$COORD_X + 500
# amphib$ymin <- amphib$COORD_Y - 500
# amphib$ymax <- amphib$COORD_Y + 500
 