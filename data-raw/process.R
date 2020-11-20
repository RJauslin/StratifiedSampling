library(sf)
# library(ggplot2)
# pathIni <- getwd()
# 
# setwd("C:/Users/jauslinr/switchdrive/StratifiedSampling/StratifiedSampling/data-raw")

SwissCommune_LV95 <- read_sf("k4g19_LV95.shp")
SwissCanton_LV95 <- read_sf("k4k19_LV95.shp")
SwissRegion_LV95 <- read_sf("k4r19_LV95.shp")
Swiss_LV95 <- read_sf("k4l19_LV95.shp")
SwissLake_LV95 <- read_sf("k4s19_LV95.shp")


SwissCommune_LV03 <- read_sf("k4g19_LV03.shp")
SwissCanton_LV03 <- read_sf("k4k19_LV03.shp")
SwissRegion_LV03 <- read_sf("k4r19_LV03.shp")
Swiss_LV03 <- read_sf("k4l19_LV03.shp")
SwissLake_LV03 <- read_sf("k4s19_LV03.shp")




# ggplot()+
#   geom_sf(data = Swiss_LV95,fill = "transparent",color = "black",size = 0.5)+
#   geom_sf(data = SwissCommune_LV95,fill = "transparent",color = "black",size = 0.5)+
#   geom_sf(data = SwissRegion_LV95,fill = "transparent",color = "grey",size = 0.5)+
#   geom_sf(data = SwissLake_LV95,fill = "grey50")
#   # geom_sf(data = Swiss_LV03,fill = "transparent",color = "black",size = 0.5)+
#   # geom_sf(data = SwissCommune_LV03,fill = "transparent",color = "black",size = 0.5)+
#   # geom_sf(data = SwissRegion_LV03,fill = "transparent",color = "grey",size = 0.5)+
#   # geom_sf(data = SwissLake_LV03,fill = "grey50")

# setwd(pathIni)


# This should be the last line.
# Note that names are unquoted.
# I like using overwrite = T so everytime I run the script the 
# updated objects are saved, but the default is overwrite = F
usethis::use_data(SwissCommune_LV95,
                  SwissCommune_LV03,
                  SwissCanton_LV95,
                  SwissCanton_LV03,
                  SwissRegion_LV95,
                  SwissRegion_LV03,
                  Swiss_LV95,
                  Swiss_LV03,
                  SwissLake_LV95,
                  SwissLake_LV03,
                  overwrite = T)

