# source load data

source("./inst/balseq/amphibians/data_loading.R")


# ggplot 

p1 <- ggplot() + 
  geom_sf(data = SwissBioRegion,aes(fill = as.factor(FRRegionNa)),color = "black",size = 0.01)+
  # geom_sf(data = SwissCanton,color = "black",size = 0.2,fill = "grey95")+
  geom_sf(data = SwissLake,fill = "white",color = "black",size = 0.2)+
  # geom_rect(data = amphib,
  # aes(xmin = xmin,xmax = xmax,ymin = ymin,ymax = ymax),
  # size = 0.4)+
  # geom_point(data = amphib,aes(x = COORD_X,y = COORD_Y,size = AREA),shape = 1)+
  # geom_point(data = amphib,aes(x = COORD_X,y = COORD_Y,size = AREA),colour = "grey10",shape = 1,stroke = 0.1)+
  # scale_colour_viridis_c(option = "G",end = 0.9,begin = 0.1)+
  scale_fill_viridis_d(option = "G")+
  # scale_fill_grey("Bio Region",start = 0,end = 1)+
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
    panel.background=element_blank(),
    panel.border=element_rect(colour = "black",fill = "transparent"),
    # panel.grid.major=element_blank(),
    # panel.grid.minor=element_blank(),
    # keep edge black facet_wrap
    # strip.background = element_rect(fill="white"),
    strip.text =element_text(color = "black",size = 8)
  )


p2 <- ggplot() + 
  geom_sf(data = SwissBioRegion,fill = "white",color = "black",size = 0.01)+
  geom_sf(data = SwissLake,fill = "white",color = "black",size = 0.2)+
  geom_point(data = amphib,aes(x = COORD_X,y = COORD_Y,size = AREA,colour = diversity),stroke = 0.1)+
  geom_point(data = amphib,aes(x = COORD_X,y = COORD_Y,size = AREA),shape = 1,stroke = 0.1)+
  # scale_colour_viridis_c(option = "G")+
  scale_colour_gradient(low = "grey10",high = "grey90")+
  scale_fill_grey("Bio Region",start = 0.2,end = 0.9)+
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
    panel.background=element_blank(),
    panel.border=element_rect(colour = "black",fill = "transparent"),
    # panel.grid.major=element_blank(),
    # panel.grid.minor=element_blank(),
    # keep edge black facet_wrap
    # strip.background = element_rect(fill="white"),
    strip.text =element_text(color = "black",size = 8)
  )

tikz(file = "stratmat.tex", width =  5.78851, height = 3,standAlone = FALSE)
gridExtra::grid.arrange(p1, p2, nrow = 1)
dev.off()

