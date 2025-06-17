# new general run
rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

# Plot Trinidad and rivers
lib<-c("cowplot","ggplot2","data.table","viridis","patchwork","sf","ggpubr","raster")
lapply(lib,library,character.only=T)

setwd("~/Dropbox/Sussex_Guppies/Analyses/NFDS_analysis/map/")

# Read in the shape file
trinidad <- st_read(dsn = "tto_adm0/TTO_adm0.shp")

## Shape file downloaded from: https://data.humdata.org/dataset/cod-ab-tto?
trinidad_full_map <- ggplot(data = trinidad) +
  geom_sf(color = "black", lwd=0.5)+
  geom_sf(fill = "#848884", alpha=0.1) +
  xlim(c(-61.9, -60.9)) + 
  ylim(c(10.0, 10.85))+
  theme_void() 

ggsave(trinidad_full_map, width = 12, height = 12, units = "cm", dpi = 400,
       file = "../figs/trinidad_full_map.pdf",device=cairo_pdf)

## Now for a zoom in map with the rivers included and points plotted

### Two shape files tried:
### Shape file 1 
all_rivers <- st_read(dsn = "Rivers_(Trinidad)_441799794901450010/Rivers_(Trinidad).shp")

## rivers downloaded from: https://data-ttacis.hub.arcgis.com/datasets/c3c857dbb3e6481f9452c5f95a792a42_0/explore 

### Shape file 2
all_rivers2 <- st_read(dsn = "Trinidat_Tobago_SRTM_Streamlines/Export_Output_2.shp")
all_rivers2 <- st_zm(all_rivers2)

# from: https://www.researchgate.net/post/Can-anyone-help-me-with-a-GIS-map-of-Trinidad-displaying-the-rivers 

## new lat long coords
study_map <- ggplot() +
  geom_sf(data = trinidad, color = "black", lwd=0.5,fill = "#848884", alpha=0.1) +
  geom_sf(data = all_rivers2,colour="#6F8FAF")+
  geom_point(aes(x=-61.25,y=10.65),fill="#3B2778",colour="black",size=2.5, pch=21, stroke=0.7)+  ### GHP POINT
  geom_point(aes(x=-61.26,y=10.71),fill="#9772C3",colour="black",size=2.5, pch=24,stroke=0.7)+  ### GLP POINT
  geom_point(aes(x=-61.23,y=10.65),fill="#5F7F3C",colour="black",size=2.5, pch=21,stroke=0.7)+  ### APHP POINT
  geom_point(aes(x=-61.23,y=10.67),fill="#C6D573",colour="black",size=2.5, pch=24,stroke=0.7)+  ### APLP POINT
  geom_point(aes(x=-61.31,y=10.77),fill="#9D2919",colour="black",size=2.5, pch=21,stroke=0.7)+  ### MHP POINT
  geom_point(aes(x=-61.29,y=10.75),fill="#BE261B",colour="black",size=2.5, pch=24,stroke=0.7)+  ### MLP POINT
  geom_point(aes(x=-61.27,y=10.75),fill="#CB6078",colour="black",size=2.5, pch=24,stroke=0.7)+  ### P POINT
  geom_point(aes(x=-61.17,y=10.66),fill="#1D4A8A",colour="black",size=2.5, pch=21,stroke=0.7)+  ### TUHP POINT
  geom_point(aes(x=-61.17,y=10.68),fill="#5F84B2",colour="black",size=2.5, pch=24,stroke=0.7)+  ### TULP POINT
  geom_point(aes(x=-61.2676,y=10.6599),fill="#F4B740",colour="black",size=2.5, pch=21,stroke=0.7)+  ### ECHP POINT
  geom_point(aes(x=-61.2658,y=10.6636),fill="#D7AD7F",colour="black",size=2.5, pch=24,stroke=0.7)+  ### ECLP POINT
  xlim(c(-61.85, -60.9)) +  # Se
  ylim(c(10.59, 10.85)) +
  theme_void()

ggsave(study_map, width = 12, height = 12, units = "cm", dpi = 400,
       file = "../figs/rivers_map.pdf",device=cairo_pdf)
