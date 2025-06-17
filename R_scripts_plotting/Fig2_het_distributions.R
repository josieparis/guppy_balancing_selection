# new general run
rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

lib<-c("dplyr","ggplot2", "cowplot", "ggrepel", "gridExtra","tidyverse")
lapply(lib,library,character.only=T)

setwd("<dir>/Analyses/NFDS_analysis/popgen_stats/")

data <- read.csv("data/ALL_het_50k_winds.txt",h=T,sep=" ")  

GHP_data <- data  %>% dplyr::filter(comp=="GHP")

data_filtered <- data %>%
  mutate(mean_het = ifelse(is.na(mean_het), 0, mean_het))  # Convert NA to 0

# Create the histogram plot
GHP_plot <- ggplot(data_filtered, aes(x = mean_het)) +
  geom_histogram(binwidth = 0.01, fill = "black", color = "black") +  # Histogram bars
  theme_bw() +
  labs(x = "Het/Kb", y = "# windows") +
  scale_x_continuous(breaks=seq(0,0.6,0.2),limits=c(-0.01,0.6))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=6,family="Times"),
        axis.text.x=element_text(size=6,family="Times"),
        axis.title=element_text(size=16,family="Times"))

ggsave("../figs/MLP_theta_windows.png", GHP_plot, width = 4.2, height = 4, units="cm",dpi = 400)
