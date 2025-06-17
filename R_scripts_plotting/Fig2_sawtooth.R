# new general run
rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

lib<-c("dplyr","ggplot2", "cowplot", "ggrepel", "gridExtra","tidyverse")
lapply(lib,library,character.only=T)

## working dir
setwd("<dir>/Analyses/NFDS_analysis/popgen_stats/")

# Read in your data
data <- read.csv("data/ALL_het_50k_winds.txt",h=T,sep=" ")  

GHP_data <- data %>% dplyr::filter(comp=="GHP")
GLP_data <- data %>% dplyr::filter(comp=="GLP")
APHP_data <- data %>% dplyr::filter(comp=="APHP")
APLP_data <- data %>% dplyr::filter(comp=="APLP")
ECHP_data <- data %>% dplyr::filter(comp=="ECHP")
ECLP_data <- data %>% dplyr::filter(comp=="ECLP")
TUHP_data <- data %>% dplyr::filter(comp=="TUHP")
TULP_data <- data %>% dplyr::filter(comp=="TULP")
MHP_data <- data %>% dplyr::filter(comp=="MHP")
MLP_data <- data %>% dplyr::filter(comp=="MLP")
PLP_data <- data %>% dplyr::filter(comp=="P")

## colour codes:
AD <- c("#5F7F3C")
AU <- c("#C6D573")
MD <- c("#9D2919")
MU <- c("#BE261B")
PU <- c("#CB6078")
GD<- c("#3B2778")
GU <- c("#9772C3")
TD <- c("#1D4A8A")
TU <- c("#5F84B2")
ED <- c("#F4B740")
EU <- c("#D7AD7F")


## Filter for only chromosomes with "chr" prefix
GHP_data <- GHP_data %>%
  filter(grepl("^chr", chrom)) %>%
  filter(!is.na(mean_het))

## Order chromosomes numerically
GHP_data$chrom <- factor(GHP_data$chrom, 
                         levels = paste0("chr", 1:23))

## create window positions within each chromosome
GHP_data$window_position <- as.numeric(factor(GHP_data$window))


## GHP plot
GHP_plot <- ggplot(GHP_data, aes(x = chrom, y = mean_het, group = window_position)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9, fill = "#3B2778") +  # Bars for heterozygosity
  theme_bw() +
  labs(x = "Chromosome", y = "Het/50 Kb") +
  # theme(axis.text.x = element_text(angle = 0, hjust = 1)) +  # Rotate x labels
  scale_x_discrete(labels = function(x) gsub("chr", "", x)) +  # Clean x-axis labels
  scale_y_continuous(breaks=seq(0,1,0.2),limits=c(0,0.5),labels=seq(0,1,0.2))+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(size=6,family="Times"),
        axis.text.y=element_text(size=6,family="Times"),
        axis.title=element_text(size=8,family="Times"))


ggsave("../figs/GHP_sawtooth.png", GHP_plot, width = 8, height = 2, units="cm",dpi = 400)
