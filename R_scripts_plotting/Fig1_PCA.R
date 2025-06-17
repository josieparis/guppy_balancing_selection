# new general run
rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

# Plot PCA
lib<-c("ggplot2","gridExtra","grid","extrafont")
lapply(lib,library,character.only=T)

setwd("<dir>/Analyses/NFDS_analysis/pop_structure/")

## eigenvec_table
eigenvec_table <- read.table('./data/holi_11_plink_out_NoLD_PCA.eigenvec', header = FALSE)

## add populations
eigenvec_table$Populations <- as.character(c(rep("AHP", 19), rep("ALP", 18), rep("EHP", 19), rep("ELP", 19),
                                             rep("GHP", 19), rep("GLP", 18),rep("MHP", 20),rep ("PLP", 10),
                                             rep("THP", 19),rep("TLP", 17),rep("MLP",17)))

## add shapes
eigenvec_table$shape<-(c(rep("21",19), rep("24",18), rep("21",19),rep("24",19),
                         rep("21",19),rep("24",18), rep("21",20), rep("24",10),
                         rep("21",19), rep("24",17), rep("24",17)))


## 21 is circle (downstream poplns)
## 24 is triangle (upstream poplns)

## check for order of pops and samples:
eigenvec_table <- eigenvec_table %>%
  dplyr::select(Populations, everything()) 

## remove repeated sample name column
eigenvec_table <- eigenvec_table[-2]

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

## order
pop_order <- c("AHP","ALP","EHP","ELP","GHP","GLP","MHP","MLP","THP","TLP","PLP")

## palette
palette <- c("#5F7F3C","#C6D573","#F4B740","#D7AD7F","#3B2778","#9772C3","#9D2919","#BE261B","#1D4A8A","#5F84B2","#CB6078")


names(palette)<-pop_order

## paste PCAs
colnames(eigenvec_table)[3:12] <- paste0("PC", 1:10)

## read in eigenval 
eigenval <- read.table('data/holi_11_plink_out_NoLD_PCA.eigenval', header = F)
percentage <- round(eigenval$V1/sum(eigenval$V1)*100,2)
percentage <- paste0(colnames(eigenvec_table)[3:12]," (",paste(as.character(percentage),"%)"))

## factor
eigenvec_table$Populations <- factor(eigenvec_table$Populations, levels=pop_order)

## plot PC1 vs PC2
PC1_PC2 <-ggplot(eigenvec_table,aes(x=PC1,y=PC2,fill=Populations,shape=shape),colour="black")+
  geom_point(size=2.8,stroke=0.3)+ 
  # stat_ellipse(type = "norm",level=0.99) +
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.title=element_text(size=16, family = "Times"),
        axis.text=element_text(size=12, family = "Times"),
        axis.line = element_line(colour = "black"),
        #      legend.position=("right"),
        #       legend.title=element_text(size=14, family = "Avenir"),
        legend.position="none",
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(size=11, family = "Times"))+
  # scale_colour_manual(values=palette)+
  scale_fill_manual(values = alpha(palette,0.9))+
  scale_shape_manual(values = c(21, 24)) + 
  xlab(percentage[1])+
  ylab(percentage[2])

PC1_PC2

ggsave("./figs/PCA_PC1_PC2.png", PC1_PC2, width = 8, height = 8, units="cm",dpi = 400)

## plot PC1 vs PC3
PC1_PC3 <-ggplot(eigenvec_table,aes(x=PC1,y=PC3,fill=Populations,shape=shape),colour="black")+
  geom_point(size=2.8,stroke=0.3)+ 
  # stat_ellipse(type = "norm",level=0.99) +
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.title=element_text(size=16, family = "Times"),
        axis.text=element_text(size=12, family = "Times"),
        axis.line = element_line(colour = "black"),
        #      legend.position=("right"),
        #       legend.title=element_text(size=14, family = "Avenir"),
        legend.position="none",
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(size=11, family = "Times"))+
  # scale_colour_manual(values=palette)+
  scale_fill_manual(values = alpha(palette,0.9))+
  scale_shape_manual(values = c(21, 24)) + 
  xlab(percentage[1])+
  ylab(percentage[3])

PC1_PC3

ggsave("./figs/PCA_PC1_PC3.png", PC1_PC3, width = 8, height = 8, units="cm",dpi = 400)

