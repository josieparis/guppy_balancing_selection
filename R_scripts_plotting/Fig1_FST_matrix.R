# new general run
rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

## load libs
lib<-c("dplyr","ggplot2","tidyverse","viridis")
lapply(lib,library,character.only=T)

setwd("<dir>/Analyses/NFDS_analysis/pop_structure")

fst <- read_tsv("pairwise_FST_matrix_R.tsv")
fst$pop1 <- factor(fst$pop1,levels=c("EU","ED","GD","GU","TD","TU","AD","AU",
                                     "MD","MU","PU"))
fst$pop2 <- factor(fst$pop2,levels=c("EU","ED","GD","GU","TD","TU","AD","AU",
                                     "MD","MU","PU"))

# Filter for upper triangle (pop1 level number < pop2 level number)
fst_upper <- fst[as.numeric(fst$pop1) < as.numeric(fst$pop2), ]

heatmap_fst <- ggplot(fst_upper, aes(pop1, pop2)) +
  theme_bw() +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0, vjust = 0),
        axis.text.x.top = element_text(margin = margin(b = 15)),
        axis.text.y = element_text(margin = margin(r = 15)),
        panel.border=element_blank(),
        panel.grid=element_blank()) +
  xlab("") + ylab("") +
  geom_tile(aes(fill = FST), color='white') +
  scale_fill_viridis(option = "plasma", direction = 1,
                     name = expression(italic(F)[ST])) +
  scale_x_discrete(position = "top") +
  coord_fixed(ratio = 1)

heatmap_fst

ggsave("../figs/FST_matrix.pdf", heatmap_fst, device="pdf", units="cm", width=15, height=12, limitsize=FALSE)

