# new general run
rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

setwd("<dir>/Analyses/NFDS_analysis/popgen_stats/")

lib<-c("dplyr","ggplot2","tidyverse","viridis")
lapply(lib,library,character.only=T)

fst <- read.table("data/ALL.50kb.fst.popgenome.out")
colnames(fst) <- c("chrom", "window",  "window_start",  "window_end", "GHP_GLP", "GHP_TUHP", "GHP_TULP", "GHP_APHP", "GHP_APLP", "GHP_ECHP", "GHP_ECLP", "GHP_MHP", "GHP_MLP", "GHP_P", "GLP_TUHP", "GLP_TULP", "GLP_APHP", "GLP_APLP", "GLP_ECHP", "GLP_ECLP", "GLP_MHP", "GLP_MLP", "GLP_P", "TUHP_TULP", "TUHP_APHP", "TUHP_APLP", "TUHP_ECHP", "TUHP_ECLP", "TUHP_MHP", "TUHP_MLP", "TUHP_P", "TULP_APHP", "TULP_APLP", "TULP_ECHP", "TULP_ECLP", "TULP_MHP", "TULP_MLP", "TULP_P", "APHP_APLP", "APHP_ECHP", "APHP_ECLP", "APHP_MHP", "APHP_MLP", "APHP_P", "APLP_ECHP", "APLP_ECLP", "APLP_MHP", "APLP_MLP", "APLP_P", "ECHP_ECLP", "ECHP_MHP", "ECHP_MLP", "ECHP_P", "ECLP_MHP", "ECLP_MLP", "ECLP_P", "MHP_MLP", "MHP_P", "MLP_P")
summary(fst)

## make long
fst_long <- fst %>%
  gather(key = "Comparison", value = "FST", -chrom, -window, -window_start, -window_end)

## calculate mean FST for each population pair
fst_summary <- fst_long %>%
  group_by(Comparison) %>%
  summarise(mean_FST = mean(FST, na.rm = TRUE))

write.table(fst_summary,file="data/MEAN_FST_summary.tsv",quote=FALSE,sep=",",row.names=F)
