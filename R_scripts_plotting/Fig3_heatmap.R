# new general run
rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

## load libs
lib<-c("dplyr","ggplot2","tidyverse","forcats")
lapply(lib,library,character.only=T)

setwd("~/Dropbox/Sussex_Guppies/Analyses/NFDS_analysis/heatmap_B2_outliers/")

## read in B2 results for all chromosomes 
B2 <- read.csv("~/Dropbox/Sussex_Guppies/Analyses/NFDS_analysis/Ballermix_bonnie_results/merged_ALL_data_50kb_no_up.txt",sep="\t",h=T)

outliers <- read.csv("outlier_IDs_23_windows.txt",h=T,sep="\t")

## get only the windows we want
B2_filtered <- B2 %>%
  semi_join(outliers, by = "ID")

pops <- c("GHP", "GLP", "ECHP", "ECLP", "APHP", "APLP", "TUHP", "TULP", "MHP", "MLP", "P")
pop_rename <- c("GD", "GU", "ED", "EU", "AD", "AU", "TD", "TU", "MD", "MU", "PU")
pop_map <- setNames(pop_rename, pops)

# make it long 
logS_long <- B2_filtered %>%
  # Ensure all logS columns are numeric
  mutate(across(matches("_logS$"), ~ as.numeric(.))) %>%
  dplyr::select(ID, matches("_logS$")) %>%
  pivot_longer(
    cols = -ID,
    names_to = "pop",
    names_pattern = "(.*)_logS$",
    values_to = "logS"
  )

## bin
logS_binned <- logS_long %>%
  mutate(category = case_when(
    logS < 0 ~ "<0",
    logS >= 0 & logS < 3 ~ "0–3",
    logS >= 3 & logS < 6 ~ "3–6",
    logS >= 6 & logS < 9 ~ "6–9",
    logS == 9 ~ ">9"))

logS_binned$category <- factor(
  logS_binned$category,
  levels = c("<0", "0–3", "3–6", "6–9", ">9"))

logS_binned <- logS_binned %>%
  mutate(pop = factor(pop, levels = pops, labels = pop_rename))

logS_binned <- logS_binned %>%
  mutate(
    chr = as.numeric(gsub("chr(\\d+)_.*", "\\1", ID)),
    end_pos = as.numeric(gsub(".*_(\\d+)", "\\1", ID)),
    start_pos = end_pos - 500,  # assuming 500 bp windows
    start_mb = formatC(start_pos / 1e6, format = "f", digits = 2),
    ID_label = paste0("LG", chr, ", ", start_mb, " Mb")
  ) %>%
  arrange(chr, end_pos) %>%
  mutate(ID_label = factor(ID_label, levels = unique(ID_label)))


## calculate CLR and add white circles to plot:
## use the B2 windowed results (called B2 here)

# convert CLR scores to long df
clr_long <- B2 %>%
  dplyr::select(ID, matches("_CLR$")) %>%
  pivot_longer(
    cols = -ID,
    names_to = "pop",
    names_pattern = "(.*)_CLR$",
    values_to = "CLR"
  )

## recode population names
clr_long <- clr_long %>%
  mutate(pop = factor(pop, levels = pops, labels = pop_rename))

### 99% CLR and filter
clr_outliers <- clr_long %>%
  group_by(pop) %>%
  mutate(threshold_99 = quantile(CLR, 0.99, na.rm = TRUE)) %>%
  filter(CLR >= threshold_99) %>%
  ungroup()

## add ID_label
clr_outliers <- clr_outliers %>%
  mutate(
    chr = as.numeric(gsub("chr(\\d+)_.*", "\\1", ID)),
    end_pos = as.numeric(gsub(".*_(\\d+)", "\\1", ID)),
    start_pos = end_pos - 500,
    start_mb = formatC(start_pos / 1e6, format = "f", digits = 2),
    ID_label = paste0("LG", chr, ", ", start_mb, " Mb"))

### only keep the ones also identified in picmin
clr_outliers_filtered <- clr_outliers %>%
  semi_join(logS_binned, by = c("pop", "ID_label"))

### plot
plot <- ggplot(logS_binned, aes(x = pop, y = fct_rev(ID_label), fill = category)) +
  geom_tile(color = "black",linewidth = 0.2) +
  scale_fill_viridis_d(
    name = bquote(log[10] * hat(s)),
    option = "C", direction = -1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,family="Times New Roman",size=12),
    axis.text = element_text(size = 6,family="Times New Roman"),
    axis.title = element_blank())+
  geom_point(
    data = clr_outliers_filtered,
    aes(x = pop, y = ID_label),
    shape = 21, size = 1.5, fill = "white", inherit.aes = FALSE)

ggsave("../figs/heatmap_outliers.png", plot, width = 12, height = 10, units = "cm")
