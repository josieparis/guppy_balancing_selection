# new general run
rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

## load libs
lib<-c("dplyr","ggplot2","tidyverse","gggenes","viridis")
lapply(lib,library,character.only=T)

#detach("package:plyr", unload = TRUE)

setwd("<dir>/Analyses/NFDS_analysis/candidate_loci/data/")

## define pop order and rename
pop_order <- c("GHP", "GLP", "ECHP", "ECLP", "APHP", "APLP", "TUHP", "TULP", "MHP", "MLP", "P")
pop_rename <- c("GD", "GU", "ED", "EU", "AD", "AU", "TD", "TU", "MD", "MU", "PU")
pop_map <- setNames(pop_rename, pop_order)

## read in raw ballermix data to plot
chr22_v1 <- read.csv("chr22_13bp_for_plotting_v2.txt",h=T,sep="\t")
chr22_v1$pop <- factor(chr22_v1$pop, levels = pop_order)

## log-transform the s_hat
chr22_v1 <- chr22_v1 %>%
  mutate(log10_s_hat = log10(s_hat))

### window with small extension
chr22_v1 <- chr22_v1 %>% filter(physPos >= 13045000 & physPos <=13122138)

## get max CLR 
max_clr_rows <- chr22_v1 %>%
  group_by(pop) %>%
  slice_max(order_by = CLR, n = 1, with_ties = FALSE) %>%
  ungroup()

## read in the genes
genes <-read.csv("genes.txt",h=T,sep="\t")

all_pops <- unique(chr22_v1$pop)

chr22_v1_genes <- genes %>%
  filter(window == "chr22_13100000") %>%
  crossing(pop = all_pops)

# Get positions for plotting the genes below the CLR points (subtract a percentage of the range)
gene_y_positions <- chr22_v1 %>%
  group_by(pop) %>%
  dplyr::summarise(y = min(CLR, na.rm = TRUE) - 0.15 * diff(range(CLR, na.rm = TRUE)))

chr22_v1_genes <- chr22_v1_genes %>%
  left_join(gene_y_positions, by = "pop")


## rename genes with their annotations
gene_names <- c(
  "g40432" = "iars2",
  "g40433" = "nlgn3a",
  "g40434" = "qrsl1",
  "g40435" = "rtn4ip1",
  "g40436" = "slc45a2")

chr22_v1_genes <- chr22_v1_genes %>%
  mutate(gene_ID = gene_names[gene_ID])

chr22_v1_genes$gene_ID <- as.character(chr22_v1_genes$gene_ID)

## colour palette for genes
gene_palette <- c(
  "iars2" = "orchid",
  "nlgn3a" = "brown",
  "qrsl1" = "seagreen",
  "rtn4ip1" = "skyblue",
  "slc45a2" = "tomato")

custom_breaks <- seq(13040000, 13500000, by = 10000)

## make sure everything is factored by pop_order
chr22_v1$pop <- factor(chr22_v1$pop, levels = pop_order)
chr22_v1_genes$pop <- factor(chr22_v1_genes$pop, levels = pop_order)

### add CLR thresholds:
clr_thresholds <- tibble(
  pop = c("GHP", "GLP", "ECHP", "ECLP", "APHP", "APLP", "TUHP", "TULP", "MHP", "MLP", "P"),
  clr_99 = c(142.73, 308.04, 170.27, 221.65, 85.71, 138.16, 109.03, 123.15, 160.99, 284.32, 232.83))

clr_thresholds$pop <- factor(clr_thresholds$pop, levels = pop_order)

## Update y-limits for plotting (5 above the threshold)
y_ranges_updated <- clr_thresholds %>%
  mutate(ymax = clr_99 + 5)

chr22_v1 <- chr22_v1 %>%
  left_join(y_ranges_updated %>% dplyr::select(pop, ymax), by = "pop")

## Recompute ymin (if needed)
chr22_v1 <- chr22_v1 %>%
  group_by(pop) %>%
  mutate(ymin = min(CLR, na.rm = TRUE) - 0.2 * diff(range(CLR, na.rm = TRUE))) %>%
  ungroup()

## plotting:
LG22_plot <- ggplot(chr22_v1, aes(x = physPos, y = CLR)) +
  geom_point(aes(colour = log10_s_hat), alpha = 0.8, size = 0.4) +
  geom_hline(data = clr_thresholds,
    aes(yintercept = clr_99),
    linetype = 3,
    linewidth = 0.3)+
  geom_blank(aes(y = ymin))+
  geom_blank(aes(y = ymax))+
  geom_gene_arrow(data = chr22_v1_genes,
    inherit.aes = FALSE,
    aes(xmin = gene_start,xmax = gene_end,y = y,forward = strand == "+",fill = gene_ID,group = pop),
    arrowhead_width = unit(1, "mm"),
    arrowhead_height = unit(1.2,"mm"),
    arrow_body_height = unit(0.8, "mm"))+
  geom_gene_label(data = chr22_v1_genes,
                  inherit.aes = FALSE,
                  aes(xmin = gene_start, xmax = gene_end, y = y, label = gene_ID,group = pop),
                  angle = 45, size = 3)+
  facet_grid(rows = vars(pop), scales = "free_y", labeller = as_labeller(pop_map)) +
  #scale_colour_manual(values = pop_cols) +
  scale_colour_viridis(name = bquote(log[10] * hat(s)), option = "C", direction = -1, na.value = "grey80")+
  scale_fill_manual(values = gene_palette, name = NULL,
    labels = c("iars2" = expression(italic("iars2")),
               "nlgn3a" = expression(italic("nlgn3a")),
               "qrsl1" = expression(italic("qrsl1")),
               "rtn4ip1" = expression(italic("rtn4ip1")),
               "slc45a2" = expression(italic("slc45a2")))) +
  scale_x_continuous(breaks = custom_breaks, labels = function(x) paste0(x / 1e6))+
  coord_cartesian(clip = "on")+
  theme_minimal() +
  labs(x = "LG22 (Mb)", y = expression(italic(B)[2]), colour = "Population") +
  theme(axis.text.x = element_text(size = 6,family="Times New Roman"),
        axis.text.y = element_text(size = 4,family="Times New Roman"),
        axis.title = element_text(size = 12,family="Times New Roman"),
        strip.text = element_text(size = 12,family="Times New Roman"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        axis.ticks = element_line(colour = "black"),
        strip.placement = "outside",
        axis.line.x = element_line(),
        strip.switch.pad.wrap = unit(0.2, "cm"),
        theme(legend.position = "right"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),     
        legend.key.size = unit(0.3, "cm"))

ggsave("../../figs/LG22_candidate_plot.png", LG22_plot, width = 8, height = 15, units = "cm")
