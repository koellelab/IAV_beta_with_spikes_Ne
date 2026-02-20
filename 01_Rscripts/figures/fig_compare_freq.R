# 1. Source functions ====
source("01_Rscripts/WithinhostIAV_Ne_fnuctions_log.R")

# 2. Load packages ====
library(tidyverse)
library(ggplot2)
library(patchwork)
library(showtext)
library(colorspace)
library(parallel)
library(stringr)
theme_set(theme_minimal() +
            theme(axis.ticks = element_line(colour = "grey50"),
                  axis.line = element_line(colour = "grey50"),
                  axis.text = element_text(size = 12),
                  axis.title = element_text(size = 15),
                  legend.text = element_text(size = 12)))

# 3. Input parameters ====
Ne_values <- seq(10, 100, 1)
mu <- 0
alpha <- 0.02
hourPerGen <- 6
m <- 0.004
freq_step <- 0.01
set.seed(12)

# 4. Input datasets ====
compare_snps <- read.csv("03_data/jt_mm_compare_snps.csv")
IDmap <- read_tsv("03_data/all_aln_map.tsv")

# 5. Generate data frame ====

# 6. Graphs ====
fig_compare_snps <- ggplot(compare_snps, aes(x = minor_allele_x_freq_y, y = minor_allele_x_freq_x)) +
  geom_point(alpha = 0.5, aes(color = duplicate), ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") + # x = y line
  labs(x = "Human IAV iSNV frequencies called by McCrone et al. (2018)",
       y = "Human IAV iSNV frequencies called by our pipeline") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#005D9E", "grey")) +
  xlim(0, 1) + 
  ylim(0, 1)

ggsave(fig_compare_snps, file = "fig_compare_snps.pdf", path = "02_plots", width = 8, height = 8)
