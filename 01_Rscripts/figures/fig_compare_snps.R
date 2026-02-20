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
McCrone_IAV_recall <- read_csv("03_data/table_S1.csv")
colnames(McCrone_IAV_recall) <- c("ID", "iSNV", "DPS1", "DPS2", "freq1", "freq2", "Subset")
McCrone_IAV <- read_csv("03_data/elife-35962-fig2-data5-v3.csv")
McCrone_IAV_nocut <- read_csv("03_data/no_freq_cut.qual.snv.csv") [, c("ENROLLID", "mutation", "freq.var", "onset", "collect")]
IDmap <- read_tsv("03_data/all_aln_map.tsv")

# 5. McCrone ====
## 5.1. Clean McCrone no cut data ====
McCrone_IAV_nocut <- McCrone_IAV_nocut %>%
  mutate(DPS = as.numeric(difftime(collect, onset, units = "days"))) %>%
  rename(ID = ENROLLID,
         iSNV = mutation,
         freq_mccrone = freq.var)

McCrone_obs1_recall <- McCrone_IAV_recall %>%
  select(iSNV, ID, DPS1, freq1) %>%
  rename(DPS = DPS1,
         freq_recall = freq1) %>%
  mutate(ID = as.character(ID))

McCrone_obs2_recall <- McCrone_IAV_recall %>%
  select(iSNV, ID, DPS2, freq2) %>%
  rename(DPS = DPS2,
         freq_recall = freq2) %>%
  mutate(ID = as.character(ID))

McCrone_recall <- rbind(McCrone_obs1_recall, McCrone_obs2_recall)

## 5.2. Match and combine data frames ====
# Function to extract numeric part of iSNV
extract_numeric <- function(isnv_string) {
  str_extract(isnv_string, "(?<=_)[A-Z](\\d+)[A-Z]") %>%  # Extract letter + number + letter
    str_extract("\\d+")  # Extract only the number
}

# Process both datasets
McCrone_IAV_nocut <- McCrone_IAV_nocut %>%
  mutate(iSNV_match = str_replace(iSNV, "(_.).*(.)$", "\\1\\2"),
         iSNV_num = extract_numeric(iSNV),
         iSNV_num = as.numeric(iSNV_num))  

McCrone_recall <- McCrone_recall %>%
  mutate(iSNV_match = str_replace(iSNV, "(_.).*(.)$", "\\1\\2"),
         iSNV_num = extract_numeric(iSNV),
         iSNV_num = as.numeric(iSNV_num))  

# Create all possible pairs within each ID and iSNV_match group
McCrone_tot_match <- McCrone_IAV_nocut %>%
  right_join(McCrone_recall, by = c("ID", "iSNV_match"), suffix = c("_mccrone", "_recall"), relationship = "many-to-many") %>%
  mutate(num_diff = abs(iSNV_num_mccrone - iSNV_num_recall)) %>%
  group_by(ID, iSNV_match, iSNV_mccrone) %>%
  slice_min(num_diff, n = 1, with_ties = FALSE) %>%  # Pick closest match
  ungroup() %>%
  select(iSNV_mccrone, iSNV_recall, ID, DPS_mccrone,
         freq_recall, freq_mccrone) %>%
  distinct() %>%
  rename(DPS = DPS_mccrone) %>%
  mutate(freq_mccrone = replace_na(freq_mccrone, -0.1),
         color_type = if_else(freq_mccrone == -0.1, "grey", "blue"))

# 6. Graphs ====
figS1 <- ggplot(McCrone_tot_match, aes(x = freq_mccrone, y = freq_recall)) +
  geom_point(aes(color = color_type), alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") + # x = y line
  scale_x_continuous(limits = c(-0.1, 1),
                     breaks = c(-0.1, 0, 0.25, 0.5, 0.75, 1.0),
                     labels = c('NA', 0, 0.25, 0.5, 0.75, 1.0),
                     name = 'Human IAV iSNV frequencies called by McCrone et al. (2018)') + 
  scale_y_continuous(limits = c(0, 1),
                     breaks = c(-0.1, 0, 0.25, 0.5, 0.75, 1.0),
                     labels = c('NA', 0, 0.25, 0.5, 0.75, 1.0),
                     name = 'Human IAV iSNV frequencies called by our pipeline') +
  scale_color_manual(values = c("#005D9E", "grey")) +
  theme(legend.position = "none")

ggsave(figS1, file = "figS1.pdf", path = "02_plots", width = 8, height = 6)
