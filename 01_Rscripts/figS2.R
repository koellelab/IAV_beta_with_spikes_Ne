# 1. Source functions ====
source("01_Rscripts/WithinhostIAV_Ne_fnuctions.R")

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
dv <- read_tsv("03_data/data/swineIAV_iSNV_table.run2_single.txt", show_col_types = FALSE,  col_types = cols(pigID = 'c'))
mm <- read_csv("03_data/data/processed_vcf.csv", show_col_types=FALSE,  col_types = cols(ID = 'c'))
dates <- read_tsv("03_data/data/all_dates.tsv", show_col_types=FALSE)
pos_map <- read_tsv("03_data/data/genomes/all_aln_map.tsv", show_col_types=FALSE) %>% 
  mutate(seg = replace_na(seg, "NA"),
         ID = as.character(ID)) %>%
  rename(c("segment" = "seg"))
dp <- read_csv("03_data/data/use_dp.csv", show_col_types=FALSE) %>%
  mutate(segment = if_else(segment == "MP", "M", segment),
         ID = as.character(ID),
         segment = replace_na(segment, "NA"))

# 5. FORMAT MM AND DV DATA ====
# in MM data NA AFs are assigned a -1 value
mm = mm %>%
  filter(data == "swine") %>%
  mutate(
    ID = as.character(ID),
    segment = case_when(
      segment == "MP" ~ "M",
      TRUE ~ segment),
    segment = replace_na(segment, "NA"),
    AF = case_when(
      AF == -1 ~ NA,
      TRUE ~ AF))

# DV data is in wide format so need to convert to long
dv = dv %>% 
  pivot_longer(
    -c('pigID', 'subtype', 'segment', 'ref_nt',
       'iSNV_loc', 'iSNV_nt', 'NS-S')) %>%
  filter(!is.na(value)) %>%
  mutate(
    iSNV_loc = iSNV_loc + 1,
    ID = as.character(pigID),
    type = str_split(name, "_", simplify=TRUE)[,2],
    d = gsub('d', '', str_split(name, "_", simplify=TRUE)[,1])) %>%
  select(-name, -pigID) %>%
  # wider format to calculate allele frequency 
  pivot_wider(names_from=type, values_from=value) %>%
  mutate(
    AF = sub/tot,
    d = as.integer(d),
    reference = segment,
    segment = str_split(segment, "_", simplify=TRUE)[,2],
    segment = if_else(segment == "NEP", "NS", segment)) %>%
  # add in dates since they're not in this file by default
  # inner merge to drop observations in the DV file not in our analysis
  inner_join(dates,
             by=c('ID', 'd'))

# 6. ADD UNIVERSAL POS AND MERGE ====
mm = mm %>%
  left_join(pos_map,
            by=c(
              'data'='data',
              'ID' = 'ID',
              'segment' = 'segment',
              'iSNV_loc' = 'pos',
              'acc' = 'acc',
              'date' = 'date')) %>%
  mutate(mm = TRUE)

dv = dv %>%
  left_join(pos_map %>% select(-date, -acc),
            by=c(
              'data'='data',
              'reference' = 'ID',
              'segment' = 'segment',
              'iSNV_loc' = 'pos')) %>%
  mutate(dv = TRUE)

mm_dv = mm %>% 
  select(ID, date, date, d, segment, upos, ref_nt, iSNV_nt, AF, SB, mm) %>%
  full_join(dv %>%
              select(ID, date, d, segment, upos, ref_nt, iSNV_nt, AF, dv),
            by=c('ID', 'date', 'd', 'segment', 'upos')) %>%
  # re-add positions in my coordinate system for variants only in dv
  left_join(pos_map %>% select(segment, ID, date, upos, pos),
            by=c('segment', 'ID', 'date', 'upos'))


# add in all depth and set my AF to 0 if suitable depth
# assuming DV file is already filtered for depth
mm_dv = mm_dv %>%
  left_join(dp %>% filter(!is.na(acc)) %>% select(ID, date, acc) %>% unique(),
            by=c('ID', 'date')) %>%
  left_join(dp,
            by=c('segment', 'ID', 'date', 'acc', 'pos')) %>%
  mutate(AF.x = if_else(is.na(AF.x) & depth > 500, 0, AF.x))

# 7. SWAP REFERENCE ALLELES IF THEY DO NOT MATCH ====
swaps = ((!is.na(mm_dv$ref_nt.x)) & 
           (!is.na(mm_dv$ref_nt.y)) & 
           (mm_dv$ref_nt.x != mm_dv$ref_nt.y) & 
           (mm_dv$iSNV_nt.x == mm_dv$ref_nt.y))

mm_dv[swaps,]$iSNV_nt.y =  mm_dv[swaps,]$iSNV_nt.x
mm_dv[swaps,]$ref_nt.y = mm_dv[swaps,]$ref_nt.x
mm_dv[swaps,]$AF.y = 1-mm_dv[swaps,]$AF.y 

stopifnot(!any(((!is.na(mm_dv$ref_nt.x)) & 
                  (!is.na(mm_dv$ref_nt.y)) & 
                  (mm_dv$ref_nt.x != mm_dv$ref_nt.y) & 
                  (mm_dv$iSNV_nt.x == mm_dv$ref_nt.y))))


# 8.  PLOT FIGURES ====
# set NA values to -0.25 for plottinng
tmp = mm_dv %>%
  mutate(AF.x = replace_na(AF.x, -0.25),
         AF.y = replace_na(AF.y, -0.25))

# drop if NA or 0 in both
tmp = tmp %>% 
  filter((AF.x > 0) | (AF.y > 0)) %>%
  mutate(color_type = if_else(AF.y == -0.25, "grey", "blue"))


figS2 <- ggplot(tmp %>% filter(depth >= 500), aes(x=AF.y, y=AF.x)) + 
  geom_point(aes(color = color_type), alpha = 0.2) + 
  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") + # x = y line
  scale_x_continuous(limits = c(-0.25, 1),
                     breaks=c(-0.25, 0, 0.25, 0.5, 0.75, 1.0),
                     labels = c('NA', 0, 0.25, 0.5, 0.75, 1.0),
                     name='Swine IAV iSNV frequencies called by VanInsberghe et al. (2024)') + 
  scale_y_continuous(limits = c(0, 1),
                     breaks=c(-0.25, 0, 0.25, 0.5, 0.75, 1.0),
                     labels = c('NA', 0, 0.25, 0.5, 0.75, 1.0),
                     name='Swine IAV iSNV frequencies called by our pipeline') +
  scale_color_manual(values = c("#005D9E", "grey")) +
  theme(legend.position = "none")

ggsave(figS2, file = "figS2.pdf", path = "02_plots", width = 8, height = 6)
