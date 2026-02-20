# 1. Source ====
source("01_Rscripts/WithinhostIAV_Ne_fnuctions_log.R")
source("01_Rscripts/fig4.R")

# 2. Load packages ====
library(tidyverse)
library(purrr)
library(ggplot2)
library(patchwork)
library(showtext)
library(colorspace)
library(parallel)
library(readxl)
library(reshape)
library(ggrepel)
theme_set(theme_minimal() +
            theme(axis.ticks = element_line(colour = "grey50"),
                  axis.line = element_line(colour = "grey50"),
                  axis.text = element_text(size = 12),
                  axis.title = element_text(size = 15),
                  legend.text = element_text(size = 12),
                  plot.title = element_text(size = 20)))

# 3. Input parameters ====
Ne_values <- seq(3, 2000, 5)
mu <- 0
alpha <- 0.02
hourPerGen <- 6
m <- 0.004

# 4. Input datasets ====
swineData <- read_csv("03_data/table_S2.csv") %>%
  filter(Subset == 1) %>%
  set_names(c("pigID", "segment", "DPS1", "DPS2", "freq1", "freq2", "Subset")) %>%
  mutate(gen = (DPS2 - DPS1)*(24/hourPerGen))

# outliers
outliers <- swineData %>%
  filter(pigID %in% c("20578", "37279", "37317", "37356"))

# 5. Generate data frame ====
## 5.1. Calculate MLE for each pair ====
MLE_swine_inliers_table <- swineData %>%
  filter(!pigID %in% c("20578", "37279", "37317", "37356")) %>%
  mutate(pig_segment = paste0(pigID, "_", segment)) %>%
  split(seq_len(nrow(.))) %>% 
  map_dfr(function(df1) {
    LLdf <- LL_Ne.df(df1, Ne_values, mu, alpha, hourPerGen, m)
    mle_ne <- LLdf$Ne[which.max(LLdf$prob)]
    tibble(
      pig_segment = paste0(df1$pigID, "_", df1$segment, "_", df1$DPS1),
      MLE_Ne = mle_ne
    )
  })

MLE_swine_outliers_table <- outliers %>%
  mutate(pig_segment = paste0(pigID, "_", segment)) %>%
  split(seq_len(nrow(.))) %>% 
  map_dfr(function(df1) {
    LLdf <- LL_Ne.df(df1, seq(2200, 4500, 5), mu, alpha, hourPerGen, m)
    mle_ne <- LLdf$Ne[which.max(LLdf$prob)]
    tibble(
      pig_segment = paste0(df1$pigID, "_", df1$segment, "_", df1$DPS1),
      MLE_Ne = mle_ne
    )
  })

MLE_swine_table <- rbind(MLE_swine_inliers_table, MLE_swine_outliers_table)

# Identify outliers by the standard 1.5*IQR rule (based on the full distribution)
q1  <- quantile(MLE_swine_table$MLE_Ne, 0.25, na.rm = TRUE)
q3  <- quantile(MLE_swine_table$MLE_Ne, 0.75, na.rm = TRUE)
iqr <- q3 - q1
lo  <- q1 - 1.5 * iqr
hi  <- q3 + 1.5 * iqr

MLE_swine_table_labeled <- MLE_swine_table %>%
  arrange(MLE_Ne) %>%
  mutate(
    is_outlier = (MLE_Ne < lo) | (MLE_Ne > hi),
    x_index = row_number(),
    color_class = case_when(
      MLE_Ne == 3        ~ "MLE_eq_3",
      is_outlier         ~ "outlier",
      TRUE               ~ "normal"
    )
  )

## 5.2. Adjust data to be plotted in iSNV frequency pair format ====
# Outliers
swine_outliers_id <- MLE_swine_table_labeled %>%
  filter(is_outlier == TRUE) %>%
  pull(pig_segment)

swineData_outliers <- swineData %>%
  mutate(pig_segment = paste0(pigID, "_", segment, "_", DPS1)) %>%
  filter(pig_segment %in% swine_outliers_id) %>%
  select(-c(pig_segment, gen))

# data frame only including first round of observations
swineData_outliers_obs1 <- swineData_outliers %>%
  select(pigID, DPS1, freq1)
swineData_outliers_obs1$type <- c("obs1")
names(swineData_outliers_obs1) <- c("pigID", "DPS", "freq", "type")
# data frame only including second round of observations
swineData_outliers_obs2 <- swineData_outliers %>%
  select(pigID, DPS2, freq2)
swineData_outliers_obs2$type <- c("obs2")
names(swineData_outliers_obs2) <- c("pigID", "DPS", "freq", "type")
# merge the data together
swineData_outliers_tot_obs <- rbind(swineData_outliers_obs1, swineData_outliers_obs2) %>%
  mutate(
    freq = if_else(freq < 0.02, 0.02, freq),
    type = factor(type, labels = c("observation 1", "observation 2")),
    link_type = case_when(
      type == "observation 1" ~ "obs1",
      type %in% "observation 2" ~ "obs2"
    )
  )
# adjust the linking index to make sure observation pairs are correctly linked
swineData_outliers_tot_obs <- swineData_outliers_tot_obs %>%
  arrange(pigID, DPS, desc(link_type)) %>%
  group_by(pigID) %>%
  mutate(
    pair_id = cumsum(link_type == "obs1")
  ) %>%
  ungroup()


# Low MLEs
swine_lows_id <- MLE_swine_table_labeled %>%
  filter(color_class == "MLE_eq_3") %>%
  pull(pig_segment)

swineData_lows <- swineData %>%
  mutate(pig_segment = paste0(pigID, "_", segment, "_", DPS1)) %>%
  filter(pig_segment %in% swine_lows_id) %>%
  select(-c(pig_segment, gen))

# data frame only including first round of observations
swineData_lows_obs1 <- swineData_lows %>%
  select(pigID, DPS1, freq1)
swineData_lows_obs1$type <- c("obs1")
names(swineData_lows_obs1) <- c("pigID", "DPS", "freq", "type")
# data frame only including second round of observations
swineData_lows_obs2 <- swineData_lows %>%
  select(pigID, DPS2, freq2)
swineData_lows_obs2$type <- c("obs2")
names(swineData_lows_obs2) <- c("pigID", "DPS", "freq", "type")
# merge the data together
swineData_lows_tot_obs <- rbind(swineData_lows_obs1, swineData_lows_obs2) %>%
  mutate(
    freq = if_else(freq < 0.02, 0.02, freq),
    type = factor(type, labels = c("observation 1", "observation 2")),
    link_type = case_when(
      type == "observation 1" ~ "obs1",
      type %in% "observation 2" ~ "obs2"
    )
  )
# adjust the linking index to make sure observation pairs are correctly linked
swineData_lows_tot_obs <- swineData_lows_tot_obs %>%
  arrange(pigID, DPS, desc(link_type)) %>%
  group_by(pigID) %>%
  mutate(
    pair_id = cumsum(link_type == "obs1")
  ) %>%
  ungroup()

# 6. Graphs ====
## 6.1. MLE scatter plot ====
fig_swine_outlier_points <- ggplot(MLE_swine_table_labeled, aes(x = x_index, y = MLE_Ne)) +
  geom_point(aes(color = color_class), size = 2, alpha = 0.85) +
  scale_color_manual(
    values = c(
      normal   = "black",
      outlier  = "#B22222",
      MLE_eq_3 = "#005D9E"
    ),
    guide = "none"
  ) +
  labs(x = "iSNV pair", y = expression("MLE of " * N[E])) +
  theme(
    axis.text.x  = element_blank(), 
    axis.ticks.x = element_blank(),panel.grid.minor = element_blank(),
    plot.margin = margin(5.5, 2, 5.5, 5.5) # reduce right margin to sit tight to boxplot
  )

## 6.2. Log 10 MLE scatter plot ====
fig_swine_outlier_points_log <- fig_swine_outlier_points +
  scale_y_log10() +
  labs(x = "iSNV pair", y = expression("MLE of " * N[E] * " (log10 scale)"))

## 6.3. iSNV frequency with outliers shown ====
fig_swine_outlier_freq <- ggplot() +
  # Background: all observations in gray
  geom_hline(yintercept = 0.02, color = "#B22222", linetype = "dashed") +
  geom_line(
    data = swineData_tot_obs,
    aes(x = DPS, y = freq, group = interaction(pigID, pair_id)),
    color = "gray",
    alpha = 0.3
  ) +
  geom_point(
    data = swineData_tot_obs,
    aes(x = DPS, y = freq, shape = type),
    color = "gray",
    alpha = 0.4
  ) +
  
  # Foreground: lows in blue
  geom_line(
    data = swineData_lows_tot_obs,
    aes(x = DPS, y = freq, group = interaction(pigID, pair_id)),
    color = "#005D9E",
    alpha = 0.4
  ) +
  geom_point(
    data = swineData_lows_tot_obs,
    aes(x = DPS, y = freq, shape = type),
    color = "#005D9E",
    alpha = 0.6
  ) +
  
  # Foreground: outliers in red
  geom_line(
    data = swineData_outliers_tot_obs,
    aes(x = DPS, y = freq, group = interaction(pigID, pair_id)),
    color = "#B22222",
    alpha = 0.4
  ) +
  geom_point(
    data = swineData_outliers_tot_obs,
    aes(x = DPS, y = freq, shape = type),
    color = "#B22222",
    alpha = 0.6
  ) +
  
  scale_shape_manual(values = c(19, 17, 2)) +
  labs(
    x = "Day of fair",
    y = "iSNV frequency",
    shape = NULL
  ) +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  #scale_y_continuous(limits = c(0, 1)) +
  theme(legend.position = "none")

# 7. Put graphs together ====
fig_swine_outlier_combined <- (fig_swine_outlier_points + fig_swine_outlier_points_log) / fig_swine_outlier_freq +
  plot_layout(heights = c(1, 1.5)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 18))

ggsave(fig_swine_outlier_combined, file = "fig_swine_outlier.pdf", 
       path = "02_plots", width = 10, height = 10)

write.csv(swineData_outliers, "03_data/table_swine_outliers.csv", row.names = FALSE, quote = FALSE)
write.csv(swineData_lows, "03_data/table_swine_lows.csv", row.names = FALSE, quote = FALSE)

