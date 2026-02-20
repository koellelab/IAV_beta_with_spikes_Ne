# 1. Source ====
source("01_Rscripts/WithinhostIAV_Ne_fnuctions_log.R")
source("01_Rscripts/fig2.R")

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
Ne_values <- seq(3, 700, 5)
mu <- 0
alpha <- 0.02
hourPerGen <- 6
m <- 0.004

# 4. Input datasets ====
McCrone_IAV <- read_csv("03_data/table_S1.csv") %>%
  set_names(c("ENROLLID", "mutation", "DPS1", "DPS2", "freq1", "freq2", "Subset")) %>%
  mutate(gen = (DPS2 - DPS1)*(24/hourPerGen)) %>%
  filter(Subset == 1)

# outliers
outliers_human <- McCrone_IAV %>%
  filter(ENROLLID %in% c("50418", "50727", "50246"))

# 5. Generate data frame ====
## 5.1. Calculate MLE for each pair ====
MLE_human_inliers_table <- McCrone_IAV %>%
  filter(!ENROLLID %in% c("50418", "50727", "50246")) %>%
  mutate(patient_iSNV = paste0(ENROLLID, "_", mutation)) %>%
  split(seq_len(nrow(.))) %>% 
  map_dfr(function(df1) {
    LLdf <- LL_Ne.df(df1, Ne_values, mu, alpha, hourPerGen, m)
    mle_ne <- LLdf$Ne[which.max(LLdf$prob)]
    tibble(
      patient_iSNV = paste0(df1$ENROLLID, "_", df1$mutation, "_", df1$DPS1),
      MLE_Ne = mle_ne
    )
  })
  
MLE_human_outliers_table <- outliers_human %>%
  mutate(patient_iSNV = paste0(ENROLLID, "_", mutation)) %>%
  split(seq_len(nrow(.))) %>% 
  map_dfr(function(df1) {
    LLdf <- LL_Ne.df(df1, seq(1300, 6500, 10), mu, alpha, hourPerGen, m)
    mle_ne <- LLdf$Ne[which.max(LLdf$prob)]
    tibble(
      patient_iSNV = paste0(df1$ENROLLID, "_", df1$mutation, "_", df1$DPS1),
      MLE_Ne = mle_ne
    )
  })

MLE_human_table <- rbind(MLE_human_inliers_table, MLE_human_outliers_table)


# Identify outliers by the standard 1.5*IQR rule (based on the full distribution)
q1  <- quantile(MLE_human_table$MLE_Ne, 0.25, na.rm = TRUE)
q3  <- quantile(MLE_human_table$MLE_Ne, 0.75, na.rm = TRUE)
iqr <- q3 - q1
lo  <- q1 - 1.5 * iqr
hi  <- q3 + 1.5 * iqr

MLE_human_table_labeled <- MLE_human_table %>%
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
human_outliers_id <- MLE_human_table_labeled %>%
  filter(is_outlier == TRUE) %>%
  pull(patient_iSNV)

human_outliers <- McCrone_IAV %>%
  mutate(patient_iSNV = paste0(ENROLLID, "_", mutation, "_", DPS1)) %>%
  filter(patient_iSNV %in% human_outliers_id) %>%
  select(-c(patient_iSNV, gen))

# data frame only including first round of observations
human_outliers_obs1 <- human_outliers %>%
  select(ENROLLID, DPS1, freq1)
human_outliers_obs1$type <- c("obs1")
names(human_outliers_obs1) <- c("ENROLLID", "DPS", "freq", "type")
# data frame only including second round of observations
human_outliers_obs2 <- human_outliers %>%
  select(ENROLLID, DPS2, freq2)
human_outliers_obs2$type <- c("obs2")
names(human_outliers_obs2) <- c("ENROLLID", "DPS", "freq", "type")
# merge the data together
human_outliers_tot_obs <- rbind(human_outliers_obs1, human_outliers_obs2) %>%
  mutate(
    freq = if_else(freq < 0.02, 0.02, freq),
    type = factor(type, labels = c("observation 1", "observation 2")),
    link_type = case_when(
      type == "observation 1" ~ "obs1",
      type %in% "observation 2" ~ "obs2"
    )
  )
# adjust the linking index to make sure observation pairs are correctly linked
human_outliers_tot_obs <- human_outliers_tot_obs %>%
  arrange(ENROLLID, DPS, desc(link_type))


# Low MLEs
human_lows_id <- MLE_human_table_labeled %>%
  filter(color_class == "MLE_eq_3") %>%
  pull(patient_iSNV)

human_lows <- McCrone_IAV %>%
  mutate(patient_iSNV = paste0(ENROLLID, "_", mutation, "_", DPS1)) %>%
  filter(patient_iSNV %in% human_lows_id) %>%
  select(-c(patient_iSNV, gen))

# data frame only including first round of observations
human_lows_obs1 <- human_lows %>%
  select(ENROLLID, DPS1, freq1)
human_lows_obs1$type <- c("obs1")
names(human_lows_obs1) <- c("ENROLLID", "DPS", "freq", "type")
# data frame only including second round of observations
human_lows_obs2 <- human_lows %>%
  select(ENROLLID, DPS2, freq2)
human_lows_obs2$type <- c("obs2")
names(human_lows_obs2) <- c("ENROLLID", "DPS", "freq", "type")
# merge the data together
human_lows_tot_obs <- rbind(human_lows_obs1, human_lows_obs2) %>%
  mutate(
    freq = if_else(freq < 0.02, 0.02, freq),
    type = factor(type, labels = c("observation 1", "observation 2")),
    link_type = case_when(
      type == "observation 1" ~ "obs1",
      type %in% "observation 2" ~ "obs2"
    )
  )
# adjust the linking index to make sure observation pairs are correctly linked
human_lows_tot_obs <- human_lows_tot_obs %>%
  arrange(ENROLLID, DPS, desc(link_type))

# 6. Graphs ====
## 6.1. MLE scatter plot ====
fig_human_outlier_points <- ggplot(MLE_human_table_labeled, aes(x = x_index, y = MLE_Ne)) +
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
fig_human_outlier_points_log <- fig_human_outlier_points +
  scale_y_log10() +
  labs(x = "iSNV pair", y = expression("MLE of " * N[E] * " (log10 scale)"))

## 6.3. iSNV frequency with outliers shown ====
fig_human_outlier_freq <- ggplot() +
  # Background: all observations in gray
  geom_hline(yintercept = 0.02, color = "#B22222", linetype = "dashed") +
  geom_line(
    data = McCrone_tot_obs,
    aes(x = DPS, y = freq, group = ENROLLID),
    color = "gray",
    alpha = 0.3
  ) +
  geom_point(
    data = McCrone_tot_obs,
    aes(x = DPS, y = freq, shape = type),
    color = "gray",
    alpha = 0.4
  ) +
  
  # Foreground: lows in blue
  geom_line(
    data = human_lows_tot_obs,
    aes(x = DPS, y = freq, group = ENROLLID),
    color = "#005D9E",
    alpha = 0.4
  ) +
  geom_point(
    data = human_lows_tot_obs,
    aes(x = DPS, y = freq, shape = type),
    color = "#005D9E",
    alpha = 0.6
  ) +
  
  # Foreground: outliers in red
  geom_line(
    data = human_outliers_tot_obs,
    aes(x = DPS, y = freq, group = ENROLLID),
    color = "#B22222",
    alpha = 0.4
  ) +
  geom_point(
    data = human_outliers_tot_obs,
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
  theme(legend.position = "none")

# 7. Put graphs together ====
fig_human_outlier_combined <- (fig_human_outlier_points + fig_human_outlier_points_log) / fig_human_outlier_freq +
  plot_layout(heights = c(1, 1.5)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 18))

ggsave(fig_human_outlier_combined, file = "fig_human_outlier.pdf", 
       path = "02_plots", width = 10, height = 10)

write.csv(human_outliers, "03_data/table_human_outliers.csv", row.names = FALSE, quote = FALSE)
write.csv(human_lows, "03_data/table_human_lows.csv", row.names = FALSE, quote = FALSE)

