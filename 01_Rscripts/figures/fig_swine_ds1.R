# 1. Source functions ====
source("01_Rscripts/WithinhostIAV_Ne_fnuctions_log.R")

# 2. Load packages ====
library(tidyverse)
library(ggplot2)
library(patchwork)
library(showtext)
library(colorspace)
library(parallel)
library(readxl)
library(reshape)
theme_set(theme_minimal() +
            theme(axis.ticks = element_line(colour = "grey50"),
                  axis.line = element_line(colour = "grey50"),
                  axis.text = element_text(size = 12),
                  axis.title = element_text(size = 15),
                  legend.text = element_text(size = 12),
                  plot.title = element_text(size = 20)))

# 3. Input parameters ====
Ne_values <- seq(3, 30, 0.1)
mu <- 0
alpha <- 0.02
hourPerGen <- 6
m <- 0.004

# 4. Input datasets ====
swineData <- read_csv("03_data/table_S2.csv") %>%
  filter(Subset == 1) %>%
  set_names(c("pigID", "segment", "DPS1", "DPS2", "freq1", "freq2", "Subset")) %>%
  mutate(gen = (DPS2 - DPS1)*(24/hourPerGen))

# 5. Generate data frame ====
## 5.1. Figure 4A ====
# data frame only including first round of observations
swineData_obs1 <- swineData %>%
  select(pigID, DPS1, freq1)
swineData_obs1$type <- c("obs1")
names(swineData_obs1) <- c("pigID", "DPS", "freq", "type")
# data frame only including second round of observations
swineData_obs2 <- swineData %>%
  select(pigID, DPS2, freq2)
swineData_obs2$type <- c("obs2")
names(swineData_obs2) <- c("pigID", "DPS", "freq", "type")
# merge the data together
swineData_tot_obs <- rbind(swineData_obs1, swineData_obs2) %>%
  mutate(
    type = if_else(freq < 0.02, "observation below 2%", type),
    freq = if_else(freq < 0.02, 0.02, freq),
    type = factor(type, labels = c("observation 1", "observation 2", "observation 2 below 2%")),
    link_type = case_when(
      type == "observation 1" ~ "obs1",
      type %in% c("observation 2", "observation 2 below 2%") ~ "obs2"
    )
  )
# adjust the linking index to make sure observation pairs are correctly linked
swineData_tot_obs <- swineData_tot_obs %>%
  arrange(pigID, DPS, desc(link_type)) %>%
  group_by(pigID) %>%
  mutate(
    pair_id = cumsum(link_type == "obs1")
  ) %>%
  ungroup()

## 5.2. Figure 4B ====
LLdf_swineData <- LL_Ne.df(swineData, Ne_values, mu, alpha, hourPerGen, m)
MLE_Ne_swineData <- LLdf_swineData$Ne[which.max(LLdf_swineData$prob)]
CI_low_Ne_swineData <- LLdf_swineData %>%
  filter(Ne < MLE_Ne_swineData) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_swineData <- LLdf_swineData %>%
  filter(Ne > MLE_Ne_swineData) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

## 5.3. Figure 4C ====
# create empty list
list_df_sim_bs_swine <- list()
for (i in 1:1000) {
  list_df_sim_bs_swine[[i]] <- data.frame(DPS1 = swineData$DPS1, 
                                    DPS2 = swineData$DPS2, 
                                    freq1 = swineData$freq1, 
                                    freq2 = 0)
}

# generate mock dataset
for (i in 1:nrow(swineData)) {
  ini_p <- swineData$freq1[i] # add noise
  gen <- swineData$gen[i]
  v <- 0.004*ini_p*(1-ini_p)
  mock_df_bs_swine <- beta_spike_threshold(MLE_Ne_swineData, ini_p, gen, mu, freq_step, v, alpha)
  for (j in 1:1000) {
    random_num <- runif(1)
    mock_beta_spike_swine <- mock_df_bs_swine[which.min(abs(mock_df_bs_swine$cdf - random_num)),]
    list_df_sim_bs_swine[[j]]$freq2[i] <- mock_beta_spike_swine$freq
  }
}

# calculate log-likelihood
list_ll_max_swine <- 1:length(list_df_sim_bs_swine) %>%
  mclapply(function(x) {
    list_ll_NE_swine <- data.frame()
    mock_ll_df_swine <- LL_Ne.df(list_df_sim_bs_swine[[x]], seq(3, 30, 0.5), mu, alpha, hourPerGen, m)
    list_ll_NE_swine <- rbind(list_ll_NE_swine, mock_ll_df_swine[which.max(mock_ll_df_swine$prob), ])
    return(list_ll_NE_swine)
  }, mc.cores = detectCores())
list_ll_max_swine <- do.call(rbind, list_ll_max_swine)

## 5.4. Sum change of frequency ====
# swine data
sum_dfreq_swine <- swineData %>%
  mutate(sum_dfreq = abs(freq2 - freq1)) %>%
  summarise(total = sum(sum_dfreq, na.rm = TRUE)) %>%
  pull(total)

# simulated data
sum_dfreq_sim_swine <- map_dbl(list_df_sim_bs_swine, function(df) {
  df %>%
    mutate(sum_dfreq = abs(freq2 - freq1)) %>%
    summarise(total = sum(sum_dfreq, na.rm = TRUE)) %>%
    pull(total)
})

## 5.5. Sum change of $(f_2 - f_1) / \sqrt{\Delta t / N_e / f_1 / (1 - f_1)}$ ====
# swine data
sum_z_swine <- sum(abs((swineData$freq2 - swineData$freq1) *
                           sqrt(Ne * swineData$freq1 * (1 - swineData$freq1) / ((swineData$DPS2 - swineData$DPS1)*(24/hourPerGen)))))

# simulated data
Ne_vec_swine <- list_ll_max_swine$Ne
sum_z_sim_swine <- numeric(length(list_df_sim_bs_swine))

for (i in seq_along(list_df_sim_bs_swine)) {
  df <- list_df_sim_bs_swine[[i]]
  Ne <- Ne_vec[i]
  dt <- (df$DPS2 - df$DPS1)*(24/hourPerGen)
  
  Z <- rep(NA_real_, nrow(df))
  Z <- (df$freq2 - df$freq1) *
    sqrt(Ne * df$freq1 * (1 - df$freq1) / dt)
  
  sum_z_sim_swine[i] <- sum(abs(Z), na.rm = TRUE)
}

## 5.6. Proportion of fixation/loss ====
# swine data
prop_fixloss_swine <- swineData %>%
  summarise(
    prop = mean(freq2 <= alpha | freq2 >= (1 - alpha), na.rm = TRUE)
  ) %>%
  pull(prop)

# simulated data
prop_fixloss_sim_swine <- map_dbl(list_df_sim_bs_swine, function(df) {
  mean(df$freq2 <= alpha | df$freq2 >= (1 - alpha), na.rm = TRUE)
})


# 6. Graphs ====
fig4a <- ggplot(swineData_tot_obs, aes(x = DPS, y = freq)) +
  geom_hline(yintercept = 0.02, color = "red", linetype = "dashed", ) +
  geom_point(aes(color = type, shape = type), alpha = 0.5) +
  geom_line(aes(group = interaction(pigID, pair_id)), alpha = 0.3) +
  scale_shape_manual(values = c(19, 17, 2))+
  scale_color_manual(values = c("#26185F", "#0095AF", "black"))+
  labs(x = "Day of fair", y = "iSNV frequency", color = NULL, shape = NULL) +
  theme(legend.position = c(0.3, 0.93)) +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  scale_y_continuous(limits = c(0, 1))

fig4b <- ggplot(LLdf_swineData) +
  geom_line(aes(x = Ne, y = prob)) +
  geom_vline(xintercept = MLE_Ne_swineData, color = "red") +
  geom_hline(yintercept = max(LLdf_swineData$prob), color = "red", linetype = "dashed") +
  annotate("text", x = 25, y = -175, label = paste("N[E] == ", MLE_Ne_swineData), parse = TRUE, color = "red") +
  annotate("rect", xmin = as.numeric(CI_low_Ne_swineData), xmax = as.numeric(CI_high_Ne_swineData), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  labs(x = expression(N[E]), y = "Log-likelihood") +
  scale_x_continuous(breaks = seq(0, 50, 5), limits = c(0, 31))

fig4c <- ggplot(list_ll_max_swine, aes(x = prob)) +
  geom_histogram(bins = 60, fill = "#005D9E") +
  geom_vline(xintercept = max(LLdf_swineData$prob), color = "red", linetype = "dashed") +
  labs(x = "Log-likelihood", y = "Count")

fig4d <- ggplot(list_ll_max_swine, aes(x = Ne)) +
  geom_bar(fill = "#005D9E") +
  geom_vline(xintercept = MLE_Ne_swineData, color = "red") +
  annotate("rect", xmin = as.numeric(CI_low_Ne_swineData), xmax = as.numeric(CI_high_Ne_swineData), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  labs(x = expression(N[E]), y = "Count") +
  scale_x_continuous(breaks = seq(0, 50, 5), limits = c(0, 31))

fig_sum_dfreq_swine <- ggplot(data.frame(sum_dfreq_sim_swine), aes(x = sum_dfreq_sim_swine)) +
  geom_histogram(bins = 60, fill = "#005D9E") +
  geom_vline(xintercept = sum_dfreq_swine, color = "red") +
  labs(x = "Total absolute frequency change", y = "Count")

fig_prop_fixloss_swine <- ggplot(data.frame(prop_fixloss_sim_swine), aes(x = prop_fixloss_sim_swine)) +
  geom_histogram(fill = "#005D9E", bins = 50) +
  geom_vline(xintercept = prop_fixloss_swine, color = "red") +
  labs(x = "Proportion fixed or lost at observation 2", y = "Count")


# 7. Put graphs together ====
fig4 <- fig4a + fig4b + fig4d + fig4c + fig_sum_dfreq_swine + fig_prop_fixloss_swine + 
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 18))

ggsave(fig4, file = "fig4.pdf", path = "02_plots", width = 15, height = 8)
