# 1. Source functions ====
source("01_Rscripts/WithinhostIAV_Ne_fnuctions_log.R")

# 2. Load packages ====
library(tidyverse)
library(ggplot2)
library(patchwork)
library(showtext)
library(colorspace)
library(parallel)
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

# 4. Input datasets ====
McCrone_IAV <- read_csv("03_data/table_S1.csv") %>%
  set_names(c("ENROLLID", "mutation", "DPS1", "DPS2", "freq1", "freq2", "Subset")) %>%
  mutate(gen = (DPS2 - DPS1)*(24/hourPerGen)) %>%
  filter(Subset == 1)

# 5. Generate data frame ====
## 5.1. Figure 2A ====
# data frame only including first round of observations
McCrone_obs1 <- McCrone_IAV %>%
  select(mutation, ENROLLID, DPS1, freq1, gen)
McCrone_obs1$type <- c("obs1")
names(McCrone_obs1) <- c("mutation", "ENROLLID", "DPS", "freq", "gen", "type")
# data frame only including second round of observations
McCrone_obs2 <- McCrone_IAV %>%
  select(mutation, ENROLLID, DPS2, freq2, gen)
McCrone_obs2$type <- c("obs2")
names(McCrone_obs2) <- c("mutation", "ENROLLID", "DPS", "freq", "gen", "type")
# merge the data together
McCrone_tot_obs <- rbind(McCrone_obs1, McCrone_obs2)
McCrone_tot_obs[McCrone_tot_obs$freq < 0.02, "type"] <- "observation 2 below 2%"
McCrone_tot_obs[McCrone_tot_obs$freq < 0.02, "freq"] <- 0.02
McCrone_tot_obs$type <- factor(McCrone_tot_obs$type, labels = c("observation 1", "observation 2", "observation 2 below 2%"))

## 5.2. Figure 2B ====
LLdf_McCrone <- LL_Ne.df(McCrone_IAV, Ne_values, mu, alpha, hourPerGen, m)
MLE_Ne_McCrone <- LLdf_McCrone$Ne[which.max(LLdf_McCrone$prob)]
CI_low_Ne_McCrone <- LLdf_McCrone %>%
  filter(Ne < MLE_Ne_McCrone) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_McCrone <- LLdf_McCrone %>%
  filter(Ne > MLE_Ne_McCrone) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

## 5.3. Figure 2C ====
# create empty list
list_df_sim_bs <- list()
for (i in 1:1000) {
  list_df_sim_bs[[i]] <- data.frame(DPS1 = McCrone_IAV$DPS1, 
                                    DPS2 = McCrone_IAV$DPS2, 
                                    freq1 = McCrone_IAV$freq1, 
                                    freq2 = 0)
}

# generate mock dataset
for (i in 1:nrow(McCrone_IAV)) {
  ini_p <- McCrone_IAV$freq1[i] 
  gen <- McCrone_IAV$gen[i]
  v <- 0.004*ini_p*(1-ini_p)
  mock_df_bs <- beta_spike_threshold(MLE_Ne_McCrone, ini_p, gen, mu, freq_step, v, alpha)
  for (j in 1:1000) {
    random_num <- runif(1)
    mock_beta_spike <- mock_df_bs[which.min(abs(mock_df_bs$cdf - random_num)),]
    list_df_sim_bs[[j]]$freq2[i] <- mock_beta_spike$freq
  }
}

# calculate log-likelihood
list_ll_max <- 1:length(list_df_sim_bs) %>%
  mclapply(function(x) {
    list_ll_NE <- data.frame()
    mock_ll_df <- LL_Ne.df(list_df_sim_bs[[x]], Ne_values, mu, alpha, hourPerGen, m)
    list_ll_NE <- rbind(list_ll_NE, mock_ll_df[which.max(mock_ll_df$prob), ])
    return(list_ll_NE)
  }, mc.cores = detectCores())
list_ll_max <- do.call(rbind, list_ll_max)

## 5.4. Sum change of frequency ====
# McCrone data
sum_dfreq_McCrone <- McCrone_IAV %>%
  mutate(sum_dfreq = abs(freq2 - freq1)) %>%
  summarise(total = sum(sum_dfreq, na.rm = TRUE)) %>%
  pull(total)

# simulated data
sum_dfreq_sim <- map_dbl(list_df_sim_bs, function(df) {
  df %>%
    mutate(sum_dfreq = abs(freq2 - freq1)) %>%
    summarise(total = sum(sum_dfreq, na.rm = TRUE)) %>%
    pull(total)
})

## 5.5. Sum change of $(f_2 - f_1) / \sqrt{\Delta t / N_e / f_1 / (1 - f_1)}$ ====
# McCrone data
sum_z_McCrone <- sum(abs((McCrone_IAV$freq2 - McCrone_IAV$freq1) *
                           sqrt(Ne * McCrone_IAV$freq1 * (1 - McCrone_IAV$freq1) / ((McCrone_IAV$DPS2 - McCrone_IAV$DPS1)*(24/hourPerGen)))))

# simulated data
Ne_vec <- list_ll_max$Ne
sum_z_sim <- numeric(length(list_df_sim_bs))

for (i in seq_along(list_df_sim_bs)) {
  df <- list_df_sim_bs[[i]]
  Ne <- Ne_vec[i]
  dt <- (df$DPS2 - df$DPS1)*(24/hourPerGen)
  
  Z <- rep(NA_real_, nrow(df))
  Z <- (df$freq2 - df$freq1) *
    sqrt(Ne * df$freq1 * (1 - df$freq1) / dt)
  
  sum_z_sim[i] <- sum(abs(Z), na.rm = TRUE)
}

## 5.6. Proportion of fixation/loss ====
# McCrone data
prop_fixloss_McCrone <- McCrone_IAV %>%
  summarise(
    prop = mean(freq2 <= alpha | freq2 >= (1 - alpha), na.rm = TRUE)
  ) %>%
  pull(prop)

# simulated data
prop_fixloss_sim <- map_dbl(list_df_sim_bs, function(df) {
  mean(df$freq2 <= alpha | df$freq2 >= (1 - alpha), na.rm = TRUE)
})


# 6. Graphs ====
fig2a <- ggplot(McCrone_tot_obs, aes(x = DPS, y = freq)) +
  geom_hline(yintercept = 0.02, color = "red", linetype = "dashed", ) +
  geom_point(aes(color = type, shape = type), alpha = 0.5) +
  geom_line(aes(group = ENROLLID), alpha = 0.3) +
  scale_shape_manual(values = c(19, 17, 2))+
  scale_color_manual(values = c("#26185F", "#0095AF", "black"))+
  labs(x = "Days post symptom onset", y = "iSNV frequency", color = NULL, shape = NULL) +
  theme(legend.position = c(0.3, 0.93)) +
  ylim(0, 0.7)

fig2b <- ggplot(LLdf_McCrone) +
  geom_line(aes(x = Ne, y = prob)) +
  geom_vline(xintercept = MLE_Ne_McCrone, color = "red") +
  geom_hline(yintercept = max(LLdf_McCrone$prob), color = "red", linetype = "dashed") +
  annotate("text", x = 70, y = -40, label = paste("N[E] == ", MLE_Ne_McCrone), parse = TRUE, color = "red") +
  annotate("rect", xmin = as.numeric(CI_low_Ne_McCrone), xmax = as.numeric(CI_high_Ne_McCrone), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  labs(x = expression(N[E]), y = "Log-likelihood") +
  scale_x_continuous(breaks = seq(0, 200, 20), limits = c(0, 100))

fig2c <- ggplot(list_ll_max, aes(x = prob)) +
  geom_histogram(bins = 60, fill = "#005D9E") +
  geom_vline(xintercept = max(LLdf_McCrone$prob), color = "red", linetype = "dashed") +
  labs(x = "Log-likelihood", y = "Count")

fig2d <- ggplot(list_ll_max, aes(x = Ne)) +
  geom_bar(fill = "#005D9E") +
  geom_vline(xintercept = MLE_Ne_McCrone, color = "red") +
  annotate("rect", xmin = as.numeric(CI_low_Ne_McCrone), xmax = as.numeric(CI_high_Ne_McCrone), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  labs(x = expression(N[E]), y = "Count") +
  scale_x_continuous(breaks = seq(0, 200, 20), limits = c(0, 100))

fig_sum_dfreq <- ggplot(data.frame(sum_dfreq_sim), aes(x = sum_dfreq_sim)) +
  geom_histogram(bins = 60, fill = "#005D9E") +
  geom_vline(xintercept = sum_dfreq_McCrone, color = "red") +
  labs(x = "Total absolute frequency change", y = "Count")

fig_prop_fixloss <- ggplot(data.frame(prop_fixloss_sim), aes(x = prop_fixloss_sim)) +
  geom_histogram(fill = "#005D9E", bins = 30) +
  geom_vline(xintercept = prop_fixloss_McCrone, color = "red") +
  labs(x = "Proportion fixed or lost at observation 2", y = "Count")

# 7. Put graphs together ====
fig2 <- fig2a + fig2b + fig2d + fig2c + fig_sum_dfreq + fig_prop_fixloss + 
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 18))

ggsave(fig2, file = "fig_human_ds1.pdf", path = "02_plots", width = 15, height = 8)
