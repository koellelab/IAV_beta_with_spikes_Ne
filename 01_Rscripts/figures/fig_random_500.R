# 1. Source functions ====
source("01_Rscripts/WithinhostIAV_Ne_fnuctions_log.R")

# 2. Load packages ====
library(tidyverse)
library(ggplot2)
library(patchwork)
library(showtext)
library(parallel)

theme_set(theme_minimal() +
            theme(axis.ticks = element_line(colour = "grey50"),
                  axis.line = element_line(colour = "grey50"),
                  axis.text = element_text(size = 12),
                  axis.title = element_text(size = 15),
                  legend.text = element_text(size = 12)))

# 3. Input parameters and datasets ====
Ne_values <- seq(10, 200, 2)
mu <- 0
alpha <- 0.02
hourPerGen <- 6
m <- 0.004
freq_step <- 0.01
n_iterations <- 500
set.seed(27)

flumina_0.02_origin <- read_csv("03_data/processed_vcf_wide_minAF0.02.csv")[, c(1:6, 28, 29, 31, 33:37)]
swineData <- read_csv("03_data/table_S2.csv") %>%
  filter(Subset == 1) %>%
  set_names(c("pigID", "segment", "DPS1", "DPS2", "freq1", "freq2", "Subset")) %>%
  mutate(gen = (DPS2 - DPS1)*(24/hourPerGen))

# 4. Generate random dataset similar to Table S1 ====
flumina_0.02_origin$index <- 1:nrow(flumina_0.02_origin)
dps_values <- c(-2, 0:6)

# Clean up original data
McCrone_IAV_flumina_data <- flumina_0.02_origin %>%
  filter(data == "mccrone") %>%
  rowwise() %>%
  mutate(DPS1 = dps_values[which(!is.na(c_across(c(`d-2`, d0:d6))))][1],
         DPS2 = dps_values[which(!is.na(c_across(c(`d-2`, d0:d6))))][2], 
         gen = (DPS2 - DPS1)*(24/hourPerGen),
         iSNV = paste0(segment, "_", ref_nt, as.character(iSNV_loc), iSNV_nt),
         non_na_values = list(na.omit(c_across(c(`d-2`, d0:d6)))),
         freq1 = non_na_values[1], 
         freq2 = non_na_values[2]) %>% 
  select(-non_na_values) %>%
  filter(!is.na(freq1), !is.na(freq2)) %>%
  ungroup()

# Create random datasets
list_McCrone_IAV_flumina_random <- map(1:n_iterations, ~generate_random_dataset(McCrone_IAV_flumina_data, .x))
# cleanup random datasets
for (i in 1:n_iterations) {
  list_McCrone_IAV_flumina_random[[i]] <- list_McCrone_IAV_flumina_random[[i]] %>%
    set_names(c("ENROLLID", "mutation", "DPS1", "DPS2", "freq1", "freq2", "Subset")) %>%
    mutate(gen = (DPS2 - DPS1)*(24/hourPerGen)) %>%
    filter(Subset == 1) %>%
    select(DPS1, DPS2, freq1, freq2)
}

# 5. Generate 500 MLE estimations for each random dataset obtained ====
list_ne_est <- 1:length(list_McCrone_IAV_flumina_random) %>%
  mclapply(function(x) {
    ll_df_random <- LL_Ne.df(list_McCrone_IAV_flumina_random[[x]], Ne_values, mu, alpha, hourPerGen, m)
    mle_row <- ll_df_random[which.max(ll_df_random$prob), ]
    return(list(full_ll_df = ll_df_random, mle_row = mle_row))
  }, mc.cores = detectCores())

# extract log-likelihood curve for each random sample
list_ll_random <- map_dfr(list_ne_est, 1, .id = "iteration")
list_ll_random$iteration <- as.factor(list_ll_random$iteration)

# extract MLE for each random sample
list_mle_max_random <- map_dfr(list_ne_est, 2)

# 6. Calculated CIs for 50 random samples ====
# Pre-allocate results
ci_df <- tibble(
  iteration = integer(n_iterations),
  mle_Ne    = numeric(n_iterations),
  ci_low    = numeric(n_iterations),
  ci_high   = numeric(n_iterations)
)

for (i in seq_len(n_iterations)) {
  ll_df <- list_ne_est[[i]]$full_ll_df %>%
    arrange(Ne)
  # MLE
  mle_idx <- which.max(ll_df$prob)
  mle_Ne  <- ll_df$Ne[mle_idx]
  max_ll  <- ll_df$prob[mle_idx]
  target <- max_ll - 1.92
  
  # CI lower
  left_df <- ll_df %>% filter(Ne < mle_Ne)
  ci_low <- if (nrow(left_df) == 0) {
    NA_real_
  } else {
    left_df$Ne[which.min(abs(left_df$prob - target))]
  }
  
  # CI upper
  right_df <- ll_df %>% filter(Ne > mle_Ne)
  ci_high <- if (nrow(right_df) == 0) {
    NA_real_
  } else {
    right_df$Ne[which.min(abs(right_df$prob - target))]
  }
  
  # Store results
  ci_df$iteration[i] <- i
  ci_df$mle_Ne[i]    <- mle_Ne
  ci_df$ci_low[i]    <- ci_low
  ci_df$ci_high[i]   <- ci_high
}

# randomly select samples to be plotted
ci_df_sample <- ci_df %>% slice_sample(n = 50)

# add the 50% method value to the dataset
close_to_50_row <- tibble(
  iteration = 0,   # placeholder; we will override x position below
  mle_Ne    = as.numeric(MLE_Ne_McCrone),
  ci_low    = as.numeric(CI_low_Ne_McCrone),
  ci_high   = as.numeric(CI_high_Ne_McCrone)
)

# adjust the df to be better for plotting
gap_size <- 5
true_x   <- 1
rand_x0  <- true_x + gap_size

rand_df <- ci_df_sample %>%
  mutate(x_pos = rand_x0 + row_number() - 1)

close_to_50_df <- close_to_50_row %>%
  mutate(x_pos = true_x)

rand_mid <- (min(rand_df$x_pos) + max(rand_df$x_pos)) / 2

# 7. Calculated if CI for 50% method is smaller than random samples ====
ci_df$ci_range <- ci_df$ci_high - ci_df$ci_low
percent50_CI <- as.numeric(CI_high_Ne_McCrone) - as.numeric(CI_low_Ne_McCrone)

percent_CI_smaller <- sum(percent50_CI < ci_df$ci_range) / nrow(ci_df)

# 8. Re-do everything with swine data ====
Ne_values <- seq(4, 80, 1)

swine_flumina_data <- flumina_0.02_origin %>%
  filter(data == "swine") %>%
  rowwise() %>%
  mutate(DPS1 = dps_values[which(!is.na(c_across(c(`d-2`, d0:d6))))][1],
         DPS2 = dps_values[which(!is.na(c_across(c(`d-2`, d0:d6))))][2], 
         gen = (DPS2 - DPS1)*(24/hourPerGen),
         iSNV = paste0(segment, "_", ref_nt, as.character(iSNV_loc), iSNV_nt),
         non_na_values = list(na.omit(c_across(c(`d-2`, d0:d6)))),
         freq1 = non_na_values[1], 
         freq2 = non_na_values[2]) %>% 
  select(-non_na_values) %>%
  ungroup()

LLdf_swineData_Ne80 <- LL_Ne.df(swineData, Ne_values, mu, alpha, hourPerGen, m)
MLE_Ne_swineData_Ne80 <- LLdf_swineData_Ne80$Ne[which.max(LLdf_swineData_Ne80$prob)]
CI_low_Ne_swineData_Ne80 <- LLdf_swineData_Ne80 %>%
  filter(Ne < MLE_Ne_swineData_Ne80) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_swineData_Ne80 <- LLdf_swineData_Ne80 %>%
  filter(Ne > MLE_Ne_swineData_Ne80) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)


# Create random datasets
list_swine_flumina_random <- map(1:n_iterations, ~generate_random_dataset(swine_flumina_data, .x))
# cleanup random datasets
for (i in 1:n_iterations) {
  list_swine_flumina_random[[i]] <- list_swine_flumina_random[[i]] %>%
    set_names(c("ENROLLID", "mutation", "DPS1", "DPS2", "freq1", "freq2", "Subset")) %>%
    mutate(gen = (DPS2 - DPS1)*(24/hourPerGen)) %>%
    filter(Subset == 1) %>%
    select(DPS1, DPS2, freq1, freq2)
}

## 8.1. Generate 500 MLE estimations for each random dataset obtained ====
list_ne_est_swine <- 1:n_iterations %>%
  mclapply(function(x) {
    ll_df_random <- LL_Ne.df(list_swine_flumina_random[[x]], Ne_values, mu, alpha, hourPerGen, m)
    mle_row <- ll_df_random[which.max(ll_df_random$prob), ]
    return(list(full_ll_df = ll_df_random, mle_row = mle_row))
  }, mc.cores = detectCores())

# extract log-likelihood curve for each random sample
list_ll_random_swine <- map_dfr(list_ne_est_swine, 1, .id = "iteration")
list_ll_random_swine$iteration <- as.factor(list_ll_random_swine$iteration)

# extract MLE for each random sample
list_mle_max_random_swine <- map_dfr(list_ne_est_swine, 2)

## 8.2. Calculated CIs for 100 random samples ====
# Pre-allocate results
ci_df_swine <- tibble(
  iteration = integer(n_iterations),
  mle_Ne    = numeric(n_iterations),
  ci_low    = numeric(n_iterations),
  ci_high   = numeric(n_iterations)
)

for (i in seq_len(n_iterations)) {
  ll_df <- list_ne_est_swine[[i]]$full_ll_df %>%
    arrange(Ne)
  # MLE
  mle_idx <- which.max(ll_df$prob)
  mle_Ne  <- ll_df$Ne[mle_idx]
  max_ll  <- ll_df$prob[mle_idx]
  target <- max_ll - 1.92
  
  # CI lower
  left_df <- ll_df %>% filter(Ne < mle_Ne)
  ci_low <- if (nrow(left_df) == 0) {
    NA_real_
  } else {
    left_df$Ne[which.min(abs(left_df$prob - target))]
  }
  
  # CI upper
  right_df <- ll_df %>% filter(Ne > mle_Ne)
  ci_high <- if (nrow(right_df) == 0) {
    NA_real_
  } else {
    right_df$Ne[which.min(abs(right_df$prob - target))]
  }
  
  # Store results
  ci_df_swine$iteration[i] <- i
  ci_df_swine$mle_Ne[i]    <- mle_Ne
  ci_df_swine$ci_low[i]    <- ci_low
  ci_df_swine$ci_high[i]   <- ci_high
}

# randomly select samples to be plotted
ci_df_sample_swine <- ci_df_swine %>% slice_sample(n = 50)

# add the 50% method value to the dataset
close_to_50_row_swine <- tibble(
  iteration = 0,   # placeholder; we will override x position below
  mle_Ne    = as.numeric(MLE_Ne_swineData),
  ci_low    = as.numeric(CI_low_Ne_swineData),
  ci_high   = as.numeric(CI_high_Ne_swineData)
)

# adjust the df to be better for plotting
gap_size <- 5
true_x   <- 1
rand_x0  <- true_x + gap_size

rand_df_swine <- ci_df_sample_swine %>%
  mutate(x_pos = rand_x0 + row_number() - 1)

close_to_50_df_swine <- close_to_50_row_swine %>%
  mutate(x_pos = true_x)

rand_mid_swine <- (min(rand_df_swine$x_pos) + max(rand_df_swine$x_pos)) / 2

## 8.3. Calculated if CI for 50% method is smaller than random samples ====
ci_df_swine$ci_range <- ci_df_swine$ci_high - ci_df_swine$ci_low
percent50_CI_swine <- as.numeric(CI_high_Ne_swineData) - as.numeric(CI_low_Ne_swineData)

percent_CI_smaller_swine <- sum(percent50_CI_swine < ci_df_swine$ci_range) / nrow(ci_df_swine)

# 9. Graphs ====
fig_traj_random_human <- ggplot() +
  geom_line(data = list_ll_random, aes(x = Ne, y = prob, group = iteration), alpha = 0.1, color = "darkgrey") + 
  geom_point(data = list_mle_max_random, aes(x = Ne, y = prob), color = "#0095AF", size = 1, alpha = 0.3) +
  geom_line(data = LLdf_McCrone, aes(x = Ne, y = prob), color = "#005D9E", linewidth = 1) +
  geom_vline(xintercept = MLE_Ne_McCrone, color = "red") +
  labs(x = expression(N[E]), y = "Log-likelihood") +
  scale_x_continuous(breaks = seq(0, 200, 20), limits = c(10, 150))

fig_CI_sample_human <- ggplot() +
  geom_errorbar(data = rand_df,
                aes(x = x_pos, ymin = ci_low, ymax = ci_high),
                width = 1, color = "darkgrey") +
  geom_point(data = rand_df,
             aes(x = x_pos, y = mle_Ne),
             size = 1.8, color = "#0095AF", alpha = 0.6) +
  geom_errorbar(data = close_to_50_df,
                aes(x = x_pos, ymin = ci_low, ymax = ci_high),
                width = 1, color = "#005D9E") +
  geom_point(data = close_to_50_df,
             aes(x = x_pos, y = mle_Ne),
             size = 2.2, color = "red") +
  scale_x_continuous(
    breaks = c(true_x, rand_mid),
    labels = c("closest-to-50%", "random samples")
  ) +
  labs(x = NULL, y = expression("MLE of " * N[E] * " with 95% CI")) 

fig_traj_random_swine <- ggplot() +
  geom_line(data = list_ll_random_swine, aes(x = Ne, y = prob, group = iteration), alpha = 0.1, color = "darkgrey") + 
  geom_point(data = list_mle_max_random_swine, aes(x = Ne, y = prob), color = "#0095AF", size = 1, alpha = 0.3) +
  geom_line(data = LLdf_swineData_Ne80, aes(x = Ne, y = prob), color = "#005D9E", linewidth = 1) +
  geom_vline(xintercept = MLE_Ne_swineData_Ne80, color = "red") +
  labs(x = expression(N[E]), y = "Log-likelihood")
  #scale_x_continuous(breaks = seq(0, 200, 20), limits = c(10, 100))

fig_CI_sample_swine <- ggplot() +
  geom_errorbar(data = rand_df_swine,
                aes(x = x_pos, ymin = ci_low, ymax = ci_high),
                width = 1, color = "darkgrey") +
  geom_point(data = rand_df_swine,
             aes(x = x_pos, y = mle_Ne),
             size = 1.8, color = "#0095AF", alpha = 0.6) +
  geom_errorbar(data = close_to_50_df_swine,
                aes(x = x_pos, ymin = ci_low, ymax = ci_high),
                width = 1, color = "#005D9E") +
  geom_point(data = close_to_50_df_swine,
             aes(x = x_pos, y = mle_Ne),
             size = 2.2, color = "red") +
  scale_x_continuous(
    breaks = c(true_x, rand_mid_swine),
    labels = c("closest-to-50% method", "random samples")
  ) +
  labs(x = NULL, y = expression("MLE of " * N[E] * " with 95% CI"))


# 10. Put panels together ====
fig_random_human <- fig_traj_random_human + fig_CI_sample_human + 
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 18),
        plot.title = element_text(face = "bold", size = 18))

fig_random_swine <- fig_traj_random_swine + fig_CI_sample_swine + 
  plot_layout(ncol = 2) +
  plot_annotation(title = "Swine data subset 1", tag_levels = list(c("C", "D"))) &
  theme(plot.tag = element_text(face = "bold", size = 18),
        plot.title = element_text(face = "bold", size = 18))

fig_random <- wrap_elements(fig_random_human) / wrap_elements(fig_random_swine) 

ggsave(fig_random, file = "fig_random.pdf", path = "02_plots", width = 10, height = 8)
ggsave(fig_random_human, file = "fig_random_human.pdf", path = "02_plots", width = 10, height = 4)
