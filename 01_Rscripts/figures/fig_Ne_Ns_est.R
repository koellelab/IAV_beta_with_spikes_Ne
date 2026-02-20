# 1. Source functions ====
source("01_Rscripts/WithinhostIAV_Ne_fnuctions_log.R")

# 2. Load packages ====
library(tidyverse)
library(ggplot2)
library(patchwork)
library(showtext)
library(colorspace)
theme_set(theme_minimal() +
            theme(axis.ticks = element_line(colour = "grey50"),
                  axis.line = element_line(colour = "grey50"),
                  axis.text = element_text(size = 12),
                  axis.title = element_text(size = 15),
                  legend.text = element_text(size = 12)))

# 3. Input parameters ====
mu <- 0
alpha <- 0.02
hourPerGen <- 6
m <- 0.004

# 4. Input datasets ====
#read_csv("03_data/table_S1_rep.csv") %>%
McCrone_IAV <- read_csv("03_data/table_S1.csv") %>%
  set_names(c("ENROLLID", "mutation", "DPS1", "DPS2", "freq1", "freq2", "Subset")) %>%
  mutate(gen = (DPS2 - DPS1) * (24 / hourPerGen)) %>%
  filter(Subset == 1)

swineData <- read_csv("03_data/table_S2.csv") %>%
  filter(Subset == 1) %>%
  set_names(c("pigID", "segment", "DPS1", "DPS2", "freq1", "freq2", "Subset")) %>%
  mutate(gen = (DPS2 - DPS1)*(24/hourPerGen))

# 5. Generate data frame ====
## 5.1. Calculate Ne/Ns MLE ====
### 5.1.1. human ====
NeNs_surface_df_McCrone <- LL_NeNs_surface(
  df = McCrone_IAV,
  Ne_values = seq(3, 200, 2),
  Ns_values = seq(1, 1000, 10),
  hourPerGen = hourPerGen
)

mle_row_McCrone <- NeNs_surface_df_McCrone %>%
  slice(which.max(logLik))
MLE_Ne_McCrone_joint <- mle_row_McCrone$Ne
MLE_Ns_McCrone_joint <- mle_row_McCrone$Ns
max_logLik_McCrone <- mle_row_McCrone$logLik

CI_logLik_McCrone <- max_logLik_McCrone - 2.995

### 5.1.2. swine ====
NeNs_surface_df_swine <- LL_NeNs_surface(
  df = swineData,
  Ne_values = seq(3, 500, 1),
  Ns_values = seq(1, 7.5, 0.05),
  hourPerGen = hourPerGen
)

mle_row_swine <- NeNs_surface_df_swine %>%
  slice(which.max(logLik))
MLE_Ne_swine_joint <- mle_row_swine$Ne
MLE_Ns_swine_joint <- mle_row_swine$Ns
max_logLik_swine <- mle_row_swine$logLik

CI_logLik_swine <- max_logLik_swine - 2.995

## 5.2. Calculate for multiple simulated datasets ====
mle_row_McCrone_df <- data.frame()
list_NeNs_surface_df_sim <- vector("list", length = 500)

for (i in 1:500) {
  NeNs_surface_df_sim <- LL_NeNs_surface(
    df = list_df_sim_bs[[i]],
    Ne_values = seq(1, 200, 5),
    Ns_values = seq(1, 200, 5),
    hourPerGen = hourPerGen
  )
  list_NeNs_surface_df_sim[[i]] <- NeNs_surface_df_sim
  
  max_row <- NeNs_surface_df_sim %>%
    slice(which.max(logLik))
  mle_row_McCrone_df <- rbind(mle_row_McCrone_df, max_row)
}

## 5.3. Find local and global MLEs for samples ====
peaks_df <- data.frame()

for (i in seq_along(list_NeNs_surface_df_sim)) {
  df_i <- list_NeNs_surface_df_sim[[i]] %>%
    select(Ne, Ns, logLik) %>%
    mutate(sample = i)
  
  # Global peak (single per sample)
  global_row <- df_i %>%
    slice(which.max(logLik)) %>%
    mutate(type = "global")
  
  # Local peaks per Ns (in Ne direction)
  local_rows_allNs <- data.frame()
  Ns_vals <- sort(unique(df_i$Ns))
  for (ns in Ns_vals) {
    df_ns <- df_i %>%
      filter(Ns == ns) %>%
      arrange(Ne)
    
    # Need at least 3 Ne points to define a local interior peak
    if (nrow(df_ns) < 3) next
    
    # Check strict local maxima: logLik[j] > logLik[j-1] AND > logLik[j+1]
    ll <- df_ns$logLik
    is_local <- rep(FALSE, length(ll))
    is_local[2:(length(ll) - 1)] <- (ll[2:(length(ll) - 1)] > ll[1:(length(ll) - 2)]) &
      (ll[2:(length(ll) - 1)] > ll[3:length(ll)])
    
    local_rows <- df_ns[is_local, , drop = FALSE]
    
    if (nrow(local_rows) > 0) {
      local_rows <- local_rows %>% filter(Ne != global_row$Ne)
      local_rows$type <- "local"
      local_rows_allNs <- rbind(local_rows_allNs, local_rows)
    }
  }
  
  # Combine for this sample
  peaks_i <- rbind(global_row, local_rows_allNs)
  peaks_df <- rbind(peaks_df, peaks_i)
} 

peaks_df$type <- factor(peaks_df$type, levels = c("global", "local"))

peaks_df <- peaks_df %>%
  arrange(sample, Ne, Ns, logLik, type) %>%
  mutate(type = as.character(type))


## 5.4. Find local MLEs for McCrone and swine data ====
# human
local_Ns_McCrone <- data.frame()
Ns_vals_McCrone <- sort(unique(NeNs_surface_df_McCrone$Ns))
for (ns in Ns_vals) {
  df_ns <- NeNs_surface_df_McCrone %>%
    filter(Ns == ns) %>%
    arrange(Ne)
  
  # Need at least 3 Ne points to define a local interior peak
  if (nrow(df_ns) < 3) next
  
  # Check strict local maxima: logLik[j] > logLik[j-1] AND > logLik[j+1]
  ll <- df_ns$logLik
  is_local <- rep(FALSE, length(ll))
  is_local[2:(length(ll) - 1)] <- (ll[2:(length(ll) - 1)] > ll[1:(length(ll) - 2)]) &
    (ll[2:(length(ll) - 1)] > ll[3:length(ll)])
  
  local_rows <- df_ns[is_local, , drop = FALSE]
  
  if (nrow(local_rows) > 0) {
    local_rows <- local_rows %>% filter(Ne != mle_row_McCrone$Ne)
    local_rows$type <- "local"
    local_Ns_McCrone <- rbind(local_Ns_McCrone, local_rows)
  }
}

#swine
local_Ns_swine <- data.frame()
Ns_vals_swine <- sort(unique(NeNs_surface_df_swine$Ns))
for (ns in Ns_vals) {
  df_ns <- NeNs_surface_df_swine %>%
    filter(Ns == ns) %>%
    arrange(Ne)
  
  # Need at least 3 Ne points to define a local interior peak
  if (nrow(df_ns) < 3) next
  
  # Check strict local maxima: logLik[j] > logLik[j-1] AND > logLik[j+1]
  ll <- df_ns$logLik
  is_local <- rep(FALSE, length(ll))
  is_local[2:(length(ll) - 1)] <- (ll[2:(length(ll) - 1)] > ll[1:(length(ll) - 2)]) &
    (ll[2:(length(ll) - 1)] > ll[3:length(ll)])
  
  local_rows <- df_ns[is_local, , drop = FALSE]
  
  if (nrow(local_rows) > 0) {
    local_rows <- local_rows %>% filter(Ne != mle_row_swine$Ne)
    local_rows$type <- "local"
    local_Ns_swine <- rbind(local_Ns_swine, local_rows)
  }
}

## 5.5. Estimate Ns using only human data replicates ====
McCrone_IAV_rep <- read_csv("03_data/table_S1_rep.csv") %>%
  set_names(c("ENROLLID", "mutation", "DPS1", "DPS2", "freq1", "freq2", "Subset")) %>%
  mutate(gen = (DPS2 - DPS1) * (24 / hourPerGen)) %>%
  filter(Subset == 1) %>%
  filter(gen == 1) %>%
  mutate(DPS2 = DPS1,
         gen = 0)

NeNs_surface_df_McCrone_rep <- LL_NeNs_surface(
  df = McCrone_IAV_rep,
  Ne_values = seq(3, 200, 2),
  Ns_values = seq(1, 1000, 10),
  hourPerGen = hourPerGen
)

mle_row_McCrone_rep <- NeNs_surface_df_McCrone_rep %>%
  slice(which.max(logLik))
MLE_Ne_McCrone_joint_rep <- mle_row_McCrone_rep$Ne
MLE_Ns_McCrone_joint_rep <- mle_row_McCrone_rep$Ns
max_logLik_McCrone_rep <- mle_row_McCrone_rep$logLik

CI_logLik_McCrone_rep <- max_logLik_McCrone_rep - 2.995

## 5.6. Generate expected values for scaled absolute change ====
df_expected_dq_scaled <- df_expected_dq_scaled %>%
  mutate(exp_q_scaled_71 = sqrt(2/pi) * sqrt(gen/71 + 1/141),
         exp_q_scaled_1500 = sqrt(2/pi) * sqrt(gen/1500 + 1/3.1))
         

# 6. Graphs ====
## 6.1.  human ====
NeNs_surface_df_McCrone_plot <- NeNs_surface_df_McCrone %>%
  mutate(logLik = ifelse(logLik < CI_logLik_McCrone, NA, logLik))

fig_NeNs_McCrone <- ggplot() +
  geom_tile(data = NeNs_surface_df_McCrone_plot, aes(x = Ne, y = Ns, fill = logLik)) +
  geom_vline(xintercept = 49, color = "red") +
  #geom_line(data = local_Ns_McCrone, aes(x = Ne, y = Ns), color = "red") +
  geom_vline(xintercept = 71, color = "red", linetype = "dashed") +
  geom_hline(yintercept = MLE_Ns_McCrone_joint_rep, color = "red", linetype = "dashed") +
  annotate("text", x = 150, y = 100, 
           #label = paste0("Ne = ", MLE_Ne_McCrone_joint, ", Ns = ", MLE_Ns_McCrone_joint_rep), 
           label = paste0("Ne = 71, Ns = ", MLE_Ns_McCrone_joint_rep), 
           color = "red") +
  labs(
    x = expression(N[E]),
    y = expression(N[S]),
    fill = "Log-likelihood"
  ) +
  scale_fill_continuous_sequential(palette = "YlGnBu", rev = FALSE, begin = 0, end = 0.6,
                                   na.value = "transparent")

fig_dfreq_scaled_McCrone_71 <- ggplot() +
  geom_point(data = df_scatter_McCrone, 
             aes(x = days, 
                 y = abs(delta_q_scaled), 
                 shape = above_threshold,
                 color = above_threshold), 
             alpha = 0.6, size = 1.5) +
  scale_shape_manual(values = c(19, 1)) +
  scale_color_manual(values = c("#26185F","#0095AF")) +
  geom_line(data = df_expected_dq_scaled, aes(x = days, y = exp_q_scaled_McCrone), color = "red") +
  geom_line(data = df_expected_dq_scaled, aes(x = days, y = exp_q_scaled_71), color = "orange") +
  geom_line(data = mean_scaled_dfreq_McCrone[1:4,], aes(x = days, y = mean), 
            color = "darkgrey", linewidth = 2, alpha = 0.5) +
  geom_point(data = mean_scaled_dfreq_McCrone[1:4,], aes(x = days, y = mean), shape = 4, size = 3) +
  scale_x_continuous(breaks = seq(0, 100, 1)) +
  labs(x = "Days between observations",
       y = expression(abs(Delta*q) / sqrt(q[0] * (1 - q[0]))),
       color = "Second observation\nrelative to threshold",
       shape = "Second observation\nrelative to threshold")

## 6.2. swine ====
NeNs_surface_df_swine_plot <- NeNs_surface_df_swine %>%
  mutate(logLik = ifelse(logLik < CI_logLik_swine, NA, logLik))

fig_NeNs_swine <- ggplot() +
  geom_tile(data = NeNs_surface_df_swine_plot, aes(x = Ne, y = Ns, fill = logLik)) +
  geom_vline(xintercept = 14.4, color = "red") +
  #geom_line(data = local_Ns_swine, aes(x = Ne, y = Ns), color = "red") +
  #geom_vline(xintercept = MLE_Ne_swine_joint, color = "red", linetype = "dashed") +
  geom_hline(yintercept = MLE_Ns_swine_joint, color = "red", linetype = "dashed") +
  annotate("text", x = 200, y = 1.5, label = paste0("Ns = ", MLE_Ns_swine_joint), color = "red") +
  labs(
    x = expression(N[E]),
    y = expression(N[S]),
    fill = "Log-likelihood"
  ) +
  scale_fill_continuous_sequential(palette = "YlGnBu", rev = FALSE, begin = 0, end = 0.6,
                                   na.value = "transparent")

fig_dfreq_scaled_swine_1000 <- ggplot() +
  geom_point(data = df_scatter_swine, 
             aes(x = days, 
                 y = abs(delta_q_scaled), 
                 shape = above_threshold,
                 color = above_threshold),
             alpha = 0.6, size = 1.5) +
  scale_shape_manual(values = c(19, 1)) +
  scale_color_manual(values = c("#26185F","#0095AF")) +
  scale_x_continuous(breaks = seq(0, 100, 1)) +
  scale_y_continuous(breaks = 1:5, labels = scales::number_format(accuracy = 0.1)) +
  geom_line(data = df_expected_dq_scaled, aes(x = days, y = exp_q_scaled_swine), color = "red") +
  geom_line(data = df_expected_dq_scaled, aes(x = days, y = exp_q_scaled_1500), color = "orange") +
  geom_line(data = mean_scaled_dfreq_swine, aes(x = days, y = mean), 
            color = "darkgrey", linewidth = 2, alpha = 0.5) +
  geom_point(data = mean_scaled_dfreq_swine, aes(x = days, y = mean), shape = 4, size = 3) +
  labs(x = "Days between observations",
       y = expression(abs(Delta*q) / sqrt(q[0] * (1 - q[0]))),
       color = "Second observation\nrelative to threshold",
       shape = "Second observation\nrelative to threshold")


# 7. Save ====
fig_Ne_Ns_est_human <- fig_NeNs_McCrone + fig_dfreq_scaled_McCrone_71 + 
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.title = element_text(face = "bold", size = 18))) &
  theme(plot.tag = element_text(face = "bold", size = 18))

fig_Ne_Ns_est_swine <- fig_NeNs_swine + fig_dfreq_scaled_swine_1000 + 
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.title = element_text(face = "bold", size = 18))) &
  theme(plot.tag = element_text(face = "bold", size = 18))

fig_Ne_Ns_est <- wrap_elements(fig_Ne_Ns_est_human) / wrap_elements(fig_Ne_Ns_est_swine)

#fig_Ne_Ns_est <- fig_NeNs_McCrone + fig_dfreq_scaled_McCrone_150 + fig_NeNs_swine + fig_dfreq_scaled_swine_1000 +
#  plot_layout(ncol = 2, guides = "collect") +
#  plot_annotation(tag_levels = "A") &
#  theme(legend.position = "bottom", plot.tag = element_text(face = "bold", size = 18))

ggsave(fig_Ne_Ns_est_human, file = "fig_Ne_Ns_heatmap_norep_human.pdf", path = "02_plots", width = 12, height = 4)
ggsave(fig_Ne_Ns_est_swine, file = "fig_Ne_Ns_heatmap_norep_swine.pdf", path = "02_plots", width = 12, height = 4)
