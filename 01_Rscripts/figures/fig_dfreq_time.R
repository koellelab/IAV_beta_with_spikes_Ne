# 1. Source functions ====
source("01_Rscripts/WithinhostIAV_Ne_fnuctions_log.R")

# 2. Load packages ====
library(tidyverse)
library(ggplot2)
library(patchwork)
library(showtext)
library(colorspace)
library(purrr)
theme_set(theme_minimal() +
            theme(axis.ticks = element_line(colour = "grey50"),
                  axis.line = element_line(colour = "grey50"),
                  axis.text = element_text(size = 12),
                  axis.title = element_text(size = 15),
                  legend.title = element_text(size = 13),
                  legend.text = element_text(size = 12)))

# 3. Input parameters ====
mu <- 0
alpha <- 0.02
hourPerGen <- 6
m <- 0.004

# 4. Input datasets ====
McCrone_IAV <- read_csv("03_data/table_S1_rep.csv") %>%
  set_names(c("ENROLLID", "mutation", "DPS1", "DPS2", "freq1", "freq2", "Subset")) %>%
  mutate(gen = (DPS2 - DPS1)*(24/hourPerGen)) %>%
  filter(Subset == 1)

swineData <- read_csv("03_data/table_S2.csv") %>%
  filter(Subset == 1) %>%
  set_names(c("pigID", "segment", "DPS1", "DPS2", "freq1", "freq2", "Subset")) %>%
  mutate(gen = (DPS2 - DPS1)*(24/hourPerGen))

McCrone_q0 <- mean(McCrone_IAV$freq1)
swine_q0 <- mean(swineData$freq1)

# 5. Generate data frame ====
df_scatter_McCrone <- McCrone_IAV %>%
  filter(freq1 > 0) %>%
  mutate(
    days = (DPS2 - DPS1), 
    gen = (DPS2 - DPS1)*(24/hourPerGen),
    q0 = freq1,
    delta_q = freq2 - freq1,
    delta_q_scaled = delta_q / sqrt(q0 * (1 - q0)),
    above_threshold = ifelse(freq2 > alpha, "above", "below")
  )

mean_scaled_dfreq_McCrone <- df_scatter_McCrone %>%
  group_by(days) %>%
  summarise(mean = mean(abs(delta_q_scaled)))

df_scatter_swine <- swineData %>%
  filter(freq1 > 0) %>%
  mutate(
    days = (DPS2 - DPS1), 
    gen = (DPS2 - DPS1)*(24/hourPerGen),
    q0 = freq1,
    delta_q = freq2 - freq1,
    delta_q_scaled = delta_q / sqrt(q0 * (1 - q0)),
    above_threshold = ifelse(freq2 > alpha, "above", "below")
  )

mean_scaled_dfreq_swine <- df_scatter_swine %>%
  group_by(days) %>%
  summarise(mean = mean(abs(delta_q_scaled)))

df_expected_dq_scaled <- within(data.frame(gen = seq(0, 24)), {
  days <- gen / (24 / hourPerGen)
  exp_q_scaled_McCrone <- sqrt(2 * gen / (pi * MLE_Ne_McCrone))
  exp_q_scaled_swine   <- sqrt(2 * gen / (pi * MLE_Ne_swineData))
})

# 6. Graphs ====
fig_dfreq_scaled_McCrone <- ggplot() +
  geom_point(data = df_scatter_McCrone, 
             aes(x = days, 
                 y = abs(delta_q_scaled), 
                 shape = above_threshold,
                 color = above_threshold), 
             alpha = 0.6, size = 1.5) +
  scale_shape_manual(values = c(19, 1)) +
  scale_color_manual(values = c("#26185F","#0095AF")) +
  geom_line(data = df_expected_dq_scaled, aes(x = days, y = exp_q_scaled_McCrone), color = "red") +
  geom_line(data = mean_scaled_dfreq_McCrone[1:4,], aes(x = days, y = mean), 
            color = "darkgrey", linewidth = 2, alpha = 0.5) +
  geom_point(data = mean_scaled_dfreq_McCrone[1:4,], aes(x = days, y = mean), shape = 4, size = 3) +
  scale_x_continuous(breaks = seq(0, 100, 1)) +
  labs(x = "Days between observations",
       y = expression(abs(Delta*q) / sqrt(q[0] * (1 - q[0]))),
       color = "Second observation relative to threshold",
       shape = "Second observation relative to threshold")

fig_dfreq_scaled_swine <- ggplot() +
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
  geom_line(data = mean_scaled_dfreq_swine, aes(x = days, y = mean), 
            color = "darkgrey", linewidth = 2, alpha = 0.5) +
  geom_point(data = mean_scaled_dfreq_swine, aes(x = days, y = mean), shape = 4, size = 3) +
  labs(x = "Days between observations",
       y = expression(abs(Delta*q) / sqrt(q[0] * (1 - q[0]))),
       color = "Second observation relative to threshold",
       shape = "Second observation relative to threshold")

# 7. Put graphs together ====
fig_dfreq <-fig_dfreq_scaled_McCrone + fig_dfreq_scaled_swine +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A",) &
  theme(legend.position = "bottom", plot.tag = element_text(face = "bold", size = 18))

ggsave(fig_dfreq, file = "fig_dfreq.pdf", path = "02_plots", width = 10, height = 4)
