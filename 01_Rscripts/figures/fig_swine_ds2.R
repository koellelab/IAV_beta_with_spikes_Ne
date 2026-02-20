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
                  axis.text = element_text(size = 15),
                  axis.title = element_text(size = 18),
                  legend.text = element_text(size = 15),
                  plot.title = element_text(size = 20)))

# 3. Input parameters ====
mu <- 0
hourPerGen <- 6
m <- 0.004
set.seed(12)

# 4. Input datasets ====
swineData <- read_csv("03_data/table_S2.csv") %>%
  set_names(c("ENROLLID", "mutation", "DPS1", "DPS2", "freq1", "freq2", "Subset")) %>%
  mutate(gen = (DPS2 - DPS1)*(24/hourPerGen)) %>%
  filter(Subset == 2)

# 5. Generate data frame ====
## 5.1. 0.02 ====
swine_freq0.02 <- swineData
swine_freq0.02$freq1 <- 0.02
swine_freq0.02_gen1 <- swine_freq0.02
swine_freq0.02_gen1$DPS1 <- swine_freq0.02_gen1$DPS2 - 0.25
swine_tot_obs_freq0.02 <- adjust_data(swine_freq0.02)

swine_tot_obs_freq0.02_gen1 <- adjust_data(swine_freq0.02_gen1)

LLdf_swine_freq0.02 <- LL_Ne.df(swine_freq0.02, seq(5, 150, 1), mu, 0.02, hourPerGen, m)
LLdf_swine_freq0.02$prob <- LLdf_swine_freq0.02$prob - log(1 - 0.02)
MLE_Ne_swine_freq0.02 <- LLdf_swine_freq0.02$Ne[which.max(LLdf_swine_freq0.02$prob)]
CI_low_Ne_swine_freq0.02 <- LLdf_swine_freq0.02 %>%
  filter(Ne < MLE_Ne_swine_freq0.02) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_swine_freq0.02 <- LLdf_swine_freq0.02 %>%
  filter(Ne > MLE_Ne_swine_freq0.02) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

LLdf_swine_freq0.02_gen3 <- LL_Ne.df(swine_freq0.02, seq(5, 150, 2), mu, 0.02, 3, m)
LLdf_swine_freq0.02_gen3$prob <- LLdf_swine_freq0.02_gen3$prob - log(1 - 0.02)
MLE_Ne_swine_freq0.02_gen3 <- LLdf_swine_freq0.02_gen3$Ne[which.max(LLdf_swine_freq0.02_gen3$prob)]
CI_low_Ne_swine_freq0.02_gen3 <- LLdf_swine_freq0.02_gen3 %>%
  filter(Ne < MLE_Ne_swine_freq0.02_gen3) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_swine_freq0.02_gen3 <- LLdf_swine_freq0.02_gen3 %>%
  filter(Ne > MLE_Ne_swine_freq0.02_gen3) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

LLdf_swine_freq0.02_gen12 <- LL_Ne.df(swine_freq0.02, seq(5, 150, 2), mu, 0.02, 12, m)
LLdf_swine_freq0.02_gen12$prob <- LLdf_swine_freq0.02_gen12$prob - log(1 - 0.02)
MLE_Ne_swine_freq0.02_gen12 <- LLdf_swine_freq0.02_gen12$Ne[which.max(LLdf_swine_freq0.02_gen12$prob)]
CI_low_Ne_swine_freq0.02_gen12 <- LLdf_swine_freq0.02_gen12 %>%
  filter(Ne < MLE_Ne_swine_freq0.02_gen12) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_swine_freq0.02_gen12 <- LLdf_swine_freq0.02_gen12 %>%
  filter(Ne > MLE_Ne_swine_freq0.02_gen12) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

LLdf_swine_freq0.02_gen1 <- LL_Ne.df(swine_freq0.02_gen1, seq(10, 40, 0.5), mu, 0.02, hourPerGen, m)
LLdf_swine_freq0.02_gen1$prob <- LLdf_swine_freq0.02_gen1$prob - log(1 - 0.02)
MLE_Ne_swine_freq0.02_gen1 <- LLdf_swine_freq0.02_gen1$Ne[which.max(LLdf_swine_freq0.02_gen1$prob)]
CI_low_Ne_swine_freq0.02_gen1 <- LLdf_swine_freq0.02_gen1 %>%
  filter(Ne < MLE_Ne_swine_freq0.02_gen1) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_swine_freq0.02_gen1 <- LLdf_swine_freq0.02_gen1 %>%
  filter(Ne > MLE_Ne_swine_freq0.02_gen1) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

## 5.2. 0.005 ====
swine_freq0.005 <- swineData
swine_freq0.005$freq1 <- 0.005
swine_freq0.005_gen1 <- swine_freq0.005
swine_freq0.005_gen1$DPS1 <- swine_freq0.005_gen1$DPS2 - 0.25
swine_tot_obs_freq0.005 <- adjust_data(swine_freq0.005)

LLdf_swine_freq0.005 <- LL_Ne.df(swine_freq0.005, seq(10, 100, 1), mu, 0.005, hourPerGen, m)
LLdf_swine_freq0.005$prob <- LLdf_swine_freq0.005$prob - log(1 - 0.005)
MLE_Ne_swine_freq0.005 <- LLdf_swine_freq0.005$Ne[which.max(LLdf_swine_freq0.005$prob)]
CI_low_Ne_swine_freq0.005 <- LLdf_swine_freq0.005 %>%
  filter(Ne < MLE_Ne_swine_freq0.005) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_swine_freq0.005 <- LLdf_swine_freq0.005 %>%
  filter(Ne > MLE_Ne_swine_freq0.005) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

LLdf_swine_freq0.005_gen3 <- LL_Ne.df(swine_freq0.005, seq(10, 100, 1), mu, 0.005, 3, m)
LLdf_swine_freq0.005_gen3$prob <- LLdf_swine_freq0.005_gen3$prob - log(1 - 0.005)
MLE_Ne_swine_freq0.005_gen3 <- LLdf_swine_freq0.005_gen3$Ne[which.max(LLdf_swine_freq0.005_gen3$prob)]
CI_low_Ne_swine_freq0.005_gen3 <- LLdf_swine_freq0.005_gen3 %>%
  filter(Ne < MLE_Ne_swine_freq0.005_gen3) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_swine_freq0.005_gen3 <- LLdf_swine_freq0.005_gen3 %>%
  filter(Ne > MLE_Ne_swine_freq0.005_gen3) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

LLdf_swine_freq0.005_gen12 <- LL_Ne.df(swine_freq0.005, seq(9, 100, 1), mu, 0.005, 12, m)
LLdf_swine_freq0.005_gen12$prob <- LLdf_swine_freq0.005_gen12$prob - log(1 - 0.005)
MLE_Ne_swine_freq0.005_gen12 <- LLdf_swine_freq0.005_gen12$Ne[which.max(LLdf_swine_freq0.005_gen12$prob)]
CI_low_Ne_swine_freq0.005_gen12 <- LLdf_swine_freq0.005_gen12 %>%
  filter(Ne < MLE_Ne_swine_freq0.005_gen12) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_swine_freq0.005_gen12 <- LLdf_swine_freq0.005_gen12 %>%
  filter(Ne > MLE_Ne_swine_freq0.005_gen12) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

LLdf_swine_freq0.005_gen1 <- LL_Ne.df(swine_freq0.005_gen1, seq(10,30,0.5), mu, 0.005, hourPerGen, m)
LLdf_swine_freq0.005_gen1$prob <- LLdf_swine_freq0.005_gen1$prob - log(1 - 0.005)
MLE_Ne_swine_freq0.005_gen1 <- LLdf_swine_freq0.005_gen1$Ne[which.max(LLdf_swine_freq0.005_gen1$prob)]
CI_low_Ne_swine_freq0.005_gen1 <- LLdf_swine_freq0.005_gen1 %>%
  filter(Ne < MLE_Ne_swine_freq0.005_gen1) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_swine_freq0.005_gen1 <- LLdf_swine_freq0.005_gen1 %>%
  filter(Ne > MLE_Ne_swine_freq0.005_gen1) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

# 6. Graphs ====
fig5a <- ggplot(swine_tot_obs_freq0.02, aes(x = DPS, y = freq)) +
  geom_point(aes(color = type, shape = type), alpha = 0.5) +
  geom_line(aes(group = ID), alpha = 0.3) +
  geom_hline(yintercept = 0.02, color = "red", linetype = "dashed", ) +
  scale_color_discrete_sequential(palette = "YlGnBu", nmax = 3, rev = FALSE) +
  labs(x = "Day of fair", y = "iSNV frequency", color = NULL, shape = NULL) +
  theme(legend.position = c(0.25, 0.9)) +
  xlim(0, 6)

fig5b <- ggplot(LLdf_swine_freq0.02) +
  geom_line(aes(x = Ne, y = prob)) +
  geom_vline(xintercept = MLE_Ne_swine_freq0.02, color = "red") +
  geom_hline(yintercept = max(LLdf_swine_freq0.02$prob), color = "red", linetype = "dashed") +
  annotate("text", x = 75, y = -320, label = paste("N[E] == ", MLE_Ne_swine_freq0.02), parse = TRUE, color = "red", size = 6) +
  annotate("rect", xmin = as.numeric(CI_low_Ne_swine_freq0.02), xmax = as.numeric(CI_high_Ne_swine_freq0.02), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  labs(x = expression(N[E]), y = "Log-likelihood", title = expression(paste("initial frequency ", "= 2%"))) 

fig5c <- ggplot(LLdf_swine_freq0.005) +
  geom_line(aes(x = Ne, y = prob)) +
  geom_vline(xintercept = MLE_Ne_swine_freq0.005, color = "red") +
  geom_hline(yintercept = max(LLdf_swine_freq0.005$prob), color = "red", linetype = "dashed") +
  annotate("text", x = 60, y = -375, label = paste("N[E] == ", MLE_Ne_swine_freq0.005), parse = TRUE, color = "red", size = 6) +
  annotate("rect", xmin = as.numeric(CI_low_Ne_swine_freq0.005), xmax = as.numeric(CI_high_Ne_swine_freq0.005), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  labs(x = expression(N[E]), y = "Log-likelihood", title = expression(paste("initial frequency ", "= 0.5%"))) 

fig5d <- ggplot(swine_tot_obs_freq0.02_gen1, aes(x = DPS, y = freq)) +
  geom_point(aes(color = type, shape = type), alpha = 0.5) +
  geom_line(aes(group = ID), alpha = 0.3) +
  geom_hline(yintercept = 0.02, color = "red", linetype = "dashed", ) +
  scale_color_discrete_sequential(palette = "YlGnBu", nmax = 3, rev = FALSE) +
  labs(x = "Day of fair", y = "iSNV frequency", color = NULL, shape = NULL) +
  theme(legend.position = c(0.25, 0.9)) +
  xlim(0, 6)

fig5e <- ggplot(LLdf_swine_freq0.02_gen1) +
  geom_line(aes(x = Ne, y = prob)) +
  geom_vline(xintercept = MLE_Ne_swine_freq0.02_gen1, color = "red") +
  geom_hline(yintercept = max(LLdf_swine_freq0.02_gen1$prob), color = "red", linetype = "dashed") +
  annotate("text", x = 32, y = -650, label = paste("N[E] == ", MLE_Ne_swine_freq0.02_gen1), parse = TRUE, color = "red", size = 6) +
  annotate("rect", xmin = as.numeric(CI_low_Ne_swine_freq0.02_gen1), xmax = as.numeric(CI_high_Ne_swine_freq0.02_gen1), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  labs(x = expression(N[E]), y = "Log-likelihood", title = expression(paste("initial frequency ", "= 2%"))) 

fig5f <- ggplot(LLdf_swine_freq0.005_gen1) +
  geom_line(aes(x = Ne, y = prob)) +
  geom_vline(xintercept = MLE_Ne_swine_freq0.005_gen1, color = "red") +
  geom_hline(yintercept = max(LLdf_swine_freq0.005_gen1$prob), color = "red", linetype = "dashed") +
  annotate("text", x = 15, y = -9000, label = paste("N[E] == ", MLE_Ne_swine_freq0.005_gen1), parse = TRUE, color = "red", size = 6) +
  annotate("rect", xmin = as.numeric(CI_low_Ne_swine_freq0.005_gen1), xmax = as.numeric(CI_high_Ne_swine_freq0.005_gen1), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  labs(x = expression(N[E]), y = "Log-likelihood", title = expression(paste("initial frequency ", "= 0.5%")))




fig0.02_gen3_swine <- ggplot(LLdf_swine_freq0.02_gen3) +
  geom_line(aes(x = Ne, y = prob)) +
  geom_vline(xintercept = MLE_Ne_swine_freq0.02_gen3, color = "red") +
  geom_hline(yintercept = max(LLdf_swine_freq0.02_gen3$prob), color = "red", linetype = "dashed") +
  annotate("text", x = 120, y = -350, label = paste("N[E] == ", MLE_Ne_swine_freq0.02_gen3), parse = TRUE, color = "red", size = 6) +
  annotate("rect", xmin = as.numeric(CI_low_Ne_swine_freq0.02_gen3), xmax = as.numeric(CI_high_Ne_swine_freq0.02_gen3), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  labs(x = expression(N[E]), y = "Log-likelihood", title = "initial frequency = 2%\nviral generation time = 3 hours") 

fig0.02_gen12_swine <- ggplot(LLdf_swine_freq0.02_gen12) +
  geom_line(aes(x = Ne, y = prob)) +
  geom_vline(xintercept = MLE_Ne_swine_freq0.02_gen12, color = "red") +
  geom_hline(yintercept = max(LLdf_swine_freq0.02_gen12$prob), color = "red", linetype = "dashed") +
  annotate("text", x = 50, y = -425, label = paste("N[E] == ", MLE_Ne_swine_freq0.02_gen12), parse = TRUE, color = "red", size = 6) +
  annotate("rect", xmin = as.numeric(CI_low_Ne_swine_freq0.02_gen12), xmax = as.numeric(CI_high_Ne_swine_freq0.02_gen12), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  labs(x = expression(N[E]), y = "Log-likelihood", title = "initial frequency = 2%\nviral generation time = 12 hours") 

fig0.005_gen3_swine <- ggplot(LLdf_swine_freq0.005_gen3) +
  geom_line(aes(x = Ne, y = prob)) +
  geom_vline(xintercept = MLE_Ne_swine_freq0.005_gen3, color = "red") +
  geom_hline(yintercept = max(LLdf_swine_freq0.005_gen3$prob), color = "red", linetype = "dashed") +
  annotate("text", x = 50, y = -425, label = paste("N[E] == ", MLE_Ne_swine_freq0.005_gen3), parse = TRUE, color = "red", size = 6) +
  annotate("rect", xmin = as.numeric(CI_low_Ne_swine_freq0.005_gen3), xmax = as.numeric(CI_high_Ne_swine_freq0.005_gen3), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  labs(x = expression(N[E]), y = "Log-likelihood", title = "initial frequency = 0.5%\nviral generation time = 3 hours") 

fig0.005_gen12_swine <- ggplot(LLdf_swine_freq0.005_gen12) +
  geom_line(aes(x = Ne, y = prob)) +
  geom_vline(xintercept = MLE_Ne_swine_freq0.005_gen12, color = "red") +
  geom_hline(yintercept = max(LLdf_swine_freq0.005_gen12$prob), color = "red", linetype = "dashed") +
  annotate("text", x = 50, y = -1250, label = paste("N[E] == ", MLE_Ne_swine_freq0.005_gen12), parse = TRUE, color = "red", size = 6) +
  annotate("rect", xmin = as.numeric(CI_low_Ne_swine_freq0.005_gen12), xmax = as.numeric(CI_high_Ne_swine_freq0.005_gen12), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  labs(x = expression(N[E]), y = "Log-likelihood", title = "initial frequency = 0.5%\nviral generation time = 12 hours") 

fig0.02_gen6_swine <- fig5b + ggtitle("initial frequency = 2%\nviral generation time = 6 hours")

fig0.005_gen6_swine <- fig5c + ggtitle("initial frequency = 0.5%\nviral generation time = 6 hours")

# 7. Put graphs together ====
fig5 <- fig5a + fig5b + fig5c + fig5d + fig5e + fig5f +
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 18))

fig_gen_ini_freq_sens_swine <- fig0.02_gen3_swine + fig0.005_gen3_swine + 
  fig0.02_gen6_swine + fig0.005_gen6_swine + 
  fig0.02_gen12_swine + fig0.005_gen12_swine + 
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 18))

ggsave(fig5, file = "fig5.pdf", path = "02_plots", width = 15, height = 8)

ggsave(fig_gen_ini_freq_sens_swine, file = "swine_ds2_gen_ini_freq_sens.pdf", path = "02_plots", width = 10.5, height = 12)
