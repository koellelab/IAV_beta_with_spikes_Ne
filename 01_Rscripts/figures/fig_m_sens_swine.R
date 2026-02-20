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
                  plot.title = element_text(size = 18),
                  legend.text = element_text(size = 12)))

# 3. Input parameters ====
Ne_values_swine <- seq(3, 30, 0.1)
mu <- 0
alpha <- 0.02
hourPerGen <- 6
m <- 0.004
freq_step <- 0.01

# 4. Input datasets ====
swineData <- read_csv("03_data/table_S2.csv") %>%
  filter(Subset == 1) %>%
  set_names(c("pigID", "segment", "DPS1", "DPS2", "freq1", "freq2", "Subset")) %>%
  mutate(gen = (DPS2 - DPS1)*(24/hourPerGen))

# 5. Generate data frame ====
## 5.2. Figure 2B ====
LLdf_swineData_m01 <- LL_Ne.df(swineData, Ne_values_swine, mu, alpha, hourPerGen, 0.001)
MLE_Ne_swineData_m01 <- LLdf_swineData_m01$Ne[which.max(LLdf_swineData_m01$prob)]
CI_low_Ne_swineData_m01 <- LLdf_swineData_m01 %>%
  filter(Ne < MLE_Ne_swineData_m01) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_swineData_m01 <- LLdf_swineData_m01 %>%
  filter(Ne > MLE_Ne_swineData_m01) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

LLdf_swineData_m04 <- LL_Ne.df(swineData, Ne_values_swine, mu, alpha, hourPerGen, 0.004)
MLE_Ne_swineData_m04 <- LLdf_swineData_m04$Ne[which.max(LLdf_swineData_m04$prob)]
CI_low_Ne_swineData_m04 <- LLdf_swineData_m04 %>%
  filter(Ne < MLE_Ne_swineData_m04) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_swineData_m04 <- LLdf_swineData_m04 %>%
  filter(Ne > MLE_Ne_swineData_m04) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

LLdf_swineData_m1 <- LL_Ne.df(swineData, Ne_values_swine, mu, alpha, hourPerGen, 0.01)
MLE_Ne_swineData_m1 <- LLdf_swineData_m1$Ne[which.max(LLdf_swineData_m1$prob)]
CI_low_Ne_swineData_m1 <- LLdf_swineData_m1 %>%
  filter(Ne < MLE_Ne_swineData_m1) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_swineData_m1 <- LLdf_swineData_m1 %>%
  filter(Ne > MLE_Ne_swineData_m1) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

LLdf_swineData_m2 <- LL_Ne.df(swineData, Ne_values_swine, mu, alpha, hourPerGen, 0.02)
MLE_Ne_swineData_m2 <- LLdf_swineData_m2$Ne[which.max(LLdf_swineData_m2$prob)]
CI_low_Ne_swineData_m2 <- LLdf_swineData_m2 %>%
  filter(Ne < MLE_Ne_swineData_m2) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_swineData_m2 <- LLdf_swineData_m2 %>%
  filter(Ne > MLE_Ne_swineData_m2) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

# 6. Graphs ====
fig_m01_swine <- ggplot(LLdf_swineData_m01) +
  geom_line(aes(x = Ne, y = prob)) +
  geom_vline(xintercept = MLE_Ne_swineData, color = "red") +
  geom_hline(yintercept = max(LLdf_swineData_m01$prob), color = "red", linetype = "dashed") +
  annotate("text", x = 25, y = -175, label = paste("N[E] == ", MLE_Ne_swineData), parse = TRUE, color = "red") +
  annotate("rect", xmin = as.numeric(CI_low_Ne_swineData), xmax = as.numeric(CI_high_Ne_swineData), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  labs(x = expression(N[E]), y = "Log-likelihood", title = "m = 0.001") +
  scale_x_continuous(breaks = seq(0, 50, 5), limits = c(0, 31))

fig_m04_swine <- ggplot(LLdf_swineData_m04) +
  geom_line(aes(x = Ne, y = prob)) +
  geom_vline(xintercept = MLE_Ne_swineData, color = "red") +
  geom_hline(yintercept = max(LLdf_swineData_m04$prob), color = "red", linetype = "dashed") +
  annotate("text", x = 25, y = -175, label = paste("N[E] == ", MLE_Ne_swineData), parse = TRUE, color = "red") +
  annotate("rect", xmin = as.numeric(CI_low_Ne_swineData), xmax = as.numeric(CI_high_Ne_swineData), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  labs(x = expression(N[E]), y = "Log-likelihood", title = "m = 0.004 (applied)") +
  scale_x_continuous(breaks = seq(0, 50, 5), limits = c(0, 31))

fig_m1_swine <- ggplot(LLdf_swineData_m1) +
  geom_line(aes(x = Ne, y = prob)) +
  geom_vline(xintercept = MLE_Ne_swineData, color = "red") +
  geom_hline(yintercept = max(LLdf_swineData_m1$prob), color = "red", linetype = "dashed") +
  annotate("text", x = 25, y = -175, label = paste("N[E] == ", MLE_Ne_swineData), parse = TRUE, color = "red") +
  annotate("rect", xmin = as.numeric(CI_low_Ne_swineData), xmax = as.numeric(CI_high_Ne_swineData), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  labs(x = expression(N[E]), y = "Log-likelihood", title = "m = 0.01") +
  scale_x_continuous(breaks = seq(0, 50, 5), limits = c(0, 31))

fig_m2_swine <- ggplot(LLdf_swineData_m2) +
  geom_line(aes(x = Ne, y = prob)) +
  geom_vline(xintercept = MLE_Ne_swineData, color = "red") +
  geom_hline(yintercept = max(LLdf_swineData_m2$prob), color = "red", linetype = "dashed") +
  annotate("text", x = 25, y = -175, label = paste("N[E] == ", MLE_Ne_swineData), parse = TRUE, color = "red") +
  annotate("rect", xmin = as.numeric(CI_low_Ne_swineData), xmax = as.numeric(CI_high_Ne_swineData), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  labs(x = expression(N[E]), y = "Log-likelihood", title = "m = 0.02") +
  scale_x_continuous(breaks = seq(0, 50, 5), limits = c(0, 31))

# 7. Put graphs together ====
fig_m_sens_swine <- fig_m01_swine + fig_m04_swine + fig_m1_swine + fig_m2_swine +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 18))

ggsave(fig_m_sens_swine, file = "fig_m_sens_swine.pdf", path = "02_plots", width = 10, height = 8)
