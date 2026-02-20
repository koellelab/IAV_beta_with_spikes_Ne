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
                  plot.title = element_text(size = 18),
                  legend.text = element_text(size = 12)))

# 3. Input parameters ====
Ne_values_human <- seq(10, 100, 1)
Ne_values_swine <- seq(3, 30, 0.1)
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

swineData <- read_csv("03_data/table_S2.csv") %>%
  filter(Subset == 1) %>%
  set_names(c("pigID", "segment", "DPS1", "DPS2", "freq1", "freq2", "Subset")) %>%
  mutate(gen = (DPS2 - DPS1)*(24/hourPerGen))

# 5. Generate data frame ====
## 5.1. Initial distribution ====
df_bs_gen0 <- data.frame(freq = seq(0, 1, 0.01),
                         generation = 0,
                         prob_m01 = dnorm(df_bs_gen0_ne20$freq, mean = ini_p, sd = sqrt(0.001*ini_p*(1-ini_p))),
                         prob_m04 = dnorm(df_bs_gen0_ne20$freq, mean = ini_p, sd = sqrt(0.004*ini_p*(1-ini_p))),
                         prob_m1 = dnorm(df_bs_gen0_ne20$freq, mean = ini_p, sd = sqrt(0.01*ini_p*(1-ini_p))),
                         prob_m2 = dnorm(df_bs_gen0_ne20$freq, mean = ini_p, sd = sqrt(0.02*ini_p*(1-ini_p))))

df_bs_gen0 <- df_bs_gen0 %>%
  mutate(prob_m01 = prob_m01/sum(prob_m01),
         prob_m04 = prob_m04/sum(prob_m04),
         prob_m1 = prob_m1/sum(prob_m1),
         prob_m2 = prob_m2/sum(prob_m2))

## 5.2. Human ====
LLdf_McCrone_m01 <- LL_Ne.df(McCrone_IAV, Ne_values_human, mu, alpha, hourPerGen, 0.001)
MLE_Ne_McCrone_m01 <- LLdf_McCrone_m01$Ne[which.max(LLdf_McCrone_m01$prob)]
CI_low_Ne_McCrone_m01 <- LLdf_McCrone_m01 %>%
  filter(Ne < MLE_Ne_McCrone_m01) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_McCrone_m01 <- LLdf_McCrone_m01 %>%
  filter(Ne > MLE_Ne_McCrone_m01) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

LLdf_McCrone_m04 <- LL_Ne.df(McCrone_IAV, Ne_values_human, mu, alpha, hourPerGen, 0.004)
MLE_Ne_McCrone_m04 <- LLdf_McCrone_m04$Ne[which.max(LLdf_McCrone_m04$prob)]
CI_low_Ne_McCrone_m04 <- LLdf_McCrone_m04 %>%
  filter(Ne < MLE_Ne_McCrone_m04) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_McCrone_m04 <- LLdf_McCrone_m04 %>%
  filter(Ne > MLE_Ne_McCrone_m04) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

LLdf_McCrone_m1 <- LL_Ne.df(McCrone_IAV, Ne_values_human, mu, alpha, hourPerGen, 0.01)
MLE_Ne_McCrone_m1 <- LLdf_McCrone_m1$Ne[which.max(LLdf_McCrone_m1$prob)]
CI_low_Ne_McCrone_m1 <- LLdf_McCrone_m1 %>%
  filter(Ne < MLE_Ne_McCrone_m1) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_McCrone_m1 <- LLdf_McCrone_m1 %>%
  filter(Ne > MLE_Ne_McCrone_m1) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

LLdf_McCrone_m2 <- LL_Ne.df(McCrone_IAV, Ne_values_human, mu, alpha, hourPerGen, 0.02)
MLE_Ne_McCrone_m2 <- LLdf_McCrone_m2$Ne[which.max(LLdf_McCrone_m2$prob)]
CI_low_Ne_McCrone_m2 <- LLdf_McCrone_m2 %>%
  filter(Ne < MLE_Ne_McCrone_m2) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_McCrone_m2 <- LLdf_McCrone_m2 %>%
  filter(Ne > MLE_Ne_McCrone_m2) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

## 5.3. Swine ====
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
## 6.1. Initial DAF ====
fig_m01_ini <- ggplot(df_bs_gen0) +
  geom_line(aes(x = freq, y = prob_m01)) +
  geom_point(aes(x = freq, y = prob_m01)) +
  labs(y = "Probability", x = "Allele frequency", title = "m = 0.001") 

fig_m04_ini <- ggplot(df_bs_gen0) +
  geom_line(aes(x = freq, y = prob_m04)) +
  geom_point(aes(x = freq, y = prob_m04)) +
  labs(y = "Probability", x = "Allele frequency", title = "m = 0.004 (applied)")

fig_m1_ini <- ggplot(df_bs_gen0) +
  geom_line(aes(x = freq, y = prob_m1)) +
  geom_point(aes(x = freq, y = prob_m1)) +
  labs(y = "Probability", x = "Allele frequency", title = "m = 0.01")

fig_m2_ini <- ggplot(df_bs_gen0) +
  geom_line(aes(x = freq, y = prob_m2)) +
  geom_point(aes(x = freq, y = prob_m2)) +
  labs(y = "Probability", x = "Allele frequency", title = "m = 0.02")

## 6.2. Human ====
fig_m01_human <- ggplot(LLdf_McCrone_m01) +
  geom_line(aes(x = Ne, y = prob)) +
  geom_vline(xintercept = MLE_Ne_McCrone_m01, color = "red") +
  geom_hline(yintercept = max(LLdf_McCrone_m01$prob), color = "red", linetype = "dashed") +
  annotate("text", x = 70, y = -40, label = paste("N[E] == ", MLE_Ne_McCrone_m01), parse = TRUE, color = "red") +
  annotate("rect", xmin = as.numeric(CI_low_Ne_McCrone_m01), xmax = as.numeric(CI_high_Ne_McCrone_m01), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  labs(x = expression(N[E]), y = "Log-likelihood", title = "m = 0.001") +
  scale_x_continuous(breaks = seq(0, 200, 20), limits = c(0, 100))

fig_m04_human <- ggplot(LLdf_McCrone_m04) +
  geom_line(aes(x = Ne, y = prob)) +
  geom_vline(xintercept = MLE_Ne_McCrone_m04, color = "red") +
  geom_hline(yintercept = max(LLdf_McCrone_m04$prob), color = "red", linetype = "dashed") +
  annotate("text", x = 70, y = -40, label = paste("N[E] == ", MLE_Ne_McCrone_m04), parse = TRUE, color = "red") +
  annotate("rect", xmin = as.numeric(CI_low_Ne_McCrone_m04), xmax = as.numeric(CI_high_Ne_McCrone_m04), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  labs(x = expression(N[E]), y = "Log-likelihood", title = "m = 0.004 (applied)") +
  scale_x_continuous(breaks = seq(0, 200, 20), limits = c(0, 100))

fig_m1_human <- ggplot(LLdf_McCrone_m1) +
  geom_line(aes(x = Ne, y = prob)) +
  geom_vline(xintercept = MLE_Ne_McCrone_m1, color = "red") +
  geom_hline(yintercept = max(LLdf_McCrone_m1$prob), color = "red", linetype = "dashed") +
  annotate("text", x = 70, y = -40, label = paste("N[E] == ", MLE_Ne_McCrone_m1), parse = TRUE, color = "red") +
  annotate("rect", xmin = as.numeric(CI_low_Ne_McCrone_m1), xmax = as.numeric(CI_high_Ne_McCrone_m1), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  labs(x = expression(N[E]), y = "Log-likelihood", title = "m = 0.01") +
  scale_x_continuous(breaks = seq(0, 200, 20), limits = c(0, 100))

fig_m2_human <- ggplot(LLdf_McCrone_m2) +
  geom_line(aes(x = Ne, y = prob)) +
  geom_vline(xintercept = MLE_Ne_McCrone_m2, color = "red") +
  geom_hline(yintercept = max(LLdf_McCrone_m2$prob), color = "red", linetype = "dashed") +
  annotate("text", x = 70, y = -40, label = paste("N[E] == ", MLE_Ne_McCrone_m2), parse = TRUE, color = "red") +
  annotate("rect", xmin = as.numeric(CI_low_Ne_McCrone_m2), xmax = as.numeric(CI_high_Ne_McCrone_m2), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  labs(x = expression(N[E]), y = "Log-likelihood", title = "m = 0.02") +
  scale_x_continuous(breaks = seq(0, 200, 20), limits = c(0, 100))

## 6.3. Swine ====
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
row_ini <- fig_m01_ini + fig_m04_ini + fig_m1_ini + fig_m2_ini +
  plot_layout(ncol = 4) +
  plot_annotation(title = "Initial DAF", tag_levels = "A",
                  theme = theme(plot.title = element_text(face = "bold", size = 18))) &
  theme(plot.tag = element_text(face = "bold", size = 18))

row_human <- fig_m01_human + fig_m04_human + fig_m1_human + fig_m2_human +
  plot_layout(ncol = 4) +
  plot_annotation(title = "Human data subset 1", tag_levels = list(c("E", "F", "G", "H")),
                  theme = theme(plot.title = element_text(face = "bold", size = 18))) &
  theme(plot.tag = element_text(face = "bold", size = 18))

row_swine <- fig_m01_swine + fig_m04_swine + fig_m1_swine + fig_m2_swine +
  plot_layout(ncol = 4) +
  plot_annotation(title = "Swine data subset 1", tag_levels = list(c("I", "J", "K", "L")),
                  theme = theme(plot.title = element_text(face = "bold", size = 18))) &
  theme(plot.tag = element_text(face = "bold", size = 18))

fig_m_sens <- wrap_elements(row_ini) / wrap_elements(row_human) / wrap_elements(row_swine) 

ggsave(fig_m_sens, file = "fig_m_sens.pdf", path = "02_plots", width = 20, height = 15)
