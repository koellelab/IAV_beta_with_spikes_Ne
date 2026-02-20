# 1. Source functions ====
source("01_Rscripts/WithinhostIAV_Ne_fnuctions.R")

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
mu <- 0 #2 * 10^(-4)
hourPerGen <- 6
m <- 0.004
set.seed(12)

# 4. Input datasets ====
McCrone_IAV_sub2 <- read_csv("03_data/table_S1.csv") %>%
  set_names(c("ENROLLID", "mutation", "DPS1", "DPS2", "freq1", "freq2", "Subset")) %>%
  mutate(gen = (DPS2 - DPS1)*(24/hourPerGen)) %>%
  filter(Subset == 2)

# 5. Generate data frame ====
## 5.1. 0.02 ====
McCrone_IAV_freq0.02 <- McCrone_IAV_sub2
McCrone_IAV_freq0.02$freq1 <- 0.02
McCrone_IAV_freq0.02_gen1 <- McCrone_IAV_freq0.02
McCrone_IAV_freq0.02_gen1$DPS1 <- McCrone_IAV_freq0.02_gen1$DPS2 - 0.25
McCrone_tot_obs_freq0.02 <- adjust_data(McCrone_IAV_freq0.02)

McCrone_tot_obs_freq0.02_gen1 <- adjust_data(McCrone_IAV_freq0.02_gen1)

LLdf_McCrone_freq0.02 <- LL_Ne.df(McCrone_IAV_freq0.02, seq(50, 1000, 5), mu, 0.02, hourPerGen, m)
LLdf_McCrone_freq0.02$prob <- LLdf_McCrone_freq0.02$prob - log(1 - 0.02)
MLE_Ne_McCrone_freq0.02 <- LLdf_McCrone_freq0.02$Ne[which.max(LLdf_McCrone_freq0.02$prob)]
CI_low_Ne_freq0.02 <- LLdf_McCrone_freq0.02 %>%
  filter(Ne < MLE_Ne_McCrone_freq0.02) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_freq0.02 <- LLdf_McCrone_freq0.02 %>%
  filter(Ne > MLE_Ne_McCrone_freq0.02) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

LLdf_McCrone_freq0.02_gen3 <- LL_Ne.df(McCrone_IAV_freq0.02, seq(100, 2000, 5), mu, 0.02, 3, m)
LLdf_McCrone_freq0.02_gen3$prob <- LLdf_McCrone_freq0.02_gen3$prob - log(1 - 0.02)
MLE_Ne_McCrone_freq0.02_gen3 <- LLdf_McCrone_freq0.02_gen3$Ne[which.max(LLdf_McCrone_freq0.02_gen3$prob)]
CI_low_Ne_freq0.02_gen3 <- LLdf_McCrone_freq0.02_gen3 %>%
  filter(Ne < MLE_Ne_McCrone_freq0.02_gen3) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_freq0.02_gen3 <- LLdf_McCrone_freq0.02_gen3 %>%
  filter(Ne > MLE_Ne_McCrone_freq0.02_gen3) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

LLdf_McCrone_freq0.02_gen12 <- LL_Ne.df(McCrone_IAV_freq0.02, seq(100, 2000, 5), mu, 0.02, 12, m)
LLdf_McCrone_freq0.02_gen12$prob <- LLdf_McCrone_freq0.02_gen12$prob - log(1 - 0.02)
MLE_Ne_McCrone_freq0.02_gen12 <- LLdf_McCrone_freq0.02_gen12$Ne[which.max(LLdf_McCrone_freq0.02_gen12$prob)]
CI_low_Ne_freq0.02_gen12 <- LLdf_McCrone_freq0.02_gen12 %>%
  filter(Ne < MLE_Ne_McCrone_freq0.02_gen12) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_freq0.02_gen12 <- LLdf_McCrone_freq0.02_gen12 %>%
  filter(Ne > MLE_Ne_McCrone_freq0.02_gen12) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

LLdf_McCrone_freq0.02_gen1 <- LL_Ne.df(McCrone_IAV_freq0.02_gen1, seq(15, 150, 1), mu, 0.02, 24, m)
LLdf_McCrone_freq0.02_gen1$prob <- LLdf_McCrone_freq0.02_gen1$prob - log(1 - 0.02)
MLE_Ne_McCrone_freq0.02_gen1 <- LLdf_McCrone_freq0.02_gen1$Ne[which.max(LLdf_McCrone_freq0.02_gen1$prob)]
CI_low_Ne_freq0.02_gen1 <- LLdf_McCrone_freq0.02_gen1 %>%
  filter(Ne < MLE_Ne_McCrone_freq0.02_gen1) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_freq0.02_gen1 <- LLdf_McCrone_freq0.02_gen1 %>%
  filter(Ne > MLE_Ne_McCrone_freq0.02_gen1) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

## 5.2. 0.005 ====
McCrone_IAV_freq0.005 <- McCrone_IAV_sub2
McCrone_IAV_freq0.005$freq1 <- 0.005
McCrone_IAV_freq0.005_gen1 <- McCrone_IAV_freq0.005
McCrone_IAV_freq0.005_gen1$DPS2 <- McCrone_IAV_freq0.005_gen1$DPS1 + 1
McCrone_tot_obs_freq0.005 <- adjust_data(McCrone_IAV_freq0.005)

LLdf_McCrone_freq0.005 <- LL_Ne.df(McCrone_IAV_freq0.005, seq(50,400,2), mu, 0.005, hourPerGen, m)
LLdf_McCrone_freq0.005$prob <- LLdf_McCrone_freq0.005$prob - log(1 - 0.005)
MLE_Ne_McCrone_freq0.005 <- LLdf_McCrone_freq0.005$Ne[which.max(LLdf_McCrone_freq0.005$prob)]
CI_low_Ne_freq0.005 <- LLdf_McCrone_freq0.005 %>%
  filter(Ne < MLE_Ne_McCrone_freq0.005) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_freq0.005 <- LLdf_McCrone_freq0.005 %>%
  filter(Ne > MLE_Ne_McCrone_freq0.005) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

LLdf_McCrone_freq0.005_gen3 <- LL_Ne.df(McCrone_IAV_freq0.005, seq(50,800,2), mu, 0.005, 3, m)
LLdf_McCrone_freq0.005_gen3$prob <- LLdf_McCrone_freq0.005_gen3$prob - log(1 - 0.005)
MLE_Ne_McCrone_freq0.005_gen3 <- LLdf_McCrone_freq0.005_gen3$Ne[which.max(LLdf_McCrone_freq0.005_gen3$prob)]
CI_low_Ne_freq0.005_gen3 <- LLdf_McCrone_freq0.005_gen3 %>%
  filter(Ne < MLE_Ne_McCrone_freq0.005_gen3) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_freq0.005_gen3 <- LLdf_McCrone_freq0.005_gen3 %>%
  filter(Ne > MLE_Ne_McCrone_freq0.005_gen3) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

LLdf_McCrone_freq0.005_gen12 <- LL_Ne.df(McCrone_IAV_freq0.005, seq(50,800,2), mu, 0.005, 12, m)
LLdf_McCrone_freq0.005_gen12$prob <- LLdf_McCrone_freq0.005_gen12$prob - log(1 - 0.005)
MLE_Ne_McCrone_freq0.005_gen12 <- LLdf_McCrone_freq0.005_gen12$Ne[which.max(LLdf_McCrone_freq0.005_gen12$prob)]
CI_low_Ne_freq0.005_gen12 <- LLdf_McCrone_freq0.005_gen12 %>%
  filter(Ne < MLE_Ne_McCrone_freq0.005_gen12) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_freq0.005_gen12 <- LLdf_McCrone_freq0.005_gen12 %>%
  filter(Ne > MLE_Ne_McCrone_freq0.005_gen12) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

LLdf_McCrone_freq0.005_gen1 <- LL_Ne.df(McCrone_IAV_freq0.005_gen1, seq(25,40,0.5), mu, 0.005, 24, m)
LLdf_McCrone_freq0.005_gen1$prob <- LLdf_McCrone_freq0.005_gen1$prob - log(1 - 0.005)
MLE_Ne_McCrone_freq0.005_gen1 <- LLdf_McCrone_freq0.005_gen1$Ne[which.max(LLdf_McCrone_freq0.005_gen1$prob)]
CI_low_Ne_freq0.005_gen1 <- LLdf_McCrone_freq0.005_gen1 %>%
  filter(Ne < MLE_Ne_McCrone_freq0.005_gen1) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_freq0.005_gen1 <- LLdf_McCrone_freq0.005_gen1 %>%
  filter(Ne > MLE_Ne_McCrone_freq0.005_gen1) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

# 6. Graphs ====
fig3a <- ggplot(McCrone_tot_obs_freq0.02, aes(x = DPS, y = freq)) +
  geom_point(aes(color = type, shape = type), alpha = 0.5) +
  geom_line(aes(group = ID), alpha = 0.3) +
  geom_hline(yintercept = 0.02, color = "red", linetype = "dashed", ) +
  scale_color_discrete_sequential(palette = "YlGnBu", nmax = 3, rev = FALSE) +
  labs(x = "Days post symptom onset", y = "iSNV frequency", color = NULL, shape = NULL) +
  theme(legend.position = c(0.25, 0.9)) +
  ylim(0, 0.4) +
  xlim(-2, 6)

fig3b <- ggplot(LLdf_McCrone_freq0.02) +
  geom_line(aes(x = Ne, y = prob)) +
  geom_vline(xintercept = MLE_Ne_McCrone_freq0.02, color = "red") +
  geom_hline(yintercept = max(LLdf_McCrone_freq0.02$prob), color = "red", linetype = "dashed") +
  annotate("text", x = 700, y = -55, label = paste("N[E] == ", MLE_Ne_McCrone_freq0.02), parse = TRUE, color = "red", size = 6) +
  annotate("rect", xmin = as.numeric(CI_low_Ne_freq0.02), xmax = as.numeric(CI_high_Ne_freq0.02), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  labs(x = expression(N[E]), y = "Log-likelihood", title = expression(paste("initial frequency ", "= 2%"))) 

fig3c <- ggplot(LLdf_McCrone_freq0.005) +
  geom_line(aes(x = Ne, y = prob)) +
  geom_vline(xintercept = MLE_Ne_McCrone_freq0.005, color = "red") +
  geom_hline(yintercept = max(LLdf_McCrone_freq0.005$prob), color = "red", linetype = "dashed") +
  annotate("text", x = 280, y = -72, label = paste("N[E] == ", MLE_Ne_McCrone_freq0.005), parse = TRUE, color = "red", size = 6) +
  annotate("rect", xmin = as.numeric(CI_low_Ne_freq0.005), xmax = as.numeric(CI_high_Ne_freq0.005), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  labs(x = expression(N[E]), y = "Log-likelihood", title = expression(paste("initial frequency ", "= 0.5%"))) 

fig3d <- ggplot(McCrone_tot_obs_freq0.02_gen1, aes(x = DPS, y = freq)) +
  geom_point(aes(color = type, shape = type), alpha = 0.5) +
  geom_line(aes(group = ID), alpha = 0.3) +
  geom_hline(yintercept = 0.02, color = "red", linetype = "dashed", ) +
  scale_color_discrete_sequential(palette = "YlGnBu", nmax = 3, rev = FALSE) +
  labs(x = "Days post symptom onset", y = "iSNV frequency", color = NULL, shape = NULL) +
  theme(legend.position = c(0.25, 0.9)) +
  ylim(0, 0.4) +
  xlim(-2, 6)

fig3e <- ggplot(LLdf_McCrone_freq0.02_gen1) +
  geom_line(aes(x = Ne, y = prob)) +
  geom_vline(xintercept = MLE_Ne_McCrone_freq0.02_gen1, color = "red") +
  geom_hline(yintercept = max(LLdf_McCrone_freq0.02_gen1$prob), color = "red", linetype = "dashed") +
  annotate("text", x = 100, y = -80, label = paste("N[E] == ", MLE_Ne_McCrone_freq0.02_gen1), parse = TRUE, color = "red", size = 6) +
  annotate("rect", xmin = as.numeric(CI_low_Ne_freq0.02_gen1), xmax = as.numeric(CI_high_Ne_freq0.02_gen1), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  labs(x = expression(N[E]), y = "Log-likelihood", title = expression(paste("initial frequency ", "= 2%"))) 

fig3f <- ggplot(LLdf_McCrone_freq0.005_gen1) +
  geom_line(aes(x = Ne, y = prob)) +
  geom_vline(xintercept = MLE_Ne_McCrone_freq0.005_gen1, color = "red") +
  geom_hline(yintercept = max(LLdf_McCrone_freq0.005_gen1$prob), color = "red", linetype = "dashed") +
  annotate("text", x = 37, y = -370, label = paste("N[E] == ", MLE_Ne_McCrone_freq0.005_gen1), parse = TRUE, color = "red", size = 6) +
  annotate("rect", xmin = as.numeric(CI_low_Ne_freq0.005_gen1), xmax = as.numeric(CI_high_Ne_freq0.005_gen1), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  labs(x = expression(N[E]), y = "Log-likelihood", title = expression(paste("initial frequency ", "= 0.5%")))





fig0.02_gen3 <- ggplot(LLdf_McCrone_freq0.02_gen3) +
  geom_line(aes(x = Ne, y = prob)) +
  geom_vline(xintercept = MLE_Ne_McCrone_freq0.02_gen3, color = "red") +
  geom_hline(yintercept = max(LLdf_McCrone_freq0.02_gen3$prob), color = "red", linetype = "dashed") +
  annotate("text", x = 1200, y = -55, label = paste("N[E] == ", MLE_Ne_McCrone_freq0.02_gen3), parse = TRUE, color = "red", size = 6) +
  annotate("rect", xmin = as.numeric(CI_low_Ne_freq0.02_gen3), xmax = as.numeric(CI_high_Ne_freq0.02_gen3), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  labs(x = expression(N[E]), y = "Log-likelihood", title = "initial frequency = 2%\nviral generation time = 3 hours") 

fig0.02_gen12 <- ggplot(LLdf_McCrone_freq0.02_gen12) +
  geom_line(aes(x = Ne, y = prob)) +
  geom_vline(xintercept = MLE_Ne_McCrone_freq0.02_gen12, color = "red") +
  geom_hline(yintercept = max(LLdf_McCrone_freq0.02_gen12$prob), color = "red", linetype = "dashed") +
  annotate("text", x = 1000, y = -70, label = paste("N[E] == ", MLE_Ne_McCrone_freq0.02_gen12), parse = TRUE, color = "red", size = 6) +
  annotate("rect", xmin = as.numeric(CI_low_Ne_freq0.02_gen12), xmax = as.numeric(CI_high_Ne_freq0.02_gen12), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  labs(x = expression(N[E]), y = "Log-likelihood", title = "initial frequency = 2%\nviral generation time = 12 hours") 

fig0.005_gen3 <- ggplot(LLdf_McCrone_freq0.005_gen3) +
  geom_line(aes(x = Ne, y = prob)) +
  geom_vline(xintercept = MLE_Ne_McCrone_freq0.005_gen3, color = "red") +
  geom_hline(yintercept = max(LLdf_McCrone_freq0.005_gen3$prob), color = "red", linetype = "dashed") +
  annotate("text", x = 600, y = -80, label = paste("N[E] == ", MLE_Ne_McCrone_freq0.005_gen3), parse = TRUE, color = "red", size = 6) +
  annotate("rect", xmin = as.numeric(CI_low_Ne_freq0.005_gen3), xmax = as.numeric(CI_high_Ne_freq0.005_gen3), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  labs(x = expression(N[E]), y = "Log-likelihood", title = "initial frequency = 0.5%\nviral generation time = 3 hours") 

fig0.005_gen12 <- ggplot(LLdf_McCrone_freq0.005_gen12) +
  geom_line(aes(x = Ne, y = prob)) +
  geom_vline(xintercept = MLE_Ne_McCrone_freq0.005_gen12, color = "red") +
  geom_hline(yintercept = max(LLdf_McCrone_freq0.005_gen12$prob), color = "red", linetype = "dashed") +
  annotate("text", x = 400, y = -120, label = paste("N[E] == ", MLE_Ne_McCrone_freq0.005_gen12), parse = TRUE, color = "red", size = 6) +
  annotate("rect", xmin = as.numeric(CI_low_Ne_freq0.005_gen12), xmax = as.numeric(CI_high_Ne_freq0.005_gen12), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  labs(x = expression(N[E]), y = "Log-likelihood", title = "initial frequency = 0.5%\nviral generation time = 12 hours") 

fig0.02_gen6 <- fig3b + ggtitle("initial frequency = 2%\nviral generation time = 6 hours")

fig0.005_gen6 <- fig3c + ggtitle("initial frequency = 0.5%\nviral generation time = 6 hours")

# 7. Put graphs together ====
fig3 <- fig3a + fig3b + fig3c + fig3d + fig3e + fig3f +
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 18))

fig_gen_ini_freq_sens <- fig0.02_gen3 + fig0.005_gen3 + fig0.02_gen6 + fig0.005_gen6 + fig0.02_gen12 + fig0.005_gen12 + 
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 18))

ggsave(fig3, file = "fig_human_ds2pdf", path = "02_plots", width = 15, height = 8)

ggsave(fig_gen_ini_freq_sens, file = "fig_human_ds2_gen_ini_freq_sens.pdf", path = "02_plots", width = 10.5, height = 12)
