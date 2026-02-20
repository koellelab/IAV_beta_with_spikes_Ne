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
Ne_values <- seq(20, 200, 1)
mu <- 0
alpha <- 0.005
hourPerGen <- 6
m <- 0.004
freq_step <- 0.01

# 4. Input datasets ====
McCrone_IAV_0.005 <- read_csv("03_data/table_S1_05.csv") %>%
  set_names(c("ENROLLID", "mutation", "DPS1", "DPS2", "freq1", "freq2", "Subset")) %>%
  mutate(gen = (DPS2 - DPS1)*(24/hourPerGen)) %>%
  filter(Subset == 1)

# 5. Generate data frame ====
## 5.1. Figure 2A ====
# data frame only including first round of observations
McCrone_obs1_0.005 <- McCrone_IAV_0.005 %>%
  select(mutation, ENROLLID, DPS1, freq1, gen)
McCrone_obs1_0.005$type <- c("obs1")
names(McCrone_obs1_0.005) <- c("mutation", "ENROLLID", "DPS", "freq", "gen", "type")
# data frame only including second round of observations
McCrone_obs2_0.005 <- McCrone_IAV_0.005 %>%
  select(mutation, ENROLLID, DPS2, freq2, gen)
McCrone_obs2_0.005$type <- c("obs2")
names(McCrone_obs2_0.005) <- c("mutation", "ENROLLID", "DPS", "freq", "gen", "type")
# merge the data together
McCrone_tot_obs_0.005 <- rbind(McCrone_obs1_0.005, McCrone_obs2_0.005)
McCrone_tot_obs_0.005[McCrone_tot_obs_0.005$freq < alpha, "type"] <- "observation 2 below 0.5%"
McCrone_tot_obs_0.005[McCrone_tot_obs_0.005$freq < alpha, "freq"] <- alpha
McCrone_tot_obs_0.005$type <- factor(McCrone_tot_obs_0.005$type, labels = c("observation 1", "observation 2", "observation 2 below 0.5%"))

## 5.2. Figure 2B ====
LLdf_McCrone_0.005 <- LL_Ne.df(McCrone_IAV_0.005, Ne_values, mu, alpha, hourPerGen, m)
MLE_Ne_McCrone_0.005 <- LLdf_McCrone_0.005$Ne[which.max(LLdf_McCrone_0.005$prob)]
CI_low_Ne_0.005 <- LLdf_McCrone_0.005 %>%
  filter(Ne < MLE_Ne_McCrone_0.005) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_0.005 <- LLdf_McCrone_0.005 %>%
  filter(Ne > MLE_Ne_McCrone_0.005) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

# 6. Graphs ====
fig2a_0.005 <- ggplot(McCrone_tot_obs_0.005, aes(x = DPS, y = freq)) +
  geom_hline(yintercept = alpha, color = "red", linetype = "dashed", ) +
  geom_point(aes(color = type, shape = type), alpha = 0.5) +
  geom_line(aes(group = ENROLLID), alpha = 0.3) +
  scale_shape_manual(values = c(19, 17, 2))+
  scale_color_manual(values = c("#26185F", "#0095AF", "black"))+
  labs(x = "Days post symptom onset", y = "iSNV frequency", color = NULL, shape = NULL) +
  theme(legend.position = c(0.3, 0.93)) +
  ylim(0, 0.7)

fig2b_0.005 <- ggplot(LLdf_McCrone_0.005) +
  geom_line(aes(x = Ne, y = prob)) +
  geom_vline(xintercept = MLE_Ne_McCrone_0.005, color = "red") +
  geom_hline(yintercept = max(LLdf_McCrone_0.005$prob), color = "red", linetype = "dashed") +
  annotate("text", x = 60, y = -95, label = paste("N[E] == ", MLE_Ne_McCrone_0.005), parse = TRUE, color = "red") +
  annotate("rect", xmin = as.numeric(CI_low_Ne_0.005), xmax = as.numeric(CI_high_Ne_0.005), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  labs(x = expression(N[E]), y = "Log-likelihood") +
  scale_x_continuous(breaks = seq(0, 200, 20))

# 7. Put graphs together ====
fig2_0.005 <- fig2a_0.005 + fig2b_0.005 + 
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 18))
