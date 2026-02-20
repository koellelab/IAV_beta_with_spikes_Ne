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
                  legend.title = element_text(size = 13),
                  legend.text = element_text(size = 12)))

# 3. Input parameters ====
Ne_values <- seq(5, 200, 1)
mu <- 0
alpha <- 0.02
m <- 0.004
freq_step <- 0.01

# 4. Input datasets ====
McCrone_IAV <- read_csv("03_data/table_S1.csv") %>%
  set_names(c("ENROLLID", "mutation", "DPS1", "DPS2", "freq1", "freq2", "Subset")) %>%
  filter(Subset == 1)

swineData <- read_csv("03_data/table_S2.csv") %>%
  filter(Subset == 1) %>%
  set_names(c("pigID", "segment", "DPS1", "DPS2", "freq1", "freq2", "Subset")) %>%
  mutate(gen = (DPS2 - DPS1)*(24/hourPerGen))

# 5. Generate data frame ====
## 5.1. 2 hour per generation ====
hourPerGen <- 2
McCrone_IAV_2_per_gen <- McCrone_IAV %>%
  mutate(gen = (DPS2 - DPS1)*(24/hourPerGen))

LLdf_McCrone_2_per_gen <- LL_Ne.df(McCrone_IAV_2_per_gen, Ne_values, mu, alpha, hourPerGen, m)
MLE_Ne_McCrone_2_per_gen <- LLdf_McCrone_2_per_gen$Ne[which.max(LLdf_McCrone_2_per_gen$prob)]
CI_low_Ne_2_per_gen <- LLdf_McCrone_2_per_gen %>%
  filter(Ne < MLE_Ne_McCrone_2_per_gen) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_2_per_gen <- LLdf_McCrone_2_per_gen %>%
  filter(Ne > MLE_Ne_McCrone_2_per_gen) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

## 5.2. 4 hour per generation ====
hourPerGen <- 4
McCrone_IAV_4_per_gen <- McCrone_IAV %>%
  mutate(gen = (DPS2 - DPS1)*(24/hourPerGen))

LLdf_McCrone_4_per_gen <- LL_Ne.df(McCrone_IAV_4_per_gen, Ne_values, mu, alpha, hourPerGen, m)
MLE_Ne_McCrone_4_per_gen <- LLdf_McCrone_4_per_gen$Ne[which.max(LLdf_McCrone_4_per_gen$prob)]
CI_low_Ne_4_per_gen <- LLdf_McCrone_4_per_gen %>%
  filter(Ne < MLE_Ne_McCrone_4_per_gen) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_4_per_gen <- LLdf_McCrone_4_per_gen %>%
  filter(Ne > MLE_Ne_McCrone_4_per_gen) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

## 5.3. 6 hour per generation ====
hourPerGen <- 6
McCrone_IAV_6_per_gen <- McCrone_IAV %>%
  mutate(gen = (DPS2 - DPS1)*(24/hourPerGen))

LLdf_McCrone_6_per_gen <- LL_Ne.df(McCrone_IAV_6_per_gen, Ne_values, mu, alpha, hourPerGen, m)
MLE_Ne_McCrone_6_per_gen <- LLdf_McCrone_6_per_gen$Ne[which.max(LLdf_McCrone_6_per_gen$prob)]
CI_low_Ne_6_per_gen <- LLdf_McCrone_6_per_gen %>%
  filter(Ne < MLE_Ne_McCrone_6_per_gen) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_6_per_gen <- LLdf_McCrone_6_per_gen %>%
  filter(Ne > MLE_Ne_McCrone_6_per_gen) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

hourPerGen <- 6
swine_6_per_gen <- swineData %>%
  mutate(gen = (DPS2 - DPS1)*(24/hourPerGen))

LLdf_swine_6_per_gen <- LL_Ne.df(swine_6_per_gen, seq(3, 60, 1), mu, alpha, hourPerGen, m)
MLE_Ne_swine_6_per_gen <- LLdf_swine_6_per_gen$Ne[which.max(LLdf_swine_6_per_gen$prob)]
CI_low_Ne_swine_6_per_gen <- LLdf_swine_6_per_gen %>%
  filter(Ne < MLE_Ne_swine_6_per_gen) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_swine_6_per_gen <- LLdf_swine_6_per_gen %>%
  filter(Ne > MLE_Ne_swine_6_per_gen) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

## 5.4. 12 hour per generation ====
hourPerGen <- 12
McCrone_IAV_12_per_gen <- McCrone_IAV %>%
  mutate(gen = (DPS2 - DPS1)*(24/hourPerGen))

LLdf_McCrone_12_per_gen <- LL_Ne.df(McCrone_IAV_12_per_gen, Ne_values, mu, alpha, hourPerGen, m)
MLE_Ne_McCrone_12_per_gen <- LLdf_McCrone_12_per_gen$Ne[which.max(LLdf_McCrone_12_per_gen$prob)]
CI_low_Ne_12_per_gen <- LLdf_McCrone_12_per_gen %>%
  filter(Ne < MLE_Ne_McCrone_12_per_gen) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_12_per_gen <- LLdf_McCrone_12_per_gen %>%
  filter(Ne > MLE_Ne_McCrone_12_per_gen) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

hourPerGen <- 12
swine_12_per_gen <- swineData %>%
  mutate(gen = (DPS2 - DPS1)*(24/hourPerGen))

LLdf_swine_12_per_gen <- LL_Ne.df(swine_12_per_gen, seq(3, 60, 1), mu, alpha, hourPerGen, m)
MLE_Ne_swine_12_per_gen <- LLdf_swine_12_per_gen$Ne[which.max(LLdf_swine_12_per_gen$prob)]
CI_low_Ne_swine_12_per_gen <- LLdf_swine_12_per_gen %>%
  filter(Ne < MLE_Ne_swine_12_per_gen) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_swine_12_per_gen <- LLdf_swine_12_per_gen %>%
  filter(Ne > MLE_Ne_swine_12_per_gen) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

## 5.5. 24 hour per generation ====
hourPerGen <- 24
McCrone_IAV_24_per_gen <- McCrone_IAV %>%
  mutate(gen = (DPS2 - DPS1)*(24/hourPerGen))

LLdf_McCrone_24_per_gen <- LL_Ne.df(McCrone_IAV_24_per_gen, Ne_values, mu, alpha, hourPerGen, m)
MLE_Ne_McCrone_24_per_gen <- LLdf_McCrone_24_per_gen$Ne[which.max(LLdf_McCrone_24_per_gen$prob)]
CI_low_Ne_24_per_gen <- LLdf_McCrone_24_per_gen %>%
  filter(Ne < MLE_Ne_McCrone_24_per_gen) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_24_per_gen <- LLdf_McCrone_24_per_gen %>%
  filter(Ne > MLE_Ne_McCrone_24_per_gen) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

## 5.6. 3 hour per generation ====
hourPerGen <- 3
McCrone_IAV_3_per_gen <- McCrone_IAV %>%
  mutate(gen = (DPS2 - DPS1)*(24/hourPerGen))

LLdf_McCrone_3_per_gen <- LL_Ne.df(McCrone_IAV_3_per_gen, Ne_values, mu, alpha, hourPerGen, m)
MLE_Ne_McCrone_3_per_gen <- LLdf_McCrone_3_per_gen$Ne[which.max(LLdf_McCrone_3_per_gen$prob)]
CI_low_Ne_3_per_gen <- LLdf_McCrone_3_per_gen %>%
  filter(Ne < MLE_Ne_McCrone_3_per_gen) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_3_per_gen <- LLdf_McCrone_3_per_gen %>%
  filter(Ne > MLE_Ne_McCrone_3_per_gen) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

hourPerGen <- 3
swine_3_per_gen <- swineData %>%
  mutate(gen = (DPS2 - DPS1)*(24/hourPerGen))

LLdf_swine_3_per_gen <- LL_Ne.df(swine_3_per_gen, seq(3, 60, 1), mu, alpha, hourPerGen, m)
MLE_Ne_swine_3_per_gen <- LLdf_swine_3_per_gen$Ne[which.max(LLdf_swine_3_per_gen$prob)]
CI_low_Ne_swine_3_per_gen <- LLdf_swine_3_per_gen %>%
  filter(Ne < MLE_Ne_swine_3_per_gen) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)
CI_high_Ne_swine_3_per_gen <- LLdf_swine_3_per_gen %>%
  filter(Ne > MLE_Ne_swine_3_per_gen) %>%
  slice(which.min(abs(prob - (max(prob)-1.92)))) %>%
  select(Ne)

## 5.7. Comebine data ====
comb_ll_gen_sens <- bind_rows(
  LLdf_McCrone_2_per_gen  %>% mutate(gen_label = "2 hours"),
  LLdf_McCrone_4_per_gen  %>% mutate(gen_label = "4 hours"),
  LLdf_McCrone_6_per_gen  %>% mutate(gen_label = "6 hours"),
  LLdf_McCrone_12_per_gen %>% mutate(gen_label = "12 hours"),
  LLdf_McCrone_24_per_gen %>% mutate(gen_label = "24 hours")
) %>%
  mutate(gen_label = factor(gen_label, levels = c("2 hours", "4 hours", "6 hours", "12 hours", "24 hours")))

gen_sens_MLE <- data.frame(gen_label = factor(c("2 hours", "4 hours", "6 hours", "12 hours", "24 hours")),
                           Ne = c(MLE_Ne_McCrone_2_per_gen, 
                                  MLE_Ne_McCrone_4_per_gen,
                                  MLE_Ne_McCrone_6_per_gen,
                                  MLE_Ne_McCrone_12_per_gen,
                                  MLE_Ne_McCrone_24_per_gen),
                           prob = c(max(LLdf_McCrone_2_per_gen$prob),
                                    max(LLdf_McCrone_4_per_gen$prob),
                                    max(LLdf_McCrone_6_per_gen$prob),
                                    max(LLdf_McCrone_12_per_gen$prob),
                                    max(LLdf_McCrone_24_per_gen$prob)))

comb_ll_gen_sens_3612 <- bind_rows(
  LLdf_McCrone_3_per_gen  %>% mutate(gen_label = "3 hours"),
  LLdf_McCrone_6_per_gen  %>% mutate(gen_label = "6 hours"),
  LLdf_McCrone_12_per_gen %>% mutate(gen_label = "12 hours")
) %>%
  mutate(gen_label = factor(gen_label, levels = c("3 hours", "6 hours", "12 hours")))

gen_sens_MLE_3612 <- data.frame(gen_label = factor(c("3 hours", "6 hours", "12 hours")),
                           Ne = c(MLE_Ne_McCrone_3_per_gen, 
                                  MLE_Ne_McCrone_6_per_gen,
                                  MLE_Ne_McCrone_12_per_gen),
                           prob = c(max(LLdf_McCrone_3_per_gen$prob),
                                    max(LLdf_McCrone_6_per_gen$prob),
                                    max(LLdf_McCrone_12_per_gen$prob)))

comb_ll_gen_sens_3612_swine <- bind_rows(
  LLdf_swine_3_per_gen  %>% mutate(gen_label = "3 hours"),
  LLdf_swine_6_per_gen  %>% mutate(gen_label = "6 hours"),
  LLdf_swine_12_per_gen %>% mutate(gen_label = "12 hours")
) %>%
  mutate(gen_label = factor(gen_label, levels = c("3 hours", "6 hours", "12 hours")))

gen_sens_MLE_3612_swine <- data.frame(gen_label = factor(c("3 hours", "6 hours", "12 hours")),
                                Ne = c(MLE_Ne_swine_3_per_gen, 
                                       MLE_Ne_swine_6_per_gen,
                                       MLE_Ne_swine_12_per_gen),
                                prob = c(max(LLdf_swine_3_per_gen$prob),
                                         max(LLdf_swine_6_per_gen$prob),
                                         max(LLdf_swine_12_per_gen$prob)))

# 6. Graphs ====
fig_gen_sens <- ggplot() +
  geom_line(data = comb_ll_gen_sens, aes(x = Ne, y = prob, color = gen_label)) +
  geom_point(data = gen_sens_MLE, aes(x = Ne, y = prob), color = "red") +
  geom_vline(xintercept = MLE_Ne_McCrone, color = "red") +
  annotate("rect", xmin = as.numeric(CI_low_Ne), xmax = as.numeric(CI_high_Ne), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  scale_color_manual(
    name = "Generation Time", 
    values = c(
      "2 hours"  = "#26185F",
      "4 hours"  = "#005396",
      "6 hours"  = "#0088B1",
      "12 hours" = "#65CCB3",
      "24 hours" = "#AEE4C1"
    )
  ) +
  labs(x = expression(N[E]), y = "Log-likelihood") +
  scale_x_continuous(breaks = seq(0, 500, 20)) +
  ylim(-150, -30)

fig_gen_sens_3612_human <- ggplot() +
  geom_line(data = comb_ll_gen_sens_3612, aes(x = Ne, y = prob, color = gen_label)) +
  geom_point(data = gen_sens_MLE_3612, aes(x = Ne, y = prob), color = "red") +
  geom_vline(xintercept = MLE_Ne_McCrone, color = "red") +
  annotate("rect", xmin = as.numeric(CI_low_Ne_McCrone), xmax = as.numeric(CI_high_Ne_McCrone), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  scale_color_manual(
    name = "Generation Time", 
    values = c(
      "3 hours"  = "#26185F",
      "6 hours"  = "#0088B1",
      "12 hours" = "#AEE4C1"
    )
  ) +
  labs(x = expression(N[E]), y = "Log-likelihood", title = "Human data subset 1") +
  scale_x_continuous(breaks = seq(0, 500, 50))

fig_gen_sens_3612_swine <- ggplot() +
  geom_line(data = comb_ll_gen_sens_3612_swine, aes(x = Ne, y = prob, color = gen_label)) +
  geom_point(data = gen_sens_MLE_3612_swine, aes(x = Ne, y = prob), color = "red") +
  geom_vline(xintercept = 14, color = "red") +
  annotate("rect", xmin = as.numeric(CI_low_Ne_swineData), xmax = as.numeric(CI_high_Ne_swineData), 
           ymin = -Inf, ymax = Inf, alpha = 0.1,fill = "black") +
  scale_color_manual(
    name = "Generation Time", 
    values = c(
      "3 hours"  = "#26185F",
      "6 hours"  = "#0088B1",
      "12 hours" = "#AEE4C1"
    )
  ) +
  labs(x = expression(N[E]), y = "Log-likelihood", title = "Swine data subset 1") +
  scale_x_continuous(breaks = seq(0, 500, 10))

# 7. Save figure ====
fig_gen_sens_3612 <- fig_gen_sens_3612_human + fig_gen_sens_3612_swine + 
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    legend.position = "bottom"
  )

ggsave(fig_gen_sens, file = "human_ds1_gen_sens.pdf", path = "02_plots", width = 7, height = 4)
ggsave(fig_gen_sens_3612, file = "fig_gen_sens.pdf", path = "02_plots", width = 10, height = 5)
