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

# 4. Input parameters ====
mu <- 0
alpha <- 0.02
hourPerGen <- 6
m <- 0.004
freq_step <- 0.01
ini_p <- 0.325
v <- m*ini_p*(1-ini_p)

# 5. Generate data frame ====
df_bs_gen2_ne20 <- beta_spike_approx(20, ini_p, 2, mu, freq_step, v)[[1]]
df_bs_gen4_ne20 <- beta_spike_approx(20, ini_p, 4, mu, freq_step, v)[[1]]
df_bs_gen6_ne20 <- beta_spike_approx(20, ini_p, 6, mu, freq_step, v)[[1]]
df_bs_gen0_ne20 <- df_bs_gen2_ne20
df_bs_gen0_ne20$prob <- dnorm(df_bs_gen0_ne20$freq, mean = ini_p, sd = sqrt(v))
df_bs_gen0_ne20$prob <- df_bs_gen0_ne20$prob/sum(df_bs_gen0_ne20$prob)
df_bs_gen0_ne20$generation <- 0
df_bs_ne20 <- rbind(df_bs_gen0_ne20, df_bs_gen2_ne20, df_bs_gen4_ne20, df_bs_gen6_ne20)

df_bs_gen2_ne50 <- beta_spike_approx(50, ini_p, 2, mu, freq_step, v)[[1]]
df_bs_gen4_ne50 <- beta_spike_approx(50, ini_p, 4, mu, freq_step, v)[[1]]
df_bs_gen6_ne50 <- beta_spike_approx(50, ini_p, 6, mu, freq_step, v)[[1]]
df_bs_gen0_ne50 <- df_bs_gen2_ne50
df_bs_gen0_ne50$prob <- dnorm(df_bs_gen0_ne50$freq, mean = ini_p, sd = sqrt(v))
df_bs_gen0_ne50$prob <- df_bs_gen0_ne50$prob/sum(df_bs_gen0_ne50$prob)
df_bs_gen0_ne50$generation <- 0
df_bs_ne50 <- rbind(df_bs_gen0_ne50, df_bs_gen2_ne50, df_bs_gen4_ne50, df_bs_gen6_ne50)

df_bs_gen2_ne500 <- beta_spike_approx(500, ini_p, 2, mu, freq_step, v)[[1]]
df_bs_gen4_ne500 <- beta_spike_approx(500, ini_p, 4, mu, freq_step, v)[[1]]
df_bs_gen6_ne500 <- beta_spike_approx(500, ini_p, 6, mu, freq_step, v)[[1]]
df_bs_gen0_ne500 <- df_bs_gen2_ne500
df_bs_gen0_ne500$prob <- dnorm(df_bs_gen0_ne500$freq, mean = ini_p, sd = sqrt(v))
df_bs_gen0_ne500$prob <- df_bs_gen0_ne500$prob/sum(df_bs_gen0_ne500$prob)
df_bs_gen0_ne500$generation <- 0
df_bs_ne500 <- rbind(df_bs_gen0_ne500, df_bs_gen2_ne500, df_bs_gen4_ne500, df_bs_gen6_ne500)

# 6. Graphs ====
fig1a <- ggplot(df_bs_ne20, aes(x = freq, y = prob)) +
  geom_line(aes(color = generation)) +
  geom_point(aes(color = generation)) +
  scale_color_discrete_sequential(palette = "YlGnBu", 
                                  nmax = 10, 
                                  order = c(4, 6, 8, 10),
                                  labels = c("generation 0", "generation 2", "generation 4", "generation 6")) +
  theme(legend.position = "none") +
  labs(y = "Probability", x = "Allele frequency", title = expression(N[E] == 20)) +
  ylim(0, 0.3)

fig1b <- ggplot(df_bs_ne50, aes(x = freq, y = prob)) +
  geom_line(aes(color = generation)) +
  geom_point(aes(color = generation)) +
  scale_color_discrete_sequential(palette = "YlGnBu", 
                                  nmax = 10, 
                                  order = c(4, 6, 8, 10)) +
  theme(legend.position = "none") +
  labs(y = "Probability", x = "Allele frequency", title = expression(N[E] == 50)) +
  ylim(0, 0.3)


fig1c <- ggplot(df_bs_ne500, aes(x = freq, y = prob)) +
  geom_line(aes(color = generation)) +
  geom_point(aes(color = generation)) +
  scale_color_discrete_sequential(palette = "YlGnBu", 
                                  nmax = 10, 
                                  order = c(4, 6, 8, 10),
                                  labels = c("generation 0", "generation 2", "generation 4", "generation 6")) +
  theme(legend.position = c(0.75, 0.85)) +
  labs(y = "Probability", x = "Allele frequency", color = NULL, title = expression(N[E] == 500)) +
  ylim(0, 0.3)

# 7. Put graphs together ====
fig1 <- fig1a + fig1b + fig1c + 
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 18))

ggsave(fig1, file = "fig1.pdf", path = "02_plots", width = 15, height = 4)
