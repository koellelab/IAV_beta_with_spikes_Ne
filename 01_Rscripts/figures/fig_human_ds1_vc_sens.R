# 1. Source functions ====
source("01_Rscripts/fig_human_ds1_05.R")
source("01_Rscripts/WithinhostIAV_Ne_fnuctions_log.R")

# 2. Load packages ====
library(tidyverse)
library(ggplot2)
library(patchwork)
library(showtext)
library(colorspace)
library(parallel)
library(scales)
theme_set(theme_minimal() +
            theme(axis.ticks = element_line(colour = "grey50"),
                  axis.line = element_line(colour = "grey50"),
                  axis.text = element_text(size = 12),
                  axis.title = element_text(size = 15),
                  legend.text = element_text(size = 12)))

# 3. Input parameters ====
af_values <- seq(0.002, 0.05, 0.0005)
min_dp <- 500
min_pos <- 200
Ne_values <- seq(10, 180, 1)
mu <- 0
hourPerGen <- 6
m <- 0.004
freq_step <- 0.01
set.seed(27)

# 4. Input datasets ====
vcf_wide <- read_csv('03_data/processed_vcf_wide_rep.csv', show_col_types=FALSE)
af_cols      <- c('d-2.0_AF','d0.0_AF','d1.0_AF','d2.0_AF','d3.0_AF','d4.0_AF','d5.0_AF','d6.0_AF')
af_raw_cols  <- c('d-2.0_AF_raw','d0.0_AF_raw','d1.0_AF_raw','d2.0_AF_raw','d3.0_AF_raw','d4.0_AF_raw','d5.0_AF_raw','d6.0_AF_raw')
depth_cols   <- c('depth_d-2.0','depth_d0.0','depth_d1.0','depth_d2.0','depth_d3.0','depth_d4.0','depth_d5.0','depth_d6.0')
dps_values   <- c(-2, 0:6)

# 5. Calculate MLE of NE for each threshold, using mclapply ====
run_one_threshold <- function(min_af, vcf_wide, Ne_values, mu, hourPerGen, m, min_dp, min_pos, dps_values) {
  keep <- rowSums(vcf_wide[af_cols] >= min_af, na.rm = TRUE) > 0
  keep <- keep & (
    rowSums(
      (vcf_wide[af_cols] == 0) & (vcf_wide[af_raw_cols] >= min_af),
      na.rm = TRUE
    ) == 0
  )
  keep <- keep & (
    rowSums(vcf_wide[depth_cols] < min_dp, na.rm = TRUE) == 0
  )
  keep <- keep & (vcf_wide$iSNV_loc > min_pos) & ((vcf_wide$segment_len - vcf_wide$iSNV_loc) > min_pos)
  keep_vcf <- vcf_wide[keep, , drop = FALSE] %>%
    dplyr::rename(
      `d-2.0` = `d-2.0_AF`,
      d0.0    = d0.0_AF,
      d1.0    = d1.0_AF,
      d2.0    = d2.0_AF,
      d3.0    = d3.0_AF,
      d4.0    = d4.0_AF,
      d5.0    = d5.0_AF,
      d6.0    = d6.0_AF
    ) %>%
    relocate(`d-2.0`, .before = `d0.0`) %>%
    relocate(`d6.0`, .before = `d0.0_AF_raw`) %>%
    relocate(`d3.0`, .before = `d4.0`) %>%
    relocate(`d-2.0_AF_raw`, .before = `d0.0_AF_raw`) %>%
    relocate(`d3.0_AF_raw`, .before = `d4.0_AF_raw`) %>%
    relocate(`d6.0_AF_raw`, .before = `depth_d-2.0`)
  
  # Human table construction
  flumin_origin <- keep_vcf[, c(1:9, 11:14, 16)]
  flumin_origin$index <- seq_len(nrow(flumin_origin))
  McCrone_IAV_flumina_data <- flumin_origin %>%
    filter(data == "mccrone") %>%
    rowwise() %>%
    mutate(
      DPS1 = dps_values[which(!is.na(c_across(c(`d-2.0`, d0.0:d6.0))))][1],
      DPS2 = dps_values[which(!is.na(c_across(c(`d-2.0`, d0.0:d6.0))))][2],
      gen  = (DPS2 - DPS1) * (24 / hourPerGen),
      iSNV = paste0(segment, "_", ref_nt, as.character(iSNV_loc), iSNV_nt),
      non_na_values = list(na.omit(c_across(c(`d-2.0`, d0.0:d6.0)))),
      freq1 = non_na_values[1],
      freq2 = non_na_values[2]
    ) %>%
    select(-non_na_values) %>%
    filter(!is.na(freq1), !is.na(freq2)) %>%
    ungroup()
  
  # Subset 2 (rare -> above threshold)
  McCrone_IAV_flumina_data_sub2 <- McCrone_IAV_flumina_data %>%
    filter(gen > 0, freq1 < min_af, freq2 >= min_af) %>%
    group_by(ID) %>%
    slice_sample(n = 1) %>%
    ungroup()
  
  # Assign subsets
  McCrone_IAV_flumina_data <- McCrone_IAV_flumina_data %>%
    mutate(Subset = if_else(index %in% McCrone_IAV_flumina_data_sub2$index, 2L, NA_integer_)) %>%
    group_by(ID) %>%
    mutate(Subset = if_else(freq1 >= min_af & gen > 0 & row_number() == which.min(abs(freq1 - 0.5)),
                            1L, Subset)) %>%
    ungroup()
  
  McCrone_IAV <- McCrone_IAV_flumina_data[, c("ID","iSNV","DPS1","DPS2","freq1","freq2","Subset")] %>%
    set_names(c("ENROLLID","mutation","DPS1","DPS2","freq1","freq2","Subset")) %>%
    mutate(gen = (DPS2 - DPS1) * (24 / hourPerGen)) %>%
    filter(Subset == 1)
  
  # MLE + CI
  LLdf_McCrone <- LL_Ne.df(McCrone_IAV, Ne_values, mu, min_af, hourPerGen, m)
  
  MLE_Ne <- LLdf_McCrone$Ne[which.max(LLdf_McCrone$prob)]
  maxLL  <- max(LLdf_McCrone$prob)
  target <- maxLL - 1.92
  
  CI_low_Ne <- LLdf_McCrone %>%
    filter(Ne < MLE_Ne) %>%
    slice(which.min(abs(prob - target))) %>%
    pull(Ne)
  
  CI_high_Ne <- LLdf_McCrone %>%
    filter(Ne > MLE_Ne) %>%
    slice(which.min(abs(prob - target))) %>%
    pull(Ne)
  
  out_row <- tibble(
    min_af = min_af,
    Ne_MLE = MLE_Ne,
    loglik = maxLL,
    CI_low = as.numeric(CI_low_Ne),
    CI_high = as.numeric(CI_high_Ne)
  )
  
  list(
    human_tbl = McCrone_IAV,
    mle_row   = out_row
  )
}

res_list <- parallel::mclapply(
  af_values,
  FUN = run_one_threshold,
  vcf_wide = vcf_wide,
  Ne_values = Ne_values,
  mu = mu,
  hourPerGen = hourPerGen,
  m = m,
  min_dp = min_dp,
  min_pos = min_pos,
  dps_values = dps_values,
  mc.cores = detectCores()
)

# Collect outputs
list_human_iSNV <- setNames(lapply(res_list, `[[`, "human_tbl"), paste0("af_", af_values))
table_vc_sens_MLE <- bind_rows(lapply(res_list, `[[`, "mle_row"))

# 6. Plot figures ====
fig_vc_sens <- ggplot(table_vc_sens_MLE, aes(x = min_af)) +
  geom_ribbon(
    aes(ymin = CI_low, ymax = CI_high),
    fill = "black",
    alpha = 0.1
  ) +
  geom_line(aes(y = Ne_MLE), color = "black", linewidth = 1) +
  annotate("rect", ymin = as.numeric(CI_low_Ne_McCrone), ymax = as.numeric(CI_high_Ne_McCrone), 
           xmin = 0.002, xmax = 0.05, alpha = 0.1, fill = "red") +
  geom_hline(yintercept = MLE_Ne_McCrone, color = "red") +
  scale_x_continuous(labels = percent_format(accuracy = 0.1)) +
  scale_y_continuous(breaks = seq(0, 200, 30)) +
  labs(
    x = "Variant-calling threshold",
    y = expression("MLE of " * N[E])
  ) 

fig_vc_sens_comb <- (fig2a_0.005 | fig2b_0.005) / fig_vc_sens +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 18))

ggsave(fig_vc_sens_comb, file = "fig_vc_sens.pdf", path = "02_plots", width = 10, height = 8)
