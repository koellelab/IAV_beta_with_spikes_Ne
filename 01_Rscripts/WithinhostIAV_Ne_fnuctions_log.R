log10_dbeta <- function(x, alpha, beta) {
  log10_x_part   <- (alpha - 1) * log10(x)
  log10_1mx_part <- (beta  - 1) * log10(1 - x)
  log10_Beta <- (lgamma(alpha) + lgamma(beta) - lgamma(alpha + beta)) / log(10)
  log10_density <- log10_x_part + log10_1mx_part - log10_Beta
  return(log10_density)
}

log10_sum <- function(log10_values) {
  if (length(log10_values) == 0) return(-Inf)  # log10(0)
  m <- max(log10_values)
  if (!is.finite(m)) return(m)                 # all -Inf -> return -Inf
  m + log10(sum(10^(log10_values - m)))
}

beta_spike_approx <- function(Ne, ini_p, gen, mu, freq_step, v) {
  P_loss <- c(as.numeric(ini_p == 0))
  P_fix <- c(as.numeric(ini_p == 1))
  P_not <- c(1 - P_loss - P_fix)
  E <- c(ini_p)
  V <- v
  # Variance set to a small number in order for the alpha to be exist. 
  # Can be considered as measurement error
  
  if (P_loss[1] == 1){
    P_loss[2] <- (1-mu)^Ne
    P_fix[2] <- mu^Ne
  } else if (P_fix[1] == 1) {
    P_loss[2] <- mu^Ne
    P_fix[2] <- (1-mu)^Ne
  } else {
    cond_E <- c(ini_p)
    cond_V <- c((V+E^2-P_fix)/P_not - cond_E^2)
    alpha <- c(((cond_E*(1-cond_E))/cond_V - 1)*cond_E) 
    beta <- c(((cond_E*(1-cond_E))/cond_V - 1)*(1-cond_E))
    
    P_loss[2] <- P_loss[1]*(1-mu)^Ne + P_fix[1]*mu^Ne + 
      P_not[1]*(1-2*mu)^Ne*beta(alpha[1], beta[1]+Ne)/beta(alpha[1], beta[1])
    P_fix[2] <- P_loss[1]*mu^Ne + P_fix[1]*(1-mu)^Ne + 
      P_not[1]*(1-2*mu)^Ne*beta(alpha[1]+Ne, beta[1])/beta(alpha[1], beta[1])
  }
  
  # generation 1
  E[2] <- 1/2 + (1-2*mu) * (ini_p-1/2)
  V[2] <- 1/4 * (1-(1-2*mu)^2*(1-1/Ne))/(Ne-(1-2*mu)^2*(Ne-1)) - 
    (ini_p-1/2)^2 * (1-2*mu)^2 * (1-(1-1/Ne))
  P_not[2] <- 1 - P_loss[2] - P_fix[2]
  cond_E[2] <- (E[2]-P_fix[2])/P_not[2]
  cond_V[2] <- (V[2]+E[2]^2-P_fix[2])/P_not[2] - cond_E[2]^2
  alpha[2] <- ((cond_E[2]*(1-cond_E[2]))/cond_V[2] - 1)*cond_E[2]
  beta[2] <- ((cond_E[2]*(1-cond_E[2]))/cond_V[2] - 1)*(1-cond_E[2])
  
  # generation 2 and after
  if (gen >= 2){
    for (t in 2:gen) {
      E[t+1] <- 1/2 + (1-2*mu)^t * (ini_p-1/2)
      V[t+1] <- 1/4 * (1-(1-2*mu)^(2*t)*(1-1/Ne)^t)/(Ne-(1-2*mu)^2*(Ne-1)) - 
        (ini_p-1/2)^2 * (1-2*mu)^(2*t) * (1-(1-1/Ne)^t)
      
      P_loss[t+1] <- P_loss[t]*(1-mu)^Ne + P_fix[t]*mu^Ne + 
        P_not[t]*(1-2*mu)^Ne * beta(alpha[t], beta[t]+Ne)/beta(alpha[t], beta[t])
      P_fix[t+1] <- P_loss[t]*mu^Ne + P_fix[t]*(1-mu)^Ne + 
        P_not[t]*(1-2*mu)^Ne*beta(alpha[t]+Ne, beta[t])/beta(alpha[t], beta[t])
      P_not[t+1] <- 1 - P_loss[t+1] - P_fix[t+1]
      
      cond_E[t+1] <- (E[t+1]-P_fix[t+1])/P_not[t+1]
      cond_V[t+1] <- (V[t+1]+E[t+1]^2-P_fix[t+1])/P_not[t+1] - cond_E[t+1]^2
      
      alpha[t+1] <- ((cond_E[t+1]*(1-cond_E[t+1]))/cond_V[t+1] - 1)*cond_E[t+1]
      beta[t+1] <- ((cond_E[t+1]*(1-cond_E[t+1]))/cond_V[t+1] - 1)*(1-cond_E[t+1])
    }
  } else {
    t <- 1
  }
  
  df <- data.frame(freq = seq(0, 1, by = freq_step), generation = t)
  df$generation <- as.factor(df$generation)
  
  log10_beta_part <- (alpha[t+1] - 1) * log10(df$freq) +
    (beta[t+1] - 1) * log10(1 - df$freq) -
    log10(beta(alpha[t+1], beta[t+1]))
  
  log10_dens <- log10(P_not[t+1]) + log10_beta_part
  df$log10_dens <- log10_dens
  
  df_prob <- df %>%
    mutate(midpoints = df$freq - freq_step / 2, # Calculate midpoints
           midpoints = ifelse(midpoints < 0, 0, midpoints), # Handle edge case for 0
           log10_density = log10_dbeta(midpoints, alpha[t+1], beta[t+1]),
           log10_prob_mass = log10(P_not[t+1]) + log10_density + log10(freq_step))
  
  # Assign discrete probabilities in log10 space
  df$log10_prob <- NA
  # Continuous part (excluding first bin)
  df$log10_prob[2:(nrow(df))] <- df_prob$log10_prob_mass[-1]
  # First bin: P_loss
  df$log10_prob[1] <- log10(P_loss[t+1])
  # Last bin: P_fix
  if (df$freq[nrow(df)] == 1) {
    df$log10_prob[nrow(df)] <- log10(P_fix[t+1])
  } else {
    P_fix_row <- data.frame(
      freq = 1,
      generation = t,
      log10_prob = log10(P_fix[t+1]),
      log10_dens = -Inf
    )
    df <- rbind(df, P_fix_row)
  }
  
  df$dens <- P_not[t+1] * df$freq^(alpha[t+1]-1)*(1-df$freq)^(beta[t+1]-1)/beta(alpha[t+1], beta[t+1])
  df$prob <- 10^df$log10_prob
  
  df_params <- data.frame(generation = 0:gen, E = E, V = V, 
                          cond_E = cond_E, cond_V = cond_V, 
                          P_loss = P_loss, P_fix = P_fix, P_not_fix_loss = P_not, 
                          shape_alpha = alpha, shape_beta = beta)
  list <- list(df, df_params)
  return(list)
}


beta_spike_threshold <- function(Ne, ini_p, gen, mu, freq_step, v, alpha) {
  df <- beta_spike_approx(Ne, ini_p, gen, mu, freq_step, v)[[1]]

  # Merge left tail into x=0 bin in log space
  df$log10_dens[1] <- log10_sum(df$log10_dens[df$freq <= alpha])
  df$log10_prob[1] <- log10_sum(df$log10_prob[df$freq <= alpha])
  
  # Merge right tail into x=1 bin in log space
  df$log10_dens[nrow(df)] <- log10_sum(df$log10_dens[df$freq >= (1 - alpha)])
  df$log10_prob[nrow(df)] <- log10_sum(df$log10_prob[df$freq >= (1 - alpha)])
  
  df$prob <- 10^(df$log10_prob)
  df$cdf  <- cumsum(df$prob) / sum(df$prob)
  
  # Trim interior bins strictly between (0, alpha) and (1-alpha, 1)
  df$threshold <- alpha
  loc_0 <- max(as.numeric(rownames(df[df$freq < alpha, ])))
  if (nrow(df[df$freq > (1-alpha), ]) == 0) {
    df <- df[-c(2:loc_0), ]
  } else {
    loc_1 <- min(as.numeric(rownames(df[df$freq > (1-alpha), ])))
    df <- df[-c(2:loc_0, loc_1:(nrow(df)-1)), ]
  }
  return(df)
}

LL_Ne_mu <- function(df, Ne, mu, alpha, hourPerGen, m) {
  df$gen <- (df$DPS2 - df$DPS1)*(24/hourPerGen)
  df <- filter(df, gen != 0)  # exclude those who has 0 generations
  log_likelihood <- 0
  for (i in 1:nrow(df)) {
    ini_p <- df$freq1[i]
    fin_p <- df$freq2[i]
    gen <- df$gen[i]
    
    # determine the appropriate step for frequency
    for (j in 1:100) {
      if (fin_p == 0) {
        freq_step <- 0.01
      } else if (fin_p/j < 0.01) {
        freq_step <- fin_p/j
        break
      }
    }
    
    v <- m*ini_p*(1-ini_p)
    fin_vec <- beta_spike_threshold(Ne, ini_p, gen, mu, freq_step, v, alpha)
    log_prob <- fin_vec[which.min(abs(fin_vec$freq - fin_p)),"log10_prob"]
    log_likelihood <- log_likelihood + log_prob
  }
  return(log_likelihood)
}

LL_Ne.df <- function(df, Ne_values, mu, alpha, hourPerGen, m) {
  prob <- c()
  for (Ne in Ne_values) {
    prob <- c(prob, LL_Ne_mu(df, Ne, mu, alpha, hourPerGen, m))
  }
  LLdf <- data.frame(Ne = Ne_values, prob = prob)
  return(LLdf)
}

LL_NeNs <- function(df, Ne, Ns, hourPerGen) {
  df$gen <- (df$DPS2 - df$DPS1) * (24 / hourPerGen)
  q0 <- df$freq1
  dq <- df$freq2 - df$freq1
  
  sd_dq <- sqrt((q0 * (1 - q0)) * (df$gen / Ne + 1 / Ns))
  
  sum_ll <- sum(dnorm(dq, mean = 0, sd = sd_dq, log = TRUE))
  
  return(sum_ll)
}

LL_Ns_noise <- function(df, Ns, hourPerGen) {
  df$gen <- (df$DPS2 - df$DPS1) * (24 / hourPerGen)
  q0 <- df$freq1
  dq <- df$freq2 - df$freq1
  
  sd_dq <- sqrt((q0 * (1 - q0)) * (1 / Ns))
  
  sum_ll <- sum(dnorm(dq, mean = 0, sd = sd_dq, log = TRUE))
  
  return(sum_ll)
}

LL_NeNs_surface <- function(df, Ne_values, Ns_values, hourPerGen) {
  grid_df <- expand.grid(
    Ne = Ne_values,
    Ns = Ns_values
  ) %>%
    as_tibble()
  
  # Compute surface
  grid_df$logLik <- NA_real_
  for (k in seq_len(nrow(grid_df))) {
    grid_df$logLik[k] <- LL_NeNs(df, grid_df$Ne[k], grid_df$Ns[k], hourPerGen)
  }
  surface_df <- grid_df
  
  return(surface_df)
}







adjust_data <- function(df) {
  df_obs1 <- select(df, DPS1, freq1)
  df_obs1$type <- c("obs1")
  names(df_obs1) <- c("DPS", "freq", "type")
  df_obs1$ID <- (1:nrow(df_obs1))
  # data frame only including second round of observations
  df_obs2 <- select(df, DPS2, freq2)
  df_obs2$type <- c("obs2")
  names(df_obs2) <- c("DPS", "freq", "type")
  df_obs2$ID <- (1:nrow(df_obs2))
  # merge the data together
  df_tot_obs <- rbind(df_obs1, df_obs2)
  df_tot_obs$type <- factor(df_tot_obs$type, labels = c("observation 1", "observation 2"))
  return(df_tot_obs)
}

find_overlapping_pairs <- function(row) {
  non_na_indices <- which(!is.na(row))
  pairs_list <- list()
  
  if (length(non_na_indices) >= 2) {
    for (i in 1:(length(non_na_indices) - 1)) {
      first_index <- non_na_indices[i] - 1
      second_index <- non_na_indices[i + 1] - 1
      first_value <- row[first_index + 1]
      second_value <- row[second_index + 1]
      pairs_list[[length(pairs_list) + 1]] <- c(first_index, second_index, first_value, second_value)
    }
  }
  
  return(do.call(rbind, pairs_list))
}

each_iSNV_LL <- function(df, Ne_values, mu, alpha, hourPerGen, m) {
  list_LLs <- list()
  for (i in 1:nrow(df)) {
    iSNV_pair <- df[i,]
    list_LLs[[i]] <- LL_Ne.df(iSNV_pair, Ne_values, mu, alpha, hourPerGen, m)
  }
  df_LLs <- bind_rows(list_LLs, .id = "iSNV_id")
  return(df_LLs)
}






generate_random_dataset <- function(original_data, random_seed) {
  set.seed(random_seed)
  data_with_index <- original_data %>% mutate(index = row_number())
  sub1 <- data_with_index %>%
    filter(gen > 0, freq1 >= threshold) %>%
    group_by(ID) %>%
    slice_sample(n = 1) %>%
    ungroup()
  sub2 <- data_with_index %>%
    filter(gen > 0, freq1 < threshold, freq2 >= threshold) %>%
    group_by(ID) %>%
    slice_sample(n = 1) %>%
    ungroup()
  final_data <- data_with_index %>%
    mutate(Subset = case_when(
      index %in% sub1$index ~ 1,
      index %in% sub2$index ~ 2,
      TRUE ~ NA_real_
    )) %>%
    select(ID, iSNV, DPS1, DPS2, freq1, freq2, Subset) %>%
    dplyr::rename(
      "Individual ID" = ID,
      "DPS1 Frequency" = freq1,
      "DPS2 Frequency" = freq2
    )
  return(final_data)
}

make_traj_plot <- function(df, index_label) {
  df_clean <- df %>%
    rename(
      ENROLLID = `Individual ID`,
      mutation = `iSNV`,
      freq1 = `DPS1 Frequency`,
      freq2 = `DPS2 Frequency`
    ) %>% 
    mutate(gen = (DPS2 - DPS1)*(24/hourPerGen)) %>%
    filter(Subset == 1)
  d1 <- df_clean %>% 
    select(mutation, ENROLLID, DPS1, freq1, gen) %>% 
    rename(DPS = DPS1, freq = freq1) %>% 
    mutate(type = "obs1")
  d2 <- df_clean %>% 
    select(mutation, ENROLLID, DPS2, freq2, gen) %>% 
    rename(DPS = DPS2, freq = freq2) %>% 
    mutate(type = "obs2")
  tot <- bind_rows(d1, d2) %>%
    mutate(
      type_label = case_when(
        type == "obs2" & freq < 0.02 ~ "observation 2 below 2%",
        type == "obs1" ~ "observation 1",
        type == "obs2" ~ "observation 2"
      ),
      freq = if_else(type == "obs2" & freq < 0.02, 0.02, freq)
    )
  p <- ggplot(tot, aes(x = DPS, y = freq)) +
    geom_hline(yintercept = 0.02, color = "red", linetype = "dashed") +
    geom_point(aes(color = type_label, shape = type_label), alpha = 0.7, size = 2) +
    geom_line(aes(group = mutation), alpha = 0.3) + # Group by mutation to connect lines
    scale_shape_manual(values = c("observation 1" = 19, "observation 2" = 17, "observation 2 below 2%" = 2)) +
    scale_color_manual(values = c("observation 1" = "#26185F", "observation 2" = "#0095AF", "observation 2 below 2%" = "black")) +
    labs(x = "DPS ", y = "iSNV frequency") +
    ylim(0, 0.7) +
    theme(legend.position = "none")
  
  return(p)
}

calc_MLE <- function(df, index_label) {
  df_clean <- df %>%
    rename(
      ENROLLID = `Individual ID`,
      mutation = `iSNV`,
      freq1 = `DPS1 Frequency`,
      freq2 = `DPS2 Frequency`
    ) %>% 
    mutate(gen = (DPS2 - DPS1)*(24/hourPerGen)) %>%
    filter(Subset == 1)
  LLdf <- LL_Ne.df(df_clean, Ne_values, mu, alpha, hourPerGen, m)
  MLE_Ne <- LLdf$Ne[which.max(LLdf$prob)]
  max_prob <- max(LLdf$prob)
  CI_low_row <- LLdf %>% filter(Ne < MLE_Ne)
  CI_low <- if(nrow(CI_low_row) > 0) CI_low_row$Ne[which.min(abs(CI_low_row$prob - (max_prob - 1.92)))] else NA
  CI_high_row <- LLdf %>% filter(Ne > MLE_Ne)
  CI_high <- if(nrow(CI_high_row) > 0) CI_high_row$Ne[which.min(abs(CI_high_row$prob - (max_prob - 1.92)))] else NA
  list <- list(LLdf, c(max_prob, CI_low, CI_high))
  return(p)
}

make_like_plot <- function(df, index_label) {
  df_clean <- df %>%
    rename(
      ENROLLID = `Individual ID`,
      mutation = `iSNV`,
      freq1 = `DPS1 Frequency`,
      freq2 = `DPS2 Frequency`
    ) %>% 
    mutate(gen = (DPS2 - DPS1)*(24/hourPerGen)) %>%
    filter(Subset == 1)
  LLdf <- LL_Ne.df(df_clean, Ne_values, mu, alpha, hourPerGen, m)
  MLE_Ne <- LLdf$Ne[which.max(LLdf$prob)]
  max_prob <- max(LLdf$prob)
  CI_low_row <- LLdf %>% filter(Ne < MLE_Ne)
  CI_low <- if(nrow(CI_low_row) > 0) CI_low_row$Ne[which.min(abs(CI_low_row$prob - (max_prob - 1.92)))] else NA
  CI_high_row <- LLdf %>% filter(Ne > MLE_Ne)
  CI_high <- if(nrow(CI_high_row) > 0) CI_high_row$Ne[which.min(abs(CI_high_row$prob - (max_prob - 1.92)))] else NA
  
  p <- ggplot(LLdf) +
    geom_line(aes(x = Ne, y = prob)) +
    geom_vline(xintercept = MLE_Ne, color = "red") +
    geom_hline(yintercept = max_prob, color = "red", linetype = "dashed") +
    annotate("rect", xmin = as.numeric(CI_low), xmax = as.numeric(CI_high), 
             ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "black") +
    annotate("text", x = 100, y = min(LLdf$prob) + 2, 
             label = paste("N[E]==", MLE_Ne), parse = TRUE, color = "red", size = 3) +
    labs(x = expression(N[E]), y = "Log-likelihood") +
    scale_x_continuous(breaks = seq(0, 200, 50))
  
  return(p)
}
