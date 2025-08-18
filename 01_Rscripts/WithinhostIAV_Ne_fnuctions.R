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
  df$dens <- P_not[t+1] * df$freq^(alpha[t+1]-1)*(1-df$freq)^(beta[t+1]-1)/beta(alpha[t+1], beta[t+1])
  
  df_prob <- df %>%
    mutate(midpoints = df$freq - freq_step / 2, # Calculate midpoints
           midpoints = ifelse(midpoints < 0, 0, midpoints), # Handle edge case for 0
           density = dbeta(midpoints, alpha[t+1], beta[t+1]), # Density at midpoint
           probability = density * freq_step) # Probability for each interval
  
  df$prob[2:(nrow(df))] <- P_not[t+1] * df_prob$probability[-1]
  df$prob[1] <- P_loss[t+1]
  if (df$freq[nrow(df)] == 1) {
    df$prob[nrow(df)] <- P_fix[t+1]
  } else {
    P_fix_row <- data.frame(freq = 1, generation = t, dens = 0, prob = P_fix[t+1])
    df <- rbind(df, P_fix_row)
  }
  df_params <- data.frame(generation = 0:gen, E = E, V = V, 
                          cond_E = cond_E, cond_V = cond_V, 
                          P_loss = P_loss, P_fix = P_fix, P_not_fix_loss = P_not, 
                          shape_alpha = alpha, shape_beta = beta)
  list <- list(df, df_params)
  return(list)
}

beta_spike_threshold <- function(Ne, ini_p, gen, mu, freq_step, v, alpha) {
  df <- beta_spike_approx(Ne, ini_p, gen, mu, freq_step, v)[[1]]
  df$dens[1] <- sum(df[df$freq <= alpha, ]$dens)
  df$prob[1] <- sum(df[df$freq <= alpha, ]$prob)
  df$dens[nrow(df)] <- sum(df[df$freq >= (1-alpha), ]$dens)
  df$prob[nrow(df)] <- sum(df[df$freq >= (1-alpha), ]$prob)
  df$threshold <- alpha
  df$cdf <- cumsum(df$prob)
  max_cdf <- max(df$cdf)
  df$cdf <- df$cdf/max_cdf
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
  prob <- 1
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
    prob <- fin_vec[which.min(abs(fin_vec$freq - fin_p)),"prob"]
    if (prob < (.Machine$double.xmin * .Machine$double.eps)) {
      prob <- (.Machine$double.xmin * .Machine$double.eps)
    }
    log_likelihood <- log_likelihood + log10(prob)
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



