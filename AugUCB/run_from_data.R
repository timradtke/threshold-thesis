# pull_arm_from_tsdata

pull_arm_from_tsdata <- function(k, data, iter) {
  data[iter, k]
}

##############################################################

# General functions

get_min <- function(x) {
  mini <- min(x)
  ifelse(sum(x==mini) == 1, which.min(x),
         sample(which(x==mini)))
}

##############################################################

#mean_002 <- c(0.033, 0.033, 0.037, 0.037,
#              0.031, 0.031, 0.039, 0.039,
#              0.027, 0.027, 0.043, 0.043,
#              0.02, 0.02, 0.05, 0.05,
#              0.005, 0.005, 0.065, 0.065)
#tau_002 <- 0.035
#epsilon_002 <- 0.005#

#set.seed(512)
#testt <- data.frame(rep(NA, times = 500))
#for(i in 1:length(mean_002)) {
#  testt[[i]] <- as.numeric(purrr::rbernoulli(500, p  = mean_002[i]))
#}
#names(testt) <- paste0("V", 1:20)

#aug_test <- AugUCB_from_tsdata(testt, rounds = 500, tau = tau_002, 
#                               verbose = TRUE, seed = 46)

# AugUCB Algorithm

AugUCB_from_tsdata <- function(data, rounds = 5000, rho = 1/3, tau = NA,
                               verbose = FALSE, seed = NA) {
  if(!is.na(seed)) set.seed(seed)
  K <- dim(data)[2]
  
  # initialize by pulling each arm once
  arm_list <- diag(as.matrix(data[1:K,]))
  arm_list <- as.list(arm_list)
  
  # initialize mean storage, arm sequence, and the counter
  mean_storage <- matrix(unlist(lapply(arm_list, mean)), nrow = 1)
  arm_sequence <- 1:K
  
  # Initialize the remaining parameters as in the paper
  B <- 1:K # the active set of arms
  m <- 0
  epsilon <- 1
  e <- exp(1) # ????????
  M <- floor(0.5*log2(rounds/e))
  psi <- rounds*epsilon/(128 * (log(3/16 * K * log(K)))^2)
  l <- ceiling(2*psi*log(rounds * epsilon)/epsilon)
  N <- K * l
  
  if(verbose) counter <- K
  
  for(i in (K+1):rounds) {
    if(verbose) message(paste("this is round", i))
    if(length(B) > 0) {
      next_arm <- get_min(unlist(get_next_arm_augucb(arm_list, tau = tau, rho = rho,
                                                     psi = psi, rounds = rounds, 
                                                     epsilon = epsilon,
                                                     active_set = B)))
      arm_sequence <- c(arm_sequence, next_arm)
      if(verbose) message("arm selecting done")
      
      arm_list[[next_arm]] <- c(arm_list[[next_arm]], data[i, next_arm])
    }
    
    # Continue adding the same mean if the active set is empty
    mean_storage <- rbind(mean_storage, unlist(lapply(arm_list, mean)))
    if(verbose) message("arm pulling done")
    
    if(verbose) counter <- counter + 1
    if(verbose) message("Counter: ", counter)
    
    # delete arms should return a logical index for each arm in active set
    if(length(B) > 0) { # stop updating the active set if already empty
      B <- B[!unlist(delete_arms(arm_list, tau = tau, rho = rho,
                                 psi = psi, rounds = rounds, epsilon = epsilon,
                                 active_set = B))]
    }
    
    if(verbose) message("B done")
    if(verbose) message(B)
    
    if((i >= N) & (m <= M) & (length(B)>0)) {
      # reset parameters
      epsilon <- epsilon/2
      psi <- rounds * epsilon / (128*(log(3/16*K*log(K)))^2)
      l <- ceiling( (2 * psi * log(rounds*epsilon))/epsilon )
      N <- i + length(B) * l
      m <- m+1
    }
  }
  return(list(means = unlist(lapply(arm_list, mean)),
              arm_list = arm_list,
              #active_set = B,
              mean_storage = mean_storage,
              arm_sequence = arm_sequence))
}

get_next_arm_augucb <- function(armls, tau, rho, psi, rounds, epsilon, active_set) {
  get_metric <- function(x) {
    xhat <- mean(x)
    ni <- length(x)
    vhat <- ifelse(is.na(var(x)),0,var(x))*(ni-1)/ni
    si <- sqrt(rho * psi * (vhat+1) * log(rounds * epsilon)/4/ni)
    return(abs(xhat - tau)-2*si)
  }
  res <- lapply(armls, get_metric)
  res[-active_set] <- Inf
  return(res)
}

delete_arms <- function(armls, tau = tau, rho = rho,
                        psi = psi, rounds = rounds, epsilon = epsilon,
                        active_set = B) {
  eliminate <- function(x) {
    xhat <- mean(x)
    ni <- length(x)
    vhat <- ifelse(is.na(var(x)),0,var(x))*(ni-1)/ni
    si <- sqrt(rho * psi * (vhat+1) * log(rounds * epsilon)/4/ni)
    return((xhat + si < tau - si) || (xhat - si > tau + si))
  }
  return(lapply(armls[active_set], eliminate))
}

##############################################################

# APT Algorithm

get_next_arm_apt <- function(armls, tau, rounds, epsilon) {
  get_metric <- function(x) {
    xhat <- mean(x)
    ni <- length(x)
    return(sqrt(ni) * (abs(xhat-tau)+epsilon))
  }
  lapply(armls, get_metric)
}

APT_from_tsdata <- function(data, rounds = 5000, tau, epsilon,
                            verbose = FALSE, seed = NA) {
  
  if(!is.na(seed)) set.seed(seed)
  K <- dim(data)[2]
  
  # initialize by pulling each arm once
  arm_list <- diag(as.matrix(data[1:K,]))
  arm_list <- as.list(arm_list)
  
  # initialize mean storage, arm sequence, and the counter
  mean_storage <- matrix(unlist(lapply(arm_list, mean)), nrow = 1)
  arm_sequence <- 1:K
  
  for(i in (K+1):rounds) {
    if(verbose) message(paste("this is round", i))
    next_arm <- get_min(unlist(get_next_arm_apt(arm_list, tau = tau, 
                                                rounds = rounds, 
                                                epsilon = epsilon)))
    arm_sequence <- c(arm_sequence, next_arm)
    
    if(verbose) message("arm selecting done")
    arm_list[[next_arm]] <- c(arm_list[[next_arm]], data[i, next_arm])
    mean_storage <- rbind(mean_storage, unlist(lapply(arm_list, mean)))
    
    if(verbose) message("arm pulling done")
    
  }
  return(list(means = unlist(lapply(arm_list, mean)),
              arm_list = arm_list,
              input = data,
              arm_sequence = arm_sequence,
              mean_storage = mean_storage))
}

##############################################################

# KL-based Algorithm

KL_bandit_from_tsdata <- function(data, rounds = 5000, tau, epsilon,
                                  verbose = FALSE, seed = NA,
                                  at_tau = FALSE, horizon = NA) {
  
  if(!is.na(seed)) set.seed(seed)
  K <- dim(data)[2]
  
  # initialize by pulling each arm once
  arm_list <- diag(as.matrix(data[1:K,]))
  arm_list <- as.list(arm_list)
  
  # initialize mean storage and the counter
  mean_storage <- matrix(unlist(lapply(arm_list, mean)), nrow = 1)
  arm_sequence <- 1:K
  
  for(i in (K+1):rounds) {
    if(verbose) message(paste("this is round", i))
    if(at_tau) { # without epsilon interval
      next_arm <- get_min(unlist(get_next_arm_kl_at_tau(arm_list, tau = tau,
                                                        epsilon = epsilon,
                                                        current_round = i,
                                                        horizon = horizon)))
    } else {
      next_arm <- get_min(unlist(get_next_arm_kl(arm_list, tau = tau,
                                                 epsilon = epsilon,
                                                 horizon = horizon)))
    }
    arm_sequence <- c(arm_sequence, next_arm)
    
    if(verbose) message("arm selecting done")
    if(verbose) message(next_arm)
    arm_list[[next_arm]] <- c(arm_list[[next_arm]], data[i, next_arm])
    
    mean_storage <- rbind(mean_storage, unlist(lapply(arm_list, mean)))
    
    if(verbose) message("arm pulling done")
    
  }
  return(list(means = unlist(lapply(arm_list, mean)),
              arm_list = arm_list,
              arm_sequence = arm_sequence,
              mean_storage = mean_storage))
}

kl_ber <- function(x,y) {
  x*log(x/y) + (1-x)*log((1-x)/(1-y))
}

get_next_arm_kl <- function(armls, tau, epsilon, horizon) {
  get_metric <- function(x) {
    if(is.na(horizon)) {
      xhat <- mean(x)
    } else {
      xhat <- ifelse(mean(x) == 1, 1-1/horizon,
                     ifelse(mean(x) == 0, 1/horizon, mean(x)))
    }
    ni <- length(x)
    KLhat <- ifelse(xhat >= tau,
                    kl_ber(xhat, tau-epsilon),
                    kl_ber(xhat, tau+epsilon))
    return(ni * ifelse(is.na(KLhat),0,KLhat))
  }
  lapply(armls, get_metric)
}

get_next_arm_kl_at_tau <- function(armls, tau, epsilon, current_round,
                                   horizon) {
  get_metric <- function(x) {
    if(is.na(horizon)) {
      xhat <- mean(x)
    } else {
      xhat <- ifelse(mean(x) == 1, 1-1/horizon,
                     ifelse(mean(x) == 0, 1/horizon, mean(x)))
    }
    ni <- length(x)
    KLhat <- kl_ber(xhat, tau) #- 1/current_round
    return(ni * ifelse(is.na(KLhat),0,KLhat))
  }
  lapply(armls, get_metric)
}

#############################################################

# KL-UCB Algorithm

get_complexity <- function(means, tau, epsilon) {
  res <- vector(length = length(means))
  for(i in 1:length(means)) {
    res[i] <- max(kl_ber(means[i], tau), epsilon^2/2)
  }
  sum(1/res)
}

KLUCB_bandit_from_tsdata <- function(data, rounds = NA, 
                                     tau, epsilon,
                                     verbose = FALSE, seed = NA,
                                     horizon = NA, H = NA) {
  
  if(!is.na(seed)) set.seed(seed)
  K <- dim(data)[2]
  exp_factor <- horizon/H
  
  # initialize by pulling each arm once
  arm_list <- diag(as.matrix(data[1:K,]))
  arm_list <- as.list(arm_list)
  
  # initialize mean storage and the counter
  mean_storage <- matrix(unlist(lapply(arm_list, mean)), nrow = 1)
  arm_sequence <- 1:K
  
  # initialize metric storage
  metric_storage <- list()
  
  for(i in (K+1):rounds) {
    if(verbose) message(paste("this is round", i))
    metric_storage[[i-K]] <- unlist(get_next_arm_klucb(arm_list, 
                                                      tau = tau,
                                                      epsilon = epsilon,
                                                      horizon = horizon,
                                                      exp_factor = exp_factor))
    
    next_arm <- get_min(-metric_storage[[i-K]])
    arm_sequence <- c(arm_sequence, next_arm)
    
    if(verbose) message("arm selecting done")
    if(verbose) message(next_arm)
    arm_list[[next_arm]] <- c(arm_list[[next_arm]], data[i, next_arm])
    
    mean_storage <- rbind(mean_storage, unlist(lapply(arm_list, mean)))
    
    if(verbose) message("arm pulling done")
    
  }
  return(list(means = unlist(lapply(arm_list, mean)),
              arm_list = arm_list,
              arm_sequence = arm_sequence,
              mean_storage = mean_storage#,
              #metric_storage = metric_storage
              ))
}

get_klub <- function(N_arm, mu_hat, exp_factor) {
  potential_vals <- seq(mu_hat, 0.99999, by = 0.0005)
  max(potential_vals[N_arm * kl_ber(mu_hat, potential_vals) 
                            <= exp_factor])
}

get_kllb <- function(N_arm, mu_hat, exp_factor) {
  potential_vals <- seq(0.00001, mu_hat, by = 0.0005)
  min(potential_vals[N_arm * kl_ber(mu_hat, potential_vals) 
                            <= exp_factor])
}

get_next_arm_klucb <- function(armls, tau, epsilon, horizon, exp_factor) {
  get_metric <- function(x) {
    if(is.na(horizon)) {
      xhat <- mean(x)
    } else {
      xhat <- ifelse(mean(x) == 1, 1-1/horizon,
                     ifelse(mean(x) == 0, 1/horizon, mean(x)))
    }
    ni <- length(x)
    ul_tau_diff <- ifelse(xhat >= tau,
                          tau - get_kllb(ni, xhat, exp_factor),
                          get_klub(ni, xhat, exp_factor) - tau)
    return(max(epsilon, ul_tau_diff))
  }
  lapply(armls, get_metric)
}

##############################################################

# Probability of Improvement

PI_bandit_from_tsdata <- function(data, rounds, tau, epsilon,
                                  alpha, beta, verbose = FALSE, seed = NA, 
                                  tadj = FALSE) {
  
  if(!is.na(seed)) set.seed(seed)
  K <- dim(data)[2]
  
  # initialize by pulling each arm once
  arm_list <- diag(as.matrix(data[1:K,]))
  arm_list <- as.list(arm_list)
  
  # initialize mean storage and the counter
  mean_storage <- matrix(unlist(lapply(arm_list, mean)), nrow = 1)
  arm_sequence <- 1:K
  
  for(i in (K+1):rounds) {
    if(verbose) message(paste("this is round", i))
    # NOTE THE MINUS SIGN; find the arm that MAXIMIZES the error probability
    next_arm <- get_min(-unlist(get_next_arm_PI(arm_list, tau = tau,
                                                epsilon = epsilon,
                                                alpha = alpha,
                                                beta = beta, tadj = tadj)))
    arm_sequence <- c(arm_sequence, next_arm)
    
    if(verbose) message("arm selecting done")
    arm_list[[next_arm]] <- c(arm_list[[next_arm]], data[i, next_arm])
    mean_storage <- rbind(mean_storage, unlist(lapply(arm_list,
                                                      get_posterior_mean,
                                                      alpha = alpha,
                                                      beta = beta)))
    
    if(verbose) message("arm pulling done")
    
  }
  return(list(means = unlist(lapply(arm_list, mean)),
              arm_list = arm_list,
              arm_sequence = arm_sequence,
              mean_storage = mean_storage))
}

get_next_arm_PI <- function(armls, tau, epsilon,
                            alpha, beta, tadj) {
  get_metric <- function(x) {
    alpha_prime <- sum(x)+alpha
    beta_prime <- length(x)+beta-sum(x)
    posterior_mean <- alpha_prime/(alpha_prime+beta_prime)
    # depending on current mean estimate, give probability of error
    if(tadj) {
      ifelse(posterior_mean >= tau,
             pbeta(tau-epsilon, alpha_prime, beta_prime)/length(x),
             (1-pbeta(tau+epsilon, alpha_prime, beta_prime))/length(x))
    } else {
      ifelse(posterior_mean >= tau,
             pbeta(tau-epsilon, alpha_prime, beta_prime),
             1-pbeta(tau+epsilon, alpha_prime, beta_prime))
    }
  }
  lapply(armls, get_metric)
}

get_posterior_mean <- function(x, alpha, beta) {
  alpha_prime <- sum(x)+alpha
  beta_prime <- length(x)+beta-sum(x)
  alpha_prime/(alpha_prime+beta_prime)
}

##############################################################

# Threshold Thompson Sampling

TTS_from_tsdata <- function(data, rounds = 5000, tau, epsilon,
                            alpha = 1, beta = 1, verbose = FALSE, seed = NA) {
  
  if(!is.na(seed)) set.seed(seed)
  K <- dim(data)[2]
  
  # initialize by pulling each arm once
  arm_list <- diag(as.matrix(data[1:K,]))
  arm_list <- as.list(arm_list)
  
  # initialize mean storage and the counter
  mean_storage <- matrix(unlist(lapply(arm_list, mean)), nrow = 1)
  arm_sequence <- 1:K
  
  for(i in (K+1):rounds) {
    if(verbose) message(paste("this is round", i))
    # We can get the posterior error probability in the same way as for PI
    # But then we don't pick the strict minimum but proportional to probs
    next_arm <- draw_arm_tts(unlist(get_next_arm_PI(arm_list, tau = tau,
                                                    epsilon = epsilon,
                                                    alpha = alpha,
                                                    beta = beta,
                                                    tadj = FALSE)))
    arm_sequence <- c(arm_sequence, next_arm)
    
    if(verbose) message("arm selecting done")
    arm_list[[next_arm]] <- c(arm_list[[next_arm]], data[i, next_arm])
    mean_storage <- rbind(mean_storage, unlist(lapply(arm_list,
                                                      get_posterior_mean,
                                                      alpha = alpha, 
                                                      beta = beta)))
    
    if(verbose) message("arm pulling done")
  }
  return(list(means = unlist(lapply(arm_list, mean)),
              arm_list = arm_list,
              arm_sequence = arm_sequence,
              mean_storage = mean_storage))
}

draw_arm_tts <- function(x) {
  x <- x/sum(x) # make error probabilities proportional
  sample(1:length(x), 1, prob = x)
}

##############################################################

# Bayes-UCB

BayesUCB_from_tsdata <- function(data, rounds = 5000, tau, epsilon,
                                 alpha = 1, beta = 1, 
                                 with_epsilon = TRUE, const = 5,
                                 rate, verbose = FALSE, seed = NA) {
  
  if(!is.na(seed)) set.seed(seed)
  K <- dim(data)[2]
  
  # initialize by pulling each arm once
  arm_list <- diag(as.matrix(data[1:K,]))
  arm_list <- as.list(arm_list)
  
  # initialize mean storage and the counter
  mean_storage <- matrix(unlist(lapply(arm_list, mean)), nrow = 1)
  arm_sequence <- 1:K
  
  for(i in (K+1):rounds) {
    if(verbose) message(paste("this is round", i))
    
    next_arm <- get_min(-unlist(lapply(arm_list, get_BayesUCB_metric, rate = rate, 
                                       tau = tau, epsilon = epsilon,
                                       alpha = alpha, beta = beta,
                                       rounds = rounds, current_round = i,
                                       with_epsilon = with_epsilon,
                                       const = const)))
    arm_sequence <- c(arm_sequence, next_arm)
    
    if(verbose) message("arm selecting done")
    arm_list[[next_arm]] <- c(arm_list[[next_arm]], data[i, next_arm])
    mean_storage <- rbind(mean_storage, unlist(lapply(arm_list,
                                                      get_posterior_mean,
                                                      alpha = alpha, 
                                                      beta = beta)))
    if(verbose) message("arm pulling done")
  }
  return(list(means = unlist(lapply(arm_list, mean)),
              arm_list = arm_list,
              arm_sequence = arm_sequence,
              mean_storage = mean_storage))
}

get_BayesUCB_metric <- function(x, rate, tau, epsilon, 
                                alpha, beta, rounds, current_round,
                                with_epsilon, const) {
  alpha_prime <- sum(x)+alpha
  beta_prime <- length(x)+beta-sum(x)
  posterior_mean <- alpha_prime/(alpha_prime+beta_prime)
  #ni <- length(x)
  # depending on current mean estimate, give probability of error
  
  if(rate == "inverse") {
    upper_quantile <- qbeta(1-1/current_round, alpha_prime, beta_prime)
    lower_quantile <- qbeta(1/current_round, alpha_prime, beta_prime)
  }
  #if(rate == "klucb") {
  #  upper_quantile <- qbeta(1-c/ni, alpha_prime, beta_prime)
  #  lower_quantile <- qbeta(c/ni, alpha_prime, beta_prime)
  #}
  if(rate == "inverse_horizon_linear") {
    upper_quantile <- qbeta(1-1/(rounds), alpha_prime, beta_prime)
    lower_quantile <- qbeta(1/(rounds), alpha_prime, beta_prime)
  }
  if(rate == "inverse_horizon_linear_c") {
    upper_quantile <- qbeta(1-1/(const*rounds), alpha_prime, beta_prime)
    lower_quantile <- qbeta(1/(const*rounds), alpha_prime, beta_prime)
  }
  if(rate == "inverse_horizon") {
    upper_quantile <- qbeta(1-1/(current_round*log(rounds)), alpha_prime, beta_prime)
    lower_quantile <- qbeta(1/(current_round*log(rounds)), alpha_prime, beta_prime)
  }
  if(rate == "inverse_horizon_c") {
    upper_quantile <- qbeta(1-1/(current_round*log(rounds)^const), alpha_prime, beta_prime)
    lower_quantile <- qbeta(1/(current_round*log(rounds)^const), alpha_prime, beta_prime)
  }
  if(rate == "inverse_squared") {
    upper_quantile <- qbeta(1-1/current_round^2, alpha_prime, beta_prime)
    lower_quantile <- qbeta(1/current_round^2, alpha_prime, beta_prime)
  }
  if(rate == "inverse_squared_horizon_c") {
    upper_quantile <- qbeta(1-1/(current_round^2*log(rounds)^const), alpha_prime, beta_prime)
    lower_quantile <- qbeta(1/(current_round^2*log(rounds)^const), alpha_prime, beta_prime)
  }
  if(rate == "inverse_cubic") {
    upper_quantile <- qbeta(1-1/current_round^3, alpha_prime, beta_prime)
    lower_quantile <- qbeta(1/current_round^3, alpha_prime, beta_prime)
  }
  if(rate == "inverse_power5") {
    upper_quantile <- qbeta(1-1/current_round^5, alpha_prime, beta_prime)
    lower_quantile <- qbeta(1/current_round^5, alpha_prime, beta_prime)
  }
  if(rate == "inverse_log") {
    upper_quantile <- qbeta(1-1/log(current_round), alpha_prime, beta_prime)
    lower_quantile <- qbeta(1/log(current_round), alpha_prime, beta_prime)
  }
  if(rate == "inverse_sqrt") {
    upper_quantile <- qbeta(1-1/sqrt(current_round), alpha_prime, beta_prime)
    lower_quantile <- qbeta(1/sqrt(current_round), alpha_prime, beta_prime)
  }
  
  if(with_epsilon) {
    ifelse(posterior_mean >= tau,
           (tau - epsilon - lower_quantile),
           (upper_quantile - tau - epsilon))
  } else { # this alternative doesn't change anything because we always just substract epsilon from everything...
    ifelse(posterior_mean >= tau,
           (tau - lower_quantile),
           (upper_quantile - tau))
  }
}


##############################################################

uniform_bandit_from_tsdata <- function(data, rounds = 5000,
                                       verbose = FALSE, seed = NA) {
  if(!is.na(seed)) set.seed(seed)
  K <- dim(data)[2]
  
  # initialize by pulling each arm once
  arm_list <- diag(as.matrix(data[1:K,]))
  arm_list <- as.list(arm_list)
  
  # initialize mean storage and the counter
  mean_storage <- matrix(unlist(lapply(arm_list, mean)), nrow = 1)
  arm_sequence <- 1:K
  
  for(i in (K+1):rounds) {
    if(verbose) message(paste("this is round", i))
    next_arm <- get_min(unlist(get_next_arm_uniform(arm_list)))
    arm_sequence <- c(arm_sequence, next_arm)
    
    if(verbose) message("arm selecting done")
    arm_list[[next_arm]] <- c(arm_list[[next_arm]], data[i, next_arm])
    mean_storage <- rbind(mean_storage, unlist(lapply(arm_list, mean)))
    
    if(verbose) message("arm pulling done")
    
  }
  return(list(means = unlist(lapply(arm_list, mean)),
              arm_list = arm_list,
              arm_sequence = arm_sequence,
              mean_storage = mean_storage))
}

get_next_arm_uniform <- function(armls) {
  lapply(armls, length)
}

##############################################################

# Likelihood Ratio based Algorithm
# for Bernoulli distribution

LR_bandit_from_tsdata <- function(data, rounds = 5000, tau, epsilon,
                                  verbose = FALSE, seed = NA) {
  
  if(!is.na(seed)) set.seed(seed)
  K <- dim(data)[2]
  
  # initialize by pulling each arm once
  arm_list <- diag(as.matrix(data[1:K,]))
  arm_list <- as.list(arm_list)
  
  # initialize mean storage and the counter
  mean_storage <- matrix(unlist(lapply(arm_list, mean)), nrow = 1)
  arm_sequence <- 1:K
  
  for(i in (K+1):rounds) {
    if(verbose) message(paste("this is round", i))
    
    next_arm <- get_min(-unlist(get_next_arm_lr(arm_list, tau = tau,
                                                epsilon = epsilon)))
    
    arm_sequence <- c(arm_sequence, next_arm)
    
    if(verbose) message("arm selecting done")
    if(verbose) message(next_arm)
    arm_list[[next_arm]] <- c(arm_list[[next_arm]], data[i, next_arm])
    
    mean_storage <- rbind(mean_storage, unlist(lapply(arm_list, mean)))
    
    if(verbose) message("arm pulling done")
    
  }
  return(list(means = unlist(lapply(arm_list, mean)),
              arm_list = arm_list,
              arm_sequence = arm_sequence,
              mean_storage = mean_storage))
}

lr_ber <- function(x,S,N) {
  # here, S/N is the empirical mean (successes over trials)
  # while x represents some H0 mean we test against
  (x/(S/N))^(S) * ((1-x)/(1-S/N))^(N-S)
}

get_next_arm_lr <- function(armls, tau, epsilon) {
  get_metric <- function(x) {
    successes <- sum(x)
    trials <- length(x)
    LRhat <- ifelse(successes/trials >= tau,
                    lr_ber(tau-epsilon, successes, trials),
                    lr_ber(tau+epsilon, successes, trials))
    return(LRhat)
  }
  lapply(armls, get_metric)
}

##############################################################

# Likelihood Ratio based Algorithm
# for Bernoulli distribution
# with additional D-tracking rule from
# Garivier, Kaufmann (2016): Optimal Best Arm Identification with Fixed Confidence

LRD_bandit_from_tsdata <- function(data, rounds = 5000, tau, epsilon,
                                   verbose = FALSE, seed = NA) {
  
  if(!is.na(seed)) set.seed(seed)
  K <- dim(data)[2]
  
  # initialize by pulling each arm once
  arm_list <- diag(as.matrix(data[1:K,]))
  arm_list <- as.list(arm_list)
  
  # initialize mean storage and the counter
  mean_storage <- matrix(unlist(lapply(arm_list, mean)), nrow = 1)
  arm_sequence <- 1:K
  
  for(i in (K+1):rounds) {
    if(verbose) message(paste("this is round", i))
    
    # D-Tracking Rule
    Uset <- which(table(arm_sequence) < sqrt(i) - K/2)
    if(length(Uset) > 0){
      next_arm <- get_min(table(arm_sequence))
    } else {
      next_arm <- get_min(-unlist(get_next_arm_lr(arm_list, tau = tau,
                                                  epsilon = epsilon)))
    }
    
    arm_sequence <- c(arm_sequence, next_arm)
    
    if(verbose) message("arm selecting done")
    if(verbose) message(next_arm)
    arm_list[[next_arm]] <- c(arm_list[[next_arm]], data[i, next_arm])
    
    mean_storage <- rbind(mean_storage, unlist(lapply(arm_list, mean)))
    
    if(verbose) message("arm pulling done")
    
  }
  return(list(means = unlist(lapply(arm_list, mean)),
              arm_list = arm_list,
              arm_sequence = arm_sequence,
              mean_storage = mean_storage))
}

##############################################################

# for Normal distribution
# actually uses log(-LikelihoodRatio) as the test statistic since
# the LR goes to 0 too quickly
# -> kl_gaussian

LR_bandit_from_tsdata_gaussian <- function(data, rounds = 5000, tau, epsilon,
                                  verbose = FALSE, seed = NA) {
  
  if(!is.na(seed)) set.seed(seed)
  K <- dim(data)[2]
  
  # initialize by pulling each arm once
  arm_list <- diag(as.matrix(data[1:K,]))
  arm_list <- as.list(arm_list)
  
  # initialize mean storage and the counter
  mean_storage <- matrix(unlist(lapply(arm_list, mean)), nrow = 1)
  arm_sequence <- 1:K
  
  for(i in (K+1):rounds) {
    if(verbose) message(paste("this is round", i))
    
    next_arm <- get_min(unlist(get_next_arm_kl_gaussian(arm_list, tau = tau,
                                                         epsilon = epsilon)))
    
    arm_sequence <- c(arm_sequence, next_arm)
    
    if(verbose) message("arm selecting done")
    if(verbose) message(next_arm)
    arm_list[[next_arm]] <- c(arm_list[[next_arm]], data[i, next_arm])
    
    mean_storage <- rbind(mean_storage, unlist(lapply(arm_list, mean)))
    
    if(verbose) message("arm pulling done")
    
  }
  return(list(means = unlist(lapply(arm_list, mean)),
              arm_list = arm_list,
              arm_sequence = arm_sequence,
              mean_storage = mean_storage))
}

lr_gaussian <- function(x, tau, mu, n) {
  sigma_tau <- 1/n*sum((x-tau)^2)
  sigma_emp <- 1/n*sum((x-mu)^2)
  prod(dnorm(x, tau, sqrt(sigma_tau))) / prod(dnorm(x, mu, sqrt(sigma_emp)))
}

#kl_gaussian <- function(x, tau, mu, n) {
#  sigma_tau <- 1/n*sum((x-tau)^2)
#  sigma_emp <- 1/n*sum((x-mu)^2)
#  ifelse(sigma_emp == 0, -Inf, n*log(sqrt(sigma_tau)/sqrt(sigma_emp)))
#}

kl_gaussian <- function(x, tau, mu, n) {
  sigma_tau <- 1/n*sum((x-tau)^2)
  sigma_emp <- 1/n*sum((x-mu)^2)
  ifelse(sigma_emp == 0, 1/10000, n*log(sqrt(sigma_tau)/sqrt(sigma_emp)))
}

get_next_arm_lr_gaussian <- function(armls, tau, epsilon) {
  get_metric <- function(x) {
    successes <- sum(x)
    trials <- length(x)
    LRhat <- ifelse(successes/trials >= tau,
                    lr_gaussian(x, tau-epsilon, successes/trials, trials),
                    lr_gaussian(x, tau+epsilon, successes/trials, trials))
    return(LRhat)
  }
  lapply(armls, get_metric)
}

get_next_arm_kl_gaussian <- function(armls, tau, epsilon) {
  get_metric <- function(x) {
    LRhat <- ifelse(mean(x) >= tau,
                    kl_gaussian(x, tau-epsilon, mean(x), length(x)),
                    kl_gaussian(x, tau+epsilon, mean(x), length(x)))
    return(LRhat)
  }
  lapply(armls, get_metric)
}

##############################################################

LR_bandit_from_tsdata_exponential <- function(data, rounds = 5000, tau, epsilon,
                                      verbose = FALSE, seed = NA) {
  
  if(!is.na(seed)) set.seed(seed)
  K <- dim(data)[2]
  
  # initialize by pulling each arm once
  arm_list <- diag(as.matrix(data[1:K,]))
  arm_list <- as.list(arm_list)
  
  # initialize mean storage and the counter
  mean_storage <- matrix(unlist(lapply(arm_list, mean)), nrow = 1)
  arm_sequence <- 1:K
  
  for(i in (K+1):rounds) {
    if(verbose) message(paste("this is round", i))
    
    next_arm <- get_min(unlist(get_next_arm_kl_exponential(arm_list, tau = tau,
                                                           epsilon = epsilon)))
    
    arm_sequence <- c(arm_sequence, next_arm)
    
    if(verbose) message("arm selecting done")
    if(verbose) message(next_arm)
    arm_list[[next_arm]] <- c(arm_list[[next_arm]], data[i, next_arm])
    
    mean_storage <- rbind(mean_storage, unlist(lapply(arm_list, mean)))
    
    if(verbose) message("arm pulling done")
    
  }
  return(list(means = unlist(lapply(arm_list, mean)),
              arm_list = arm_list,
              arm_sequence = arm_sequence,
              mean_storage = mean_storage))
}

lr_exponential <- function(x, mu, tau) {
  prod(dexp(x, 1/tau))/prod(dexp(x, 1/mu))
}

kl_exponential <- function(n, mu, tau) {
  n*(log(tau/mu)+mu/tau-1)
}

get_next_arm_kl_exponential <- function(armls, tau, epsilon) {
  get_metric <- function(x) {
    LRhat <- ifelse(mean(x) >= tau,
                    kl_exponential(length(x), mean(x), tau-epsilon),
                    kl_exponential(length(x), mean(x), tau+epsilon))
    return(LRhat)
  }
  lapply(armls, get_metric)
}


##############################################################

LR_bandit_from_tsdata_poisson <- function(data, rounds = 5000, tau, epsilon,
                                          verbose = FALSE, seed = NA) {
  
  if(!is.na(seed)) set.seed(seed)
  K <- dim(data)[2]
  
  # initialize by pulling each arm once
  arm_list <- diag(as.matrix(data[1:K,]))
  arm_list <- as.list(arm_list)
  
  # initialize mean storage and the counter
  mean_storage <- matrix(unlist(lapply(arm_list, mean)), nrow = 1)
  arm_sequence <- 1:K
  
  for(i in (K+1):rounds) {
    if(verbose) message(paste("this is round", i))
    
    next_arm <- get_min(unlist(get_next_arm_kl_poisson(arm_list, tau = tau,
                                                       epsilon = epsilon)))
    
    arm_sequence <- c(arm_sequence, next_arm)
    
    if(verbose) message("arm selecting done")
    if(verbose) message(next_arm)
    arm_list[[next_arm]] <- c(arm_list[[next_arm]], data[i, next_arm])
    
    mean_storage <- rbind(mean_storage, unlist(lapply(arm_list, mean)))
    
    if(verbose) message("arm pulling done")
    
  }
  return(list(means = unlist(lapply(arm_list, mean)),
              arm_list = arm_list,
              arm_sequence = arm_sequence,
              mean_storage = mean_storage))
}

lr_poisson <- function(x, mu, tau) {
  prod(dpois(x, tau))/prod(dpois(x, mu))
}

kl_poisson <- function(n, mu, tau) {
  mu_noise <- ifelse(mu == 0, 0.0001/n, mu)
  n * (tau - mu + log(mu_noise/tau)*mu)
}

get_next_arm_kl_poisson <- function(armls, tau, epsilon) {
  get_metric <- function(x) {
    LRhat <- ifelse(mean(x) >= tau,
                    kl_poisson(length(x), mean(x), tau-epsilon),
                    kl_poisson(length(x), mean(x), tau+epsilon))
    return(LRhat)
  }
  lapply(armls, get_metric)
}

##############################################################

# EVT Algorithm

get_next_arm_evt <- function(armls, tau, rounds, epsilon) {
  get_metric <- function(x) {
    xhat <- mean(x)
    ni <- length(x)
    varhat <- (ni-1)/ni*var(x)
    return(sqrt(ni) * (sqrt(varhat+(abs(xhat-tau)+epsilon))-sqrt(varhat)))
  }
  lapply(armls, get_metric)
}

EVT_from_tsdata <- function(data, rounds = 5000, tau, epsilon,
                            verbose = FALSE, seed = NA) {
  
  if(!is.na(seed)) set.seed(seed)
  K <- dim(data)[2]
  
  # initialize by pulling each arm TWICE! to get a variance if not bernoulli
  arm_list_init <- rbind(diag(as.matrix(data[1:K,])), 
                         diag(as.matrix(data[(K+1):(2*K),])))
  arm_list <- list()
  for(i in 1:K) {
    arm_list[[i]] <- arm_list_init[,i]
  }
  
  # initialize mean storage, arm sequence, and the counter
  mean_storage <- matrix(unlist(lapply(arm_list, mean)), nrow = 1)
  arm_sequence <- c(1:K,1:K)
  
  for(i in (2*K+1):rounds) {
    if(verbose) message(paste("this is round", i))
    next_arm <- get_min(unlist(get_next_arm_evt(arm_list, tau = tau, 
                                                rounds = rounds, 
                                                epsilon = epsilon)))
    arm_sequence <- c(arm_sequence, next_arm)
    
    if(verbose) message("arm selecting done")
    arm_list[[next_arm]] <- c(arm_list[[next_arm]], data[i, next_arm])
    mean_storage <- rbind(mean_storage, unlist(lapply(arm_list, mean)))
    
    if(verbose) message("arm pulling done")
    
  }
  return(list(means = unlist(lapply(arm_list, mean)),
              arm_list = arm_list,
              input = data,
              arm_sequence = arm_sequence,
              mean_storage = mean_storage))
}

##############################################################


#mean_002 <- c(0.033, 0.033, 0.037, 0.037,
#              0.031, 0.031, 0.039, 0.039,
#              0.027, 0.027, 0.043, 0.043,
#              0.02, 0.02, 0.05, 0.05,
#              0.005, 0.005, 0.065, 0.065)
#tau_002 <- 0.035
#epsilon_002 <- 0.005

#set.seed(512)
#testt <- data.frame(rep(NA, times = 6000))
#for(i in 1:length(mean_002)) {
#  testt[[i]] <- as.numeric(purrr::rbernoulli(6000, p  = mean_002[i]))
#}

#apt_test <- APT_from_tsdata(testt, rounds = 5000, tau = tau_002, 
#                            epsilon = epsilon_002, seed = 46)
#uniform_test <- uniform_bandit_from_tsdata(testt, rounds = 5000, seed = 46)
#bucb_test <- BayesUCB_from_tsdata(testt, rounds = 5000, rate = "inverse",
#                                  tau = tau_002, epsilon = epsilon_002, 
#                                  alpha = tau_002, beta = 1-tau_002, seed = 46)
#tts_test <- TTS_from_tsdata(testt, rounds = 5000, 
#                            tau = tau_002, epsilon = epsilon_002, 
#                            alpha = tau_002, beta = 1 - tau_002, seed = 46)
#pi_test <- PI_bandit_from_tsdata(testt, rounds = 5000, 
#                                 tau = tau_002, epsilon = epsilon_002, 
#                                 alpha = tau_002, beta = 1 - tau_002, seed = 46)
#kl_test <- KL_bandit_from_tsdata(testt, rounds = 5000, tau = tau_002, 
#                                 epsilon = epsilon_002, seed = 46,
#                                 at_tau = FALSE)


#table(apt_test$arm_sequence)
#plot(apt_test$arm_sequence, type = "l")
#table(uniform_test$arm_sequence)
#plot(uniform_test$arm_sequence, type = "l")
#table(bucb_test$arm_sequence)
#plot(bucb_test$arm_sequence, type = "l")
#table(tts_test$arm_sequence)
#plot(tts_test$arm_sequence, type = "l")
#table(pi_test$arm_sequence)
#plot(pi_test$arm_sequence, type = "l")
#table(kl_test$arm_sequence)
#plot(kl_test$arm_sequence, type = "l")


