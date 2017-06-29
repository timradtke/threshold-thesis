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

# Probability of Improvement

BETA_from_tsdata <- function(data, rounds, tau, epsilon,
                             alpha, beta, verbose = FALSE, 
                             seed = NA) {
  
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
    next_arm <- get_min(-unlist(get_next_arm_BETA(arm_list, tau = tau,
                                                  epsilon = epsilon,
                                                  alpha = alpha,
                                                  beta = beta)))
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

get_next_arm_BETA <- function(armls, tau, epsilon,
                              alpha, beta) {
  get_metric <- function(x) {
    alpha_prime <- sum(x)+alpha
    beta_prime <- length(x)+beta-sum(x)
    posterior_mean <- alpha_prime/(alpha_prime+beta_prime)
    border <- ifelse(posterior_mean >= tau, tau-epsilon, tau+epsilon)
    alpha_alt <- border * (length(x)+alpha+beta)
    beta_alt <- (1-border) * (length(x)+alpha+beta)
    act_vec <- rbeta(100000, alpha_prime, beta_prime)
    alt_vec <- rbeta(100000, alpha_alt, beta_alt)
    if(posterior_mean >= tau) {
      mean(act_vec < alt_vec)
    } else {
      mean(act_vec > alt_vec)
    }
  }
  lapply(armls, get_metric)
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
                                                    beta = beta)))
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

# Bayes UCB

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
  # depending on current mean estimate, give probability of error
  
  if(rate == "inverse") {
    upper_quantile <- qbeta(1-1/current_round, alpha_prime, beta_prime)
    lower_quantile <- qbeta(1/current_round, alpha_prime, beta_prime)
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
  if(rate == "inverse_log") {
    upper_quantile <- qbeta(1-1/log(current_round), alpha_prime, beta_prime)
    lower_quantile <- qbeta(1/log(current_round), alpha_prime, beta_prime)
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


