# Using Bayes Quantiles, similar to
# On Bayesian Upper Confidence Bounds for Bandit Problems
# Kaufmann, Cappe, Garivier

#system.time(res_bucb1 <- BayesUCB(means = mean_twenty,
#                                   tau = 0.04, epsilon = 0.01,
#                                   alpha = 0.04, beta = 0.96,
#                                   K = 20, rounds = 2000,
#                                   variances = NA, seed = 60,
#                                   rate = "inverse"))
#system.time(res_bucb2 <- BayesUCB(means = mean_twenty,
#                                    tau = 0.04, epsilon = 0.01,
#                                    alpha = 0.04, beta = 0.96,
#                                    K = 20, rounds = 2000,
#                                    variances = NA, seed = 60,
#                                    rate = "inverse_horizon"))
#system.time(res_bucb3 <- BayesUCB(means = mean_twenty,
#                                    tau = 0.04, epsilon = 0.01,
#                                    alpha = 0.04, beta = 0.96,
#                                    K = 20, rounds = 2000,
#                                    variances = NA, seed = 60,
#                                    rate = "inverse_squared"))
#system.time(res_bucb4 <- BayesUCB(means = mean_twenty,
#                                    tau = 0.04, epsilon = 0.01,
#                                    alpha = 0.04, beta = 0.96,
#                                    K = 20, rounds = 2000,
#                                    variances = NA, seed = 60,
#                                    rate = "inverse_log"))
#lapply(res_bucb1$arm_list, length)
#data.frame(means = mean_twenty,
#           a = unlist(lapply(res_bucb1$arm_list, length)),
#           b = unlist(lapply(res_bucb2$arm_list, length)),
#           c = unlist(lapply(res_bucb3$arm_list, length)),
#           d = unlist(lapply(res_bucb4$arm_list, length)))
#lapply(res_bucb2$arm_list, length)
#res_bucb$means#

#res_bucb <- as.data.frame(res_bucb$mean_storage)
#res_bucb$round <- 1:dim(res_bucb)[1]
#res_bucb_long <- gather(res_bucb[1:700,], key = "Arm", value = "Mean", 
#                        -round)
#ggplot(res_bucb_long, aes(x = round, y = Mean, color = Arm)) +
#  geom_line()

###############################################################################

BayesUCB <- function(means, K = 4, rounds = 5000, tau = 0.5, epsilon = 0.1,
                     alpha = 1, beta = 1, mean_prior = NA, n_prior = NA,
                     rate,
                     variances = NA, verbose = FALSE, seed = NA) {
  
  if(!is.na(seed)) set.seed(seed)
  
  # initialize by pulling each arm once
  arm_list <- NA
  for(i in 1:K) {
    # Ideally it should also take an argument specifying how
    # the sample is generated (e.g., from existing data, bernoulli or other distr.)
    arm_list[[i]] <- pull_arm(k = i, means = means, variances = variances)
  }
  arm_list <- as.list(arm_list)
  # initialize mean storage and the counter
  mean_storage <- matrix(unlist(lapply(arm_list, mean)), nrow = 1)
  counter <- K
  
  for(i in (K+1):rounds) {
    if(verbose) message(paste("this is round", i))
    
    next_arm <- get_min(-unlist(lapply(arm_list, get_BayesUCB_metric, rate = rate, 
                                       tau = tau, epsilon = epsilon,
                                       alpha = alpha, beta = beta,
                                       rounds = rounds, current_round = i)))
    #next_arm <- get_min(-unlist(get_next_arm_BayesUCB(arm_list, tau = tau,
    #                                                epsilon = epsilon,
    #                                                alpha = alpha,
    #                                                beta = beta,
    #                                                rate = rate,
    #                                                rounds = rounds,
    #                                                current_round = i)))
    
    if(verbose) message("arm selecting done")
    arm_list[[next_arm]] <- c(arm_list[[next_arm]], pull_arm(k = next_arm, 
                                                             means = means,
                                                             variances = variances))
    mean_storage <- rbind(mean_storage, unlist(lapply(arm_list,
                                                      get_posterior_mean,
                                                      alpha = alpha, beta = beta,
                                                      mean_prior = mean_prior,
                                                      n_prior = n_prior)))
    counter <- counter + 1
    
    if(verbose) message("arm pulling done")
    
  }
  return(list(means = unlist(lapply(arm_list, mean)),
              arm_list = arm_list,
              true_means = means,
              counter = counter,
              mean_storage = mean_storage))
}

# Ideally it should also take an argument specifying how
# the sample is generated (e.g., from existing data, bernoulli or other distr.)
pull_arm <- function(k, means, variances = NA) {
  if(is.na(variances)) {
    rbinom(1,1,means[k])
  } else {
    rnorm(1,means[k], variances[k])
  }
}

get_BayesUCB_metric <- function(x, rate, tau, epsilon, alpha, beta, rounds, current_round) {
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
  if(rate == "inverse_squared") {
    upper_quantile <- qbeta(1-1/current_round^2, alpha_prime, beta_prime)
    lower_quantile <- qbeta(1/current_round^2, alpha_prime, beta_prime)
  }
  if(rate == "inverse_log") {
    upper_quantile <- qbeta(1-1/log(current_round), alpha_prime, beta_prime)
    lower_quantile <- qbeta(1/log(current_round), alpha_prime, beta_prime)
  }
  
  ifelse(posterior_mean >= tau,
         tau - epsilon - lower_quantile,
         upper_quantile - tau - epsilon
  )
}



get_next_arm_BayesUCB <- function(armls, tau, epsilon, current_round, rounds,
                                  alpha, beta, rate) {
  # Get how far away the quantile is from the respective tau+epsilon
  
  
  lapply(armls, get_metric, rate = rate)
}

get_min <- function(x) {
  mini <- min(x)
  ifelse(sum(x==mini) == 1, which.min(x),
         sample(which(x==mini)))
}

get_posterior_mean <- function(x, alpha, beta, mean_prior, n_prior) {
  # depending on whether variances are available,
  # do method for bernoulli or gaussian
  if(is.na(mean_prior)) {
    alpha_prime <- sum(x)+alpha
    beta_prime <- length(x)+beta-sum(x)
    alpha_prime/(alpha_prime+beta_prime)
  } else {
    # Gaussian
    n <- length(x)
    tau_posterior <- get_mean_from_gamma(x, alpha_prior = alpha,
                                         beta_prior = beta, 
                                         mean_prior = mean_prior, 
                                         n_prior = n_prior)
    mean_posterior <- n * tau_posterior / (n*tau_posterior + n_prior*tau_posterior) * mean(x) +
      n_prior * tau_posterior / (n*tau_posterior + n_prior*tau_posterior) * mean_prior
    return(mean_posterior)
  }
}

