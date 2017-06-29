# Using Bayes Quantiles, similar to
# On Bayesian Upper Confidence Bounds for Bandit Problems
# Kaufmann, Cappe, Garivier

#means <- c(0.55, 0.38, 0.95)
#variances <- c(0.5, 0.9, 1.7)
#system.time(res_tts <- TTS(means = means, alpha = 1, beta = 1,
#                           K = 3, rounds = 5000, 
#                           variances = NA, seed = 54))
# user  system elapsed 
#0.390   0.057   0.446 
#lapply(res_tts$arm_list, length)
#res_tts$means
#means
#plot(x = c(0,5000), y = c(0,2), type = "n")
#lines(res_tts$mean_storage[,1])
#lines(res_tts$mean_storage[,2], col = "blue")
#lines(res_tts$mean_storage[,3], col = "red")

#system.time(res_blucb <- BayesLUCB(means = mean_twenty,
#                                   tau = 0.04, epsilon = 0.01,
#                                   alpha = 0.04, beta = 0.96,
#                                   K = 20, rounds = 2000,
#                                   variances = NA, seed = 60))
##lapply(res_bucb$arm_list, length)
#lapply(res_bucb2$arm_list, length)
#res_bucb$means#

#res_bucb <- as.data.frame(res_bucb$mean_storage)
#res_bucb$round <- 1:dim(res_bucb)[1]
#res_bucb_long <- gather(res_bucb[1:700,], key = "Arm", value = "Mean", 
#                        -round)
#ggplot(res_bucb_long, aes(x = round, y = Mean, color = Arm)) +
#  geom_line()

###############################################################################

BayesLUCB <- function(means, K = 4, rounds = 5000, tau = 0.5, epsilon = 0.1,
                      alpha = 1, beta = 1, mean_prior = NA, n_prior = NA,
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
  
  for(i in (K+1):(ceiling(rounds/2))) {
    if(verbose) message(paste("this is round", i))
    # We can get the posterior error probability in the same way as for PI
    # But then we don't pick the strict minimum but proportional to probs
    next_lower_arm <- get_min(-unlist(get_next_arm_lower(arm_list, 
                                                         tau = tau,
                                                          epsilon = epsilon,
                                                          alpha = alpha,
                                                          beta = beta,
                                                          current_round = (2*i)-1 )))
    next_upper_arm <- get_min(-unlist(get_next_arm_upper(arm_list, 
                                                         tau = tau,
                                                         epsilon = epsilon,
                                                         alpha = alpha,
                                                         beta = beta,
                                                         current_round = (2*i)-1 )))
    
    if(verbose) message("arm selecting done")
    arm_list[[next_lower_arm]] <- c(arm_list[[next_lower_arm]], pull_arm(k = next_lower_arm, 
                                                             means = means,
                                                             variances = variances))
    arm_list[[next_upper_arm]] <- c(arm_list[[next_upper_arm]], pull_arm(k = next_upper_arm, 
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


get_next_arm_lower <- function(armls, tau, epsilon, current_round,
                               alpha, beta) {
  # Get how far away the quantile is from the respective tau+epsilon
  get_metric <- function(x) {
    alpha_prime <- sum(x)+alpha
    beta_prime <- length(x)+beta-sum(x)
    posterior_mean <- alpha_prime/(alpha_prime+beta_prime)
    # depending on current mean estimate, give probability of error
    upper_quantile <- qbeta(1-1/current_round, alpha_prime, beta_prime)
    lower_quantile <- qbeta(1/current_round, alpha_prime, beta_prime)
    ifelse(posterior_mean >= tau,
           -999999,
           upper_quantile - tau - epsilon
    )
  }
  lapply(armls, get_metric)
}

get_next_arm_upper <- function(armls, tau, epsilon, current_round,
                               alpha, beta) {
  # Get how far away the quantile is from the respective tau+epsilon
  get_metric <- function(x) {
    alpha_prime <- sum(x)+alpha
    beta_prime <- length(x)+beta-sum(x)
    posterior_mean <- alpha_prime/(alpha_prime+beta_prime)
    # depending on current mean estimate, give probability of error
    upper_quantile <- qbeta(1-1/current_round, alpha_prime, beta_prime)
    lower_quantile <- qbeta(1/current_round, alpha_prime, beta_prime)
    ifelse(posterior_mean >= tau,
           tau - epsilon - lower_quantile,
           -999999
    )
  }
  lapply(armls, get_metric)
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

