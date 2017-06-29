# Use probability of error but do proportional sampling
# Threshold Thompson Sampling-like

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

#system.time(res_tts <- TTS(means = mean_twenty,
#                          tau = 0.04, epsilon = 0.01,
#                          alpha = 0.04, beta = 0.96,
#                          K = 20, rounds = 2000,
#                          variances = NA, seed = 60))
#system.time(res_tts2 <- TTS(means = mean_twenty,
#                           tau = 0.04, epsilon = 0.01,
#                           alpha = 0.04, beta = 0.96,
#                           K = 20, rounds = 2000,
#                           variances = NA, seed = 60,
#                           interval_adj = TRUE))
#lapply(res_tts$arm_list, length)
#res_tts$means

#tts_means <- as.data.frame(res_tts$mean_storage)
#tts_means$round <- 1:dim(tts_means)[1]
#tts_means_long <- gather(tts_means[1:700,], key = "Arm", value = "Mean", -round)
#ggplot(tts_means_long, aes(x = round, y = Mean, color = Arm)) +
#  geom_line()

###############################################################################

TTS <- function(means, K = 4, rounds = 5000, tau = 0.5, epsilon = 0.1,
                alpha = 1, beta = 1, mean_prior = NA, n_prior = NA,
                variances = NA, verbose = FALSE, seed = NA,
                interval_adj = FALSE) {
  
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
    # We can get the posterior error probability in the same way as for PI
    # But then we don't pick the strict minimum but proportional to probs
    next_arm <- draw_arm_tts(unlist(get_next_arm_PI(arm_list, tau = tau,
                                                    epsilon = epsilon,
                                                    alpha = alpha,
                                                    beta = beta,
                                                    mean_prior = mean_prior,
                                                    n_prior = n_prior,
                                                    interval_adj = interval_adj)))
    
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

get_mean_from_gamma <- function(x, alpha_prior, beta_prior, mean_prior, n_prior) {
  n <- length(x)
  alpha_posterior <- alpha_prior + n/2
  beta_posterior <- beta_prior + 1/2*(sum((x-mean(x))^2)) + 
    n * n_prior / (2*(n+n_prior)) * (mean(x) - mean_prior)^2
  
  alpha_posterior/beta_posterior
}

get_mean_from_normal_posterior <- function(x, mean_prior, n_prior, tau_posterior) {
  n <- length(x)
  n * tau_posterior / (n*tau_posterior + n_prior*tau_posterior) * mean(x) +
    n_prior * tau_posterior / (n*tau_posterior + n_prior*tau_posterior) * mean_prior
}

get_next_arm_PI <- function(armls, tau, epsilon,
                            alpha, beta, mean_prior, n_prior,
                            interval_adj) {
  # depending on whether variances are available, do method for bernoulli or gaussian
  if(is.na(mean_prior)) {
    # Bernoulli
    get_metric <- function(x) {
      alpha_prime <- sum(x)+alpha
      beta_prime <- length(x)+beta-sum(x)
      posterior_mean <- alpha_prime/(alpha_prime+beta_prime)
      # depending on current mean estimate, give probability of error
      if(interval_adj == FALSE) {
        ifelse(posterior_mean >= tau,
               pbeta(tau-epsilon, alpha_prime, beta_prime),
               1-pbeta(tau+epsilon, alpha_prime, beta_prime))
      } else {
        interval_prob <- pbeta(tau+epsilon, alpha_prime, beta_prime) - pbeta(tau-epsilon, alpha_prime, beta_prime)
        ifelse(posterior_mean >= tau,
               pbeta(tau-epsilon, alpha_prime, beta_prime)/interval_prob,
               (1-pbeta(tau+epsilon, alpha_prime, beta_prime))/interval_prob)
      }
      
    }
  } else {
    # Gaussian
    get_metric <- function(x) {
      n <- length(x)
      tau_posterior <- get_mean_from_gamma(x, alpha_prior = alpha,
                                           beta_prior = beta, 
                                           mean_prior = mean_prior, 
                                           n_prior = n_prior)
      mean_posterior <- n * tau_posterior / (n*tau_posterior + n_prior*tau_posterior) * mean(x) +
        n_prior * tau_posterior / (n*tau_posterior + n_prior*tau_posterior) * mean_prior
      # depending on current mean estimate, give probability of error
      ifelse(mean_posterior >= tau,
             pnorm(tau-epsilon, mean_posterior, sqrt(1/(n_prior*tau_posterior)+1/(n*tau_posterior))),
             1-pbeta(tau+epsilon, mean_posterior, sqrt(1/(n_prior*tau_posterior)+1/(n*tau_posterior))))
    }
  }
  lapply(armls, get_metric)
}

draw_arm_tts <- function(x) {
  x <- x/sum(x) # make error probabilities proportional
  sample(1:length(x), 1, prob = x)
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

