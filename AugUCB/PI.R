# Probability of Improvement-like (actually, probability of error)

#means <- c(0.55, 0.38, 0.95)
#variances <- c(0.5, 0.9, 1.7)
#system.time(res_pi <- PI_bandit(means = means, variances = variances,
#                                alpha = 10, beta = 5000,
#                                mean_prior = 0.5, n_prior = 10,
#                                K = 3, rounds = 5000,
#                                seed = 54))
# Bernoulli version
#means <- c(0.02, 0.06, 0.065, 0.07, 0.01, 0.03)
#system.time(res_pi <- PI_bandit(means = means, alpha = 0.04*28, beta = 0.96*28,
#                                K = 6, rounds = 5000,
#                                tau = 0.04, epsilon = 0.011,
#                                variances = NA, seed = 54))

#means <- c(0.55, 0.38, 0.95)
#system.time(res_pi <- PI_bandit(means = means, alpha = 1, beta = 1,
#                                K = 3, rounds = 5000,
#                                variances = NA, seed = 54))
# user  system elapsed 
#0.401   0.079   0.481
#lapply(res_pi$arm_list, length)
#res_pi$means
#means
#plot(x = c(0,5000), y = c(0,0.1), type = "n")
#lines(res_pi$mean_storage[,1])
#lines(res_pi$mean_storage[,2], col = "blue")
#lines(res_pi$mean_storage[,3], col = "red")
#lines(res_pi$mean_storage[,4], col = "red")
#lines(res_pi$mean_storage[,5], col = "red")
#lines(res_pi$mean_storage[,6], col = "red")

#mean_twenty <- c(0.0005, 0.0005, 0.001, 0.001, 0.005,
#                 0.3, 0.4, 0.5, 0.6, 0.7,
#                 0.01, 0.02, 0.02,
#                 0.06, 0.06, 0.07,
#                 0.035, 0.035, 0.045, 0.045)
#system.time(res_pi <- PI_bandit(means = mean_twenty,
#                                tau = 0.04, epsilon = 0.01,
#                                alpha = 0.04, beta = 0.96,
#                                K = 20, rounds = 2000,
#                                variances = NA, seed = 58))
#lapply(res_pi$arm_list, length)
#res_pi$means
#pi_means <- as.data.frame(res_pi$mean_storage)
#pi_means$round <- 1:dim(pi_means)[1]
#pi_means_long <- gather(pi_means, key = "Arm", value = "Mean", -round)
#ggplot(pi_means_long, aes(x = round, y = Mean, color = Arm)) +
#  geom_line()



PI_bandit <- function(means, K, rounds, tau, epsilon,
                      alpha, beta, mean_prior = NA, n_prior = NA,
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
    # NOTE THE MINUS SIGN; find the arm that MAXIMIZES the error probability
    next_arm <- get_min(-unlist(get_next_arm_PI(arm_list, tau = tau,
                                               epsilon = epsilon,
                                               alpha = alpha,
                                               beta = beta,
                                               mean_prior = mean_prior,
                                               n_prior = n_prior)))
    
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

# for information on the conjugate prior for the normal distribution
# https://people.eecs.berkeley.edu/~jordan/courses/260-spring10/lectures/lecture5.pdf
# https://en.wikipedia.org/wiki/Gamma_distribution

# we simply use the mean of the gamma distribution as the posterior estimate for
# the variance
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
                            alpha, beta, mean_prior, n_prior) {
  # depending on whether variances are available, do method for bernoulli or gaussian
  if(is.na(mean_prior)) {
    # Bernoulli
    get_metric <- function(x) {
      alpha_prime <- sum(x)+alpha
      beta_prime <- length(x)+beta-sum(x)
      posterior_mean <- alpha_prime/(alpha_prime+beta_prime)
      # depending on current mean estimate, give probability of error
      ifelse(posterior_mean >= tau,
             pbeta(tau-epsilon, alpha_prime, beta_prime),
             1-pbeta(tau+epsilon, alpha_prime, beta_prime))
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

