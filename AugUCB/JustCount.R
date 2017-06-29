# Just count

#mean_twenty <- c(0.0005, 0.0005, 0.001, 0.001, 0.005,
#                 0.3, 0.4, 0.5, 0.6, 0.7,
#                 0.01, 0.02, 0.02,
#                 0.06, 0.06, 0.07,
#                 0.035, 0.035, 0.045, 0.045)
#system.time(res_jc <- JC_bandit(means = mean_twenty,
#                                tau = 0.04, epsilon = 0.01,
#                                alpha = 0.04, beta = 0.96,
#                                K = 20, rounds = 5000, seed = 58))
#lapply(res_jc$arm_list, length)
#res_jc$means
#jc_means <- as.data.frame(res_jc$mean_storage)
#jc_means$round <- 1:dim(jc_means)[1]
#jc_means_long <- gather(jc_means, key = "Arm", value = "Mean", -round)
#ggplot(jc_means_long, aes(x = round, y = Mean, color = Arm)) +
#  geom_line()



JC_bandit <- function(means, K, rounds, tau, epsilon,
                      alpha, beta, verbose = FALSE, seed = NA) {
  
  if(!is.na(seed)) set.seed(seed)
  
  # initialize by pulling each arm once
  arm_list <- NA
  for(i in 1:K) {
    # Ideally it should also take an argument specifying how
    # the sample is generated (e.g., from existing data, bernoulli or other distr.)
    arm_list[[i]] <- pull_arm_JC(k = i, means = means)
  }
  arm_list <- as.list(arm_list)
  # initialize mean storage and the counter
  mean_storage <- matrix(unlist(lapply(arm_list, mean)), nrow = 1)
  counter <- K
  
  for(i in (K+1):rounds) {
    if(verbose) message(paste("this is round", i))
    next_arm <- get_min(unlist(get_next_arm_JC(arm_list, tau = tau,
                                               epsilon = epsilon,
                                               alpha = alpha,
                                               beta = beta)))
    
    if(verbose) message("arm selecting done")
    arm_list[[next_arm]] <- c(arm_list[[next_arm]], pull_arm_JC(k = next_arm, 
                                                                means = means))
    mean_storage <- rbind(mean_storage, unlist(lapply(arm_list,
                                                      get_posterior_mean_JC,
                                                      alpha = alpha, beta = beta)))
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
pull_arm_JC <- function(k, means) {
  rbinom(1,1,means[k])
}

draw_arm_tts <- function(x) {
  x <- x/sum(x) # make error probabilities proportional
  sample(1:length(x), 1, prob = x)
}

get_next_arm_JC <- function(armls, tau, epsilon,
                            alpha, beta) {
  # depending on whether variances are available, do method for bernoulli or gaussian
  # only Bernoulli possible for now
  get_metric <- function(x) {
    alpha_prime <- sum(x)+alpha
    beta_prime <- length(x)+beta-sum(x)
    posterior_mean <- alpha_prime/(alpha_prime+beta_prime)
    # depending on current mean estimate, give probability of error
    ifelse(posterior_mean >= tau,
           ceiling(alpha_prime/(tau-epsilon)-beta_prime),
           ceiling(beta_prime*(tau+epsilon)-alpha_prime))
  }
  lapply(armls, get_metric)
}

get_min <- function(x) {
  mini <- min(x)
  ifelse(sum(x==mini) == 1, which.min(x),
         sample(which(x==mini)))
}

get_posterior_mean_JC <- function(x, alpha, beta) {
  # depending on whether variances are available,
  # do method for bernoulli or gaussian
  alpha_prime <- sum(x)+alpha
  beta_prime <- length(x)+beta-sum(x)
  alpha_prime/(alpha_prime+beta_prime)
}

