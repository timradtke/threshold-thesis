# KL version of the APT

#means <- c(0.55, 0.38, 0.95)
#system.time(res_kl <- KL_bandit(means = means, K = 3, rounds = 5000,
#                                tau = 0.5, epsilon = 0.1,
#                                variances = NA, seed = 56,
#                                verbose = FALSE))
#lapply(res_kl$arm_list, length)
#res_kl$means
#means
#plot(x = c(0,5000), y = c(0,1), type = "n")
#lines(res_kl$mean_storage[,1])
#lines(res_kl$mean_storage[,2], col = "blue")
#lines(res_kl$mean_storage[,3], col = "red")

#mean_twenty <- c(0.0005, 0.0005, 0.001, 0.001, 0.005,
#                 0.3, 0.4, 0.5, 0.6, 0.7,
#                 0.01, 0.02, 0.02,
#                 0.06, 0.06, 0.07,
#                 0.035, 0.035, 0.045, 0.045)
#system.time(res_kl <- KL_bandit(means = mean_twenty, K = 20, rounds = 1000,
#                                tau = 0.04, epsilon = 0.01,
#                                variances = NA, seed = 56,
#                                verbose = FALSE))
#lapply(res_kl$arm_list, length)
#res_kl$means
#means

KL_bandit <- function(means, K = 4, rounds = 5000, tau, epsilon,
                      variances = NA, verbose = FALSE, seed = NA,
                      at_tau = FALSE
) {
  
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
    if(at_tau) {
      next_arm <- get_min(unlist(get_next_arm_kl_at_tau(arm_list, tau = tau,
                                                 epsilon = epsilon,
                                                 variances = variances,
                                                 current_round = i)))
    } else {
      next_arm <- get_min(unlist(get_next_arm_kl(arm_list, tau = tau,
                                                 epsilon = epsilon,
                                                 variances = variances)))
    }
    
    
    if(verbose) message("arm selecting done")
    if(verbose) message(next_arm)
    arm_list[[next_arm]] <- c(arm_list[[next_arm]], pull_arm(k = next_arm, 
                                                             means = means,
                                                             variances = variances))
    mean_storage <- rbind(mean_storage, unlist(lapply(arm_list, mean)))
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

kl_ber <- function(x,y) {
  x*log(x/y) + (1-x)*log((1-x)/(1-y))
}

get_next_arm_kl <- function(armls, tau, epsilon, variances) {
  if(is.na(variances)) {
    # do bernoulli version
    get_metric <- function(x) {
      xhat <- mean(x)
      ni <- length(x)
      KLhat <- ifelse(xhat >= tau,
                      kl_ber(xhat, tau-epsilon),
                      kl_ber(xhat, tau+epsilon))
      return(sqrt(ni) * ifelse(is.na(KLhat),0,KLhat))
    }
  } else {
    # do gaussian version
    get_metric <- function(x) {
      xhat <- mean(x)
      ni <- length(x)
      return(sqrt(ni) * (abs(xhat-tau)+epsilon))
    }
  }
  lapply(armls, get_metric)
}

get_next_arm_kl_at_tau <- function(armls, tau, epsilon, variances, current_round) {
  if(is.na(variances)) {
    # do bernoulli version
    get_metric <- function(x) {
      xhat <- mean(x)
      ni <- length(x)
      KLhat <- kl_ber(xhat, tau) - 1/current_round
      return(ni * ifelse(is.na(KLhat),0,KLhat))
    }
  } else {
    # do gaussian version
    get_metric <- function(x) {
      xhat <- mean(x)
      ni <- length(x)
      return(sqrt(ni) * (abs(xhat-tau)+epsilon))
    }
  }
  lapply(armls, get_metric)
}

get_min <- function(x) {
  mini <- min(x)
  ifelse(sum(x==mini) == 1, which.min(x),
         sample(which(x==mini)))
}
