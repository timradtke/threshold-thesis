# UNIFORM
# Common knowledge
# Our ancestors

means <- c(0.55, 0.38, 0.95)
variances <- c(0.5, 0.9, 1.7)
system.time(res_uniform <- uniform_bandit(means = means, K = 3, rounds = 5000, 
                                          variances = variances, seed = 54))
# user  system elapsed 
#0.533   0.065   0.603 
lapply(res_uniform$arm_list, length)
res_uniform$means
means
plot(x = c(0,5000), y = c(0,2), type = "n")
lines(res_uniform$mean_storage[,1])
lines(res_uniform$mean_storage[,2], col = "blue")
lines(res_uniform$mean_storage[,3], col = "red")

uniform_bandit <- function(means, K = 4, rounds = 5000,
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
    next_arm <- get_min(unlist(get_next_arm_uniform(arm_list)))
    
    if(verbose) message("arm selecting done")
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

get_next_arm_uniform <- function(armls) {
  lapply(armls, length)
}

get_min <- function(x) {
  mini <- min(x)
  ifelse(sum(x==mini) == 1, which.min(x),
         sample(which(x==mini)))
}
