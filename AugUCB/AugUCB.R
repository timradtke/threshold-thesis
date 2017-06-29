# AugUCB
# Thresholding Bandits with Augmented UCB
# Mukherjee et al. (2017)
# arXiv:1704.02281v1

# arguments:
# budget rounds
# parameter rho
# threshold tau
# number of arms K
# sim_fun (takes as input an arm index and the arguments passed to it ...)

AugUCB <- function(means, K = 4, rounds = 5000, rho = 1/3, tau = 0.5,
                   variances = NA, verbose = FALSE, seed = NA
                   ) {
  B <- 1:K # the active set of arms
  m <- 0
  epsilon <- 1
  e <- exp(1) # ????????
  M <- floor(0.5*log2(rounds/e))
  psi <- rounds*epsilon/(128 * (log(3/16 * K * log(K)))^2)
  l <- ceiling(2*psi*log(rounds * epsilon)/epsilon)
  N <- K * l
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
    next_arm <- get_min(unlist(get_next_arm_augucb(arm_list, tau = tau, rho = rho,
                                              psi = psi, rounds = rounds, 
                                              epsilon = epsilon,
                                              active_set = B)))
    
    if(verbose) message("arm selecting done")
    arm_list[[next_arm]] <- c(arm_list[[next_arm]], pull_arm(k = next_arm, 
                                                             means = means,
                                                             variances = variances))
    mean_storage <- rbind(mean_storage, unlist(lapply(arm_list, mean)))
    counter <- counter + 1
    
    if(verbose) message("arm pulling done")
    # delete arms should return a logical index for each arm in active set
    B <- B[!unlist(delete_arms(arm_list, tau = tau, rho = rho,
                               psi = psi, rounds = rounds, epsilon = epsilon,
                               active_set = B))]
    
    if(verbose) message("B done")
    if(verbose) message(B)
    
    if((i >= N) && (m <= M)) {
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
              active_set = B,
              true_means = means,
              counter = counter,
              mean_storage = mean_storage))
}

#res <- AugUCB(means = c(0.1,0.2,0.45,0.7,0.9,0.55,0.38,0.95), K = 8, rounds = 10000)
#means <- c(0.55, 0.38, 0.95)
#variances <- c(0.5, 0.9, 1.7)
#system.time(res2 <- AugUCB(means = means, K = 3, rounds = 5000, 
#                          variances = variances, seed = 54))
#   user  system elapsed 
#  2.025   0.040   2.076 
#lapply(res2$arm_list, length)
#res$active_set
#res$means
#means
#plot(res2$mean_storage[,1], type = "l")
#plot(res2$mean_storage[,2], type = "l")
#plot(res2$mean_storage[,3], type = "l")


# Ideally it should also take an argument specifying how
# the sample is generated (e.g., from existing data, bernoulli or other distr.)
pull_arm <- function(k, means, variances = NA) {
  if(is.na(variances)) {
    rbinom(1,1,means[k])
  } else {
    rnorm(1,means[k], variances[k])
  }
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

# randomize arm if several arms have lowest value
get_min <- function(x) {
  mini <- min(x)
  ifelse(sum(x==mini) == 1, which.min(x),
         sample(which(x==mini)))
}


simulate_bandit_data_ber <- function(K, means, N = 50000, seed = NA) {
  if(!is.na(seed)) set.seed(seed)
  data <- data.frame(arms = sample(1:K, N, replace = TRUE))
  data$obs <- apply(data, 1, function(x) rbinom(1,1,means[x]))
  return(data)
}
#simdata <- simulate_bandit_data_ber(3, means)
#simdata %>% group_by(arms) %>% summarize(mean = mean(obs))

simulate_bandit_data_norm <- function(K, means, variances, N = 50000, seed = NA) {
  if(!is.na(seed)) set.seed(seed)
  data <- data.frame(arms = sample(1:K, N, replace = TRUE))
  data$obs <- apply(data, 1, function(x) rnorm(1,means[x],sqrt(variances[x])))
  return(data)
}
#variances <- c(0.5, 0.9, 1.7)
#simdata <- simulate_bandit_data_norm(3, means, variances)
#simdata %>% group_by(arms) %>% summarize(mean = mean(obs))
#simdata %>% group_by(arms) %>% summarize(var = var(obs))

# take a data set, and go chronologically through it until the next
# observation for arm k; take that observation and discard remaining data
pull_arm_from_data <- function(k, data) {
  # test whether some data is remaining
  if(is.na(data) || dim(data)[1]==0) stop("No data remaining.")
  
  #data <- data.frame("arm" = c(1,3,3,2,1,1,3), "obs" = c(1,0,0,0,1,1,1))
  # expect that data is a data frame where the first column has the arm names
  # and second column has the observation (reward)
  next_row <- which(data[1] == k)[1]
  if(!is.na(next_row)) {
    next_obs <- data[next_row, 2]
    next_data <- data[-c(1:next_row),]
  } else { # when there was no observation left for arm k
    next_obs <- NA
    next_data <- NA
  }
  
  return(list(obs = next_obs,
              data = next_data))
}

# data is a data frame with first column giving the arms, second column the observations
# both columns need to be numeric
AugUCB_from_data <- function(data, K = 4, rounds = 5000, rho = 1/3, tau = 0.5,
                             seed = NA, verbose = FALSE) {
  if(length(unique(data[[1]])) != K) stop("Data does not have K arms.")
  
  B <- 1:K # the active set of arms
  m <- 0
  epsilon <- 1
  e <- exp(1) # ????????
  M <- floor(0.5*log2(rounds/e))
  psi <- rounds*epsilon/(128 * (log(3/16 * K * log(K)))^2)
  l <- ceiling(2*psi*log(rounds * epsilon)/epsilon)
  N <- K * l
  
  # initialize by pulling each arm once
  arm_list <- NA
  for(i in 1:K) {
    # Ideally it should also take an argument specifying how
    # the sample is generated (e.g., from existing data, bernoulli or other distr.)
    new_data <- pull_arm_from_data(k = i, data = data)
    arm_list[[i]] <- new_data$obs
    data <- new_data$data
  }
  arm_list <- as.list(arm_list)
  # initialize mean storage and the counter
  mean_storage <- matrix(unlist(lapply(arm_list, mean)), nrow = 1)
  counter <- K
  
  for(i in (K+1):rounds) {
    if(verbose) message(paste("this is round", i))
    if(!is.na(seed)) set.seed(seed+i)
    next_arm <- get_min(unlist(get_next_arm_augucb(arm_list, tau = tau, rho = rho,
                                              psi = psi, rounds = rounds, 
                                              epsilon = epsilon,
                                              active_set = B)))
    
    if(verbose) message("arm selecting done")
    new_data <- pull_arm_from_data(k = next_arm, data = data)
    arm_list[[next_arm]] <- c(arm_list[[next_arm]], new_data$obs)
    data <- new_data$data
    counter <- counter + 1
    
    # keep track of the estimates at each round
    mean_storage <- rbind(mean_storage, unlist(lapply(arm_list, mean)))
    if(is.na(data) || dim(data)[1]==0) break
    
    if(verbose) message("arm pulling done")
    # delete arms should return a logical index for each arm in active set
    B <- B[!unlist(delete_arms(arm_list, tau = tau, rho = rho,
                               psi = psi, rounds = rounds, epsilon = epsilon,
                               active_set = B))]
    
    if(verbose) message("B done")
    if(verbose) message(B)
    
    if((i >= N) && (m <= M)) {
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
              active_set = B,
              counter = counter,
              mean_storage = mean_storage))
}

#means <- c(0.55, 0.38, 0.95)
#variances <- c(0.5, 0.9, 1.7)
#simdata <- simulate_bandit_data_norm(3, means, variances)
#system.time(res <- AugUCB_from_data(simdata, K = 3, rounds = 5000, 
#                                    rho = 1/3, tau = 0.5, seed = 54))
#   user  system elapsed 
#132.082  13.134 146.622 
#res$mean_storage
#res$means
#lapply(res$arm_list, length)
##############################################################
# analyze behavior of parameters

get_psi0 <- function(rounds, K) {
  rounds / (128*(log(3/16*K*log(K)))^2)
}

get_l0 <- function(rounds, psi) {
  ceiling(2*psi*log(rounds))
}

get_N0 <- function(K, l) K*l

get_si_factor <- function(psi, rounds) {
  sqrt(1/3*psi*log(rounds)/4/1)
}

#K <- 2^(1:8)
#rounds <- 2^(8:16)

#paramss <- expand.grid(K, rounds)
#paramss %>% rename(K = Var1, rounds = Var2) %>%
#  mutate(rounds_per_K = rounds/K,
#         psi0 = get_psi0(rounds, K),
#         l0 = get_l0(rounds, psi0),
#         N0 = get_N0(K, l0),
#         si = get_si_factor(psi0, rounds)) %>%
#  mutate(psi0 = round(psi0,3)) %>%
#  arrange(rounds_per_K)

#log(3/16*2*log(2))
#log(3/16*4*log(4))
#log(3/16*8*log(8))
