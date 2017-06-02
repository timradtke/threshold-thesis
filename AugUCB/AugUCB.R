# AugUCB

# arguments:
# budget rounds
# parameter rho
# threshold tau
# number of arms K
# sim_fun (takes as input an arm index and the arguments passed to it ...)

AugUCB <- function(means, K = 4, rounds = 5000, rho = 1/3, tau = 0.5#, sim_fun, ...
                   ) {
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
    arm_list[[i]] <- pull_arm(k = i, means = means)
  }
  arm_list <- as.list(arm_list)
  
  for(i in (K+1):rounds) {
    message(paste("this is round", i))
    next_arm <- get_min(unlist(get_next_arm_augucb(arm_list, tau = tau, rho = rho,
                                              psi = psi, rounds = rounds, 
                                              epsilon = epsilon,
                                              active_set = B)))
    
    message("arm selecting done")
    arm_list[[next_arm]] <- c(arm_list[[next_arm]], pull_arm(k = next_arm, 
                                                             means = means))
    
    message("arm pulling done")
    # delete arms should return a logical index for each arm in active set
    B <- B[!unlist(delete_arms(arm_list, tau = tau, rho = rho,
                               psi = psi, rounds = rounds, epsilon = epsilon,
                               active_set = B))]
    
    message("B done")
    message(B)
    
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
              true_means = means))
}

res <- AugUCB(means = c(0.1,0.2,0.45,0.7,0.9,0.55,0.38,0.95), K = 8, rounds = 10000)
lapply(res$arm_list, length)
res$active_set
res$means
means

psi

# Ideally it should also take an argument specifying how
# the sample is generated (e.g., from existing data, bernoulli or other distr.)
pull_arm <- function(k, means) {
  rbinom(1,1,means[k])
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
  return(lapply(armls[B], eliminate))
}

# randomize arm if several arms have lowest value
get_min <- function(x) {
  mini <- min(x)
  ifelse(sum(x==mini) == 1, which.min(x),
         sample(which(x==mini)))
}


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

K <- 2^(1:8)
rounds <- 2^(8:16)

paramss <- expand.grid(K, rounds)
paramss %>% rename(K = Var1, rounds = Var2) %>%
  mutate(rounds_per_K = rounds/K,
         psi0 = get_psi0(rounds, K),
         l0 = get_l0(rounds, psi0),
         N0 = get_N0(K, l0),
         si = get_si_factor(psi0, rounds)) %>%
  mutate(psi0 = round(psi0,3)) %>%
  arrange(rounds_per_K)

log(3/16*2*log(2))
log(3/16*4*log(4))
log(3/16*8*log(8))
