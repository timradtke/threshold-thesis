# simulation_functions

library(plyr); library(dplyr); library(doParallel); library(parallel)

simulate_data_bernoulli <- function(variants = 3,
                                    means = c(0.5, 0.55, 0.6), 
                                    rounds = 10000,
                                    equal_splits = TRUE,
                                    splits = NULL,
                                    seed = NULL) {
  if(variants != length(means)) {
    stop("Number of variants does not equal number of means.")
  }
  
  if(!is.null(splits)) {
    if(sum(splits != 1)) {
      stop("Splits don't sum up to 1.")
    }
  }
  
  if(is.null(splits)) {
    splits <- rep(1/variants, times = variants)
  }
  
  if(variants != length(splits)) {
    stop("Number of variants does not equal number of splits.")
  }
  
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  variant_names <- toupper(letters[1:variants])
  variant_seq <- sample(variant_names, size = rounds, 
                        replace = TRUE, prob = splits)
  
  params <- data.frame(variant = variant_names, mean = means,
                       stringsAsFactors = FALSE)
  experiment <- data.frame(variant = variant_seq,
                           stringsAsFactors = FALSE)
  experiment <- dplyr::left_join(experiment, params)
  
  experiment <- experiment %>% group_by(variant) %>%
    dplyr::mutate(success = rbinom(n(), 1, mean),
                  trials = n(),
                  successes = sum(success),
                  success_rate = successes/trials) %>%
    ungroup() %>% as.data.frame()
  
  experiment_summary <- experiment %>%
    select(variant, mean, trials, successes, success_rate) %>%
    distinct() %>%
    dplyr::arrange(variant)
  
  return(list(
    experiment = experiment,
    experiment_summary = experiment_summary
  ))
}


# assumes that the argument data is a data frame with first column containing
# the variants that are checked for;
get_next_row <- function(data, variant, current_row) {
  n <- dim(data)[1]
  # if there is an observation remaining
  if(sum(data[(current_row+1):n,1] == variant) > 0) {
    next_row <- current_row + min((1:n)[data[(current_row+1):n,1] == variant])
    relevant_data <- data[next_row:n,]
    end <- FALSE
  } else { # else (when there is no observation remaining)
    next_row <- NA
    relevant_data <- NA
    end <- TRUE
  }
  
  return(list(
    data = data,
    variant = variant,
    next_row = next_row,
    relevant_data = relevant_data,
    end = end
  ))
}

interval_prob2 <- function(tau, epsilon, a, b) {
  ifelse(a/(a+b) >= tau,
         pbeta(tau - epsilon, a, b),
         1 - pbeta(tau + epsilon, a, b))
}

fake_function <- function(x) return(x)

run_BTTS2_bandit <- function(means, tau, eps, rounds = 500, seed = NA,
                             a = 1, b = 1) {
  
  variants <- 1:length(means)
  k <- length(variants)
  storage <- data.frame(round = 1:rounds, variant = NA, obs = NA)
  
  #current_row <- 0
  successes <- rep(0, times = k)
  trials <- rep(0, times = k)
  probs <- rep(0.5, times = k)
  probs_prop <- probs/sum(probs)
  successes_storage <- matrix(ncol = k, nrow = 500)
  trials_storage <- matrix(ncol = k, nrow = 500)
  names(successes) <- variants
  names(trials) <- variants
  names(probs) <- variants
  names(probs_prop) <- variants
  
  #original_data <- data
  
  if(!is.null(seed)) set.seed(seed)
  
  for (i in 1:rounds) {
    next_variant <- sample(variants, 1, prob = probs_prop)
    next_obs <- purrr::rbernoulli(1, p = means[next_variant])
    
    storage[i,2] <- next_variant
    storage[i,3] <- next_obs
    
    # recompute probabilities
    successes[next_variant] <- successes[next_variant] + next_obs
    trials[next_variant] <- trials[next_variant] + 1
    probs[next_variant] <- interval_prob2(tau = tau, epsilon = eps,
                                          a = a + successes[next_variant],
                                          b = b + trials[next_variant] - 
                                            successes[next_variant])
    probs_prop <- probs/sum(probs)
    
    successes_storage[i,] <- successes
    trials_storage[i,] <- trials
  }
  
  return(list(successes = successes,
              trials = trials,
              probs = probs,
              means = means,
              means_est = (a+successes)/(a+b+trials),
              set = (a+successes)/(a+b+trials) >= tau,
              storage = storage,
              successes_storage = successes_storage,
              trials_storage = trials_storage,
              alpha_storage = successes_storage + a,
              beta_storage = trials_storage + b - successes_storage))
}

run_BTTS2_simulation <- function(reps, rounds, threshold, epsilon, 
                                 means, a, b, loss_function) {
  
  # prepare for parallel processing of 'foreach' loop
  cl <- makeCluster(max(1, detectCores() - 1)) # for example: if 4 cores, use 3
  registerDoParallel(cl)
  required_packages <- c("plyr", "dplyr", "purrr")
  
  ff <- foreach(i = 1:reps, 
                .packages = required_packages,
                .export = c("interval_prob2", "run_BTTS2_bandit", 
                            "fake_function")) %dopar% {
    
    btts <- run_BTTS2_bandit(means = means, tau = threshold, 
                            eps = epsilon, rounds = rounds, seed = NULL,
                            a = a, b = b)
    
    successes_storage <- btts$successes_storage
    trials_storage <- btts$trials_storage
    alpha_storage <- btts$alpha_storage
    beta_storage <- btts$beta_storage
    mean_storage <- btts$alpha_storage / (btts$alpha_storage+btts$beta_storage)
    set_storage <- (btts$alpha_storage / 
                      (btts$alpha_storage+btts$beta_storage)) >= threshold
    
    
    loss_storage <- as.vector(aaply(set_storage, 1, loss_function))
    
    result <- list(successes_storage = successes_storage,
                   trials_storage = trials_storage,
                   alpha_storage = alpha_storage,
                   beta_storage = beta_storage,
                   mean_storage = mean_storage,
                   set_storage = set_storage,
                   loss_storage = loss_storage)
  }
  
  return(ff)
}


############# EXPERIMENT 1 #####################################################

means <- c(0.1,0.1,0.1, 0.35,0.45, 0.55,0.65, 0.9,0.9,0.9)
loss <- function(x) {return(1-(sum(x[1:4]) == 0 & sum(x[7:10]) == 4))}

# test first
gc()
system.time(btts_sim <- run_BTTS2_simulation(reps = 100, rounds = 500, 
                                            threshold = 0.5, epsilon = 0.1,
                                            means = means, a = 1, b = 1,
                                            loss_function = loss))
tail(round(btts_sim[[1]]$mean_storage,3))
tail(round(btts_sim[[2]]$mean_storage,3))
btts_sim[[1]]$loss_storage
btts_sim[[2]]$loss_storage

# and now for real
library(beepr)
gc()
system.time(btts2_sim1 <- run_BTTS2_simulation(reps = 5000, rounds = 500, 
                                             threshold = 0.5, epsilon = 0.1,
                                             means = means, a = 1, b = 1,
                                             loss_function = loss))
#user  system elapsed 
#5.290   2.689 214.312
save(btts2_sim1, file = "/Users/timradtke/Dropbox/1Master/Master Thesis/1Code/Simulations/results/btts2_sim1.Rda")
beep(2)

gc()
means_arith <- c(0.2+(0:3)*0.05, 0.45, 0.55, 0.6 + 0.05*(0:3))
system.time(btts2_sim_arith <- run_BTTS2_simulation(reps = 5000, rounds = 500, 
                                                   threshold = 0.5, epsilon = 0.1,
                                                   means = means_arith, 
                                                   a = 1, b = 1,
                                                   loss_function = loss))
# user  system elapsed 
# 5.422   2.952 224.228 
save(btts2_sim_arith, file = "/Users/timradtke/Dropbox/1Master/Master Thesis/1Code/Simulations/results/btts2_sim_arith.Rda")
beep(2)

gc()
means_geom <- c(0.4-(0.2^(1:4)), 0.45, 0.55, 0.6 + 0.2^(4:1))
system.time(btts2_sim_geom <- run_BTTS2_simulation(reps = 5000, rounds = 500, 
                                                  threshold = 0.5, epsilon = 0.1,
                                                  means = means_geom, 
                                                  a = 1, b = 1,
                                                  loss_function = loss))
# user  system elapsed 
# 5.422   2.952 224.228 
save(btts2_sim_geom, file = "/Users/timradtke/Dropbox/1Master/Master Thesis/1Code/Simulations/results/btts2_sim_geom.Rda")
beep(2)


### New, own test

gc()
means_low <- c(0.01, 0.01, 0.03, 0.04, 0.05, 0.06)
loss_low <- function(x) {return(1-(sum(x[1:2]) == 0 & sum(x[4:6]) == 3))}

system.time(btts2_sim_low <- run_BTTS2_simulation(reps = 5000, rounds = 500, 
                                                   threshold = 0.03, epsilon = 0.005,
                                                   means = means_low, 
                                                   a = 1, b = 1,
                                                   loss_function = loss_low))
# user  system elapsed 
# 5.422   2.952 224.228 
save(btts2_sim_low, file = "/Users/timradtke/Dropbox/1Master/Master Thesis/1Code/Simulations/results/btts2_sim_low.Rda")
beep(2)



########## Plot the results ############################################

load(file = "/Users/timradtke/Dropbox/1Master/Master Thesis/1Code/Simulations/results/btts2_sim1.Rda")
load(file = "/Users/timradtke/Dropbox/1Master/Master Thesis/1Code/Simulations/results/btts2_sim_arith.Rda")
load(file = "/Users/timradtke/Dropbox/1Master/Master Thesis/1Code/Simulations/results/btts2_sim_geom.Rda")

## Three groups

sim1_loss_storage <- matrix(ncol = 5000, nrow = 500)
for (i in 1:5000) {
  sim1_loss_storage[,i] <- btts2_sim1[[i]]$loss_storage
}

success_prob <- rowMeans(sim1_loss_storage)

plot(1:500, success_prob, type = "l")
abline(h=0.05, lty = 2)
plot(1:500, log(rowMeans(sim1_loss_storage)), type = "l", ylim = c(-8,0),
     ylab = "log(1-est. probability of success)", lwd = 1.5)
abline(h=log(0.05), lty = 2)

## means_arith

sim_arith_loss_storage <- matrix(ncol = 5000, nrow = 500)
for (i in 1:5000) {
  sim_arith_loss_storage[,i] <- btts2_sim_arith[[i]]$loss_storage
}

plot(1:500, log(rowMeans(sim_arith_loss_storage)), type = "l", ylim = c(-7,0),
     ylab = "log(1-est. probability of success)", xlab = "horizon", lwd = 1.5)
abline(h=log(0.05), lty = 2)
plot(1:500, (rowMeans(sim_arith_loss_storage)), type = "l")
abline(h=(0.05), lty = 2)

## Geom. means

sim_geom_loss_storage <- matrix(ncol = 5000, nrow = 500)
for (i in 1:5000) {
  sim_geom_loss_storage[,i] <- btts2_sim_geom[[i]]$loss_storage
}

plot(1:500, log(rowMeans(sim_geom_loss_storage)), type = "l", ylim = c(-3,0),
     ylab = "log(1-est. probability of success)", xlab = "horizon", lwd = 1.5)
abline(h=log(0.05), lty = 2)
abline(h=-2.5, lty = 2)
plot(1:500, (rowMeans(sim_geom_loss_storage)), type = "l")
abline(h=(0.05), lty = 2)

## Low means

sim_low_loss_storage <- matrix(ncol = 5000, nrow = 500)
for (i in 1:5000) {
  sim_low_loss_storage[,i] <- btts2_sim_low[[i]]$loss_storage
}

plot(1:500, log(rowMeans(sim_low_loss_storage)), type = "l", ylim = c(-3,0),
     ylab = "log(1-est. probability of success)", xlab = "horizon", lwd = 1.5)
abline(h=log(0.05), lty = 2)
abline(h=-2.5, lty = 2)
plot(1:500, (rowMeans(sim_low_loss_storage)), type = "l")
abline(h=(0.05), lty = 2)
