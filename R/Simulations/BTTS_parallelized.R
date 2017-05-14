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

interval_prob <- function(tau, epsilon, a, b) {
  pbeta(tau+epsilon, a, b) - pbeta(tau-epsilon, a, b)
}

run_BTTS_bandit <- function(data, tau, eps, rounds = 500, seed = NA,
                            a = 1, b = 1) {
  
  variants <- sort(unique(data[,"variant"]))
  k <- length(variants)
  storage <- data.frame(variant = NA, round = NA,
                        obs_index = NA, obs = NA)
  current_row <- 0
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
  
  original_data <- data
  
  if(!is.null(seed)) set.seed(seed)
  
  for (i in 1:rounds) {
    next_variant <- sample(variants, 1, prob = probs_prop)
    
    next_data <- get_next_row(data, next_variant, current_row)
    storage[i,1] <- next_data$variant
    storage[i,2] <- i
    storage[i,3] <- next_data$next_row
    storage[i,4] <- next_data$relevant_data$success[1]
    data <- next_data$relevant_data
    
    # recompute probabilities
    successes[next_data$variant] <- successes[next_data$variant] +
      next_data$relevant_data$success[1]
    trials[next_data$variant] <- trials[next_data$variant] + 1
    probs[next_data$variant] <- interval_prob(tau = tau, epsilon = eps,
                                              a = a + successes[next_data$variant],
                                              b = b + trials[next_data$variant] - successes[next_data$variant])
    
    probs_prop <- probs/sum(probs)
    
    successes_storage[i,] <- successes
    trials_storage[i,] <- trials
  }
  
  return(list(successes = successes,
              trials = trials,
              probs = probs,
              mean = successes/trials,
              set = successes/trials > tau,
              storage = storage,
              successes_storage = successes_storage,
              trials_storage = trials_storage))
}

run_BTTS_simulation <- function(reps, rounds, threshold) {
  
  # prepare for parallel processing of 'foreach' loop
  cl <- makeCluster( max(1,detectCores()-1)) # for example: if 4 cores, use 3
  registerDoParallel(cl)
  required_packages <- c("plyr", "dplyr", "purrr")
  
  ff <- foreach(i = 1:reps, 
            .packages = required_packages,
            .export = c("simulate_data_bernoulli", "get_next_row",
                      "interval_prob", "run_BTTS_bandit", "fake_function")) %dopar% {

    #set.seed(511+i)
    experiment <- simulate_data_bernoulli(variants = 10, 
                                          means = c(0.4-(0.2^(1:4)), 0.45, 0.55, 0.6 + 0.2^(4:1)),
                                          rounds = 200000)$experiment
    
    btts <- run_BTTS_bandit(experiment, threshold, 0.1, rounds, seed = NULL,
                            a = 1, b = 1)
    
    successes_storage <- btts$successes_storage
    trials_storage <- btts$trials_storage
    mean_storage <- btts$successes_storage / btts$trials_storage
    set_storage <- (btts$successes_storage/btts$trials_storage) > threshold
    
    loss <- function(x) {return(1-(sum(x[1:4]) == 0 & sum(x[7:10]) == 4))}
    loss_storage <- as.vector(aaply(set_storage, 1, loss))
    
    result <- list(successes_storage = successes_storage,
                   trials_storage = trials_storage,
                   mean_storage = mean_storage,
                   set_storage = set_storage,
                   loss_storage = loss_storage)
  }
  
  return(ff)
}


gc()
system.time(btts_sim <- run_BTTS_simulation(10, 500, 0.5))
tail(round(btts_sim[[1]]$mean_storage,3))
tail(round(btts_sim[[2]]$mean_storage,3))
btts_sim[[1]]$loss_storage
btts_sim[[2]]$loss_storage

library(beepr)
gc()
system.time(btts_sim1 <- run_BTTS_simulation(1000, 500, 0.5))
#     user    system   elapsed 
#   21.316    29.042 13465.151
save(btts_sim1, file = "/Users/timradtke/Dropbox/1Master/Master Thesis/1Code/Simulations/results/btts_sim1.Rda")
beep(2)

gc()
system.time(btts_sim2 <- run_BTTS_simulation(500, 500, 0.5))
#    user   system  elapsed 
#  12.489   16.796 6702.540 
save(btts_sim2, file = "/Users/timradtke/Dropbox/1Master/Master Thesis/1Code/Simulations/results/btts_sim2.Rda")
beep(2)

gc()
system.time(btts_sim3 <- run_BTTS_simulation(1500, 500, 0.5))
#  user    system   elapsed 
#31.427    46.477 19744.650 
save(btts_sim3, file = "/Users/timradtke/Dropbox/1Master/Master Thesis/1Code/Simulations/results/btts_sim3.Rda")
beep(2)

gc()
system.time(btts_sim_geom <- run_BTTS_simulation(2000, 500, 0.5))
#  user    system   elapsed 
#36.235    53.529 25366.544 
save(btts_sim_geom, file = "/Users/timradtke/Dropbox/1Master/Master Thesis/1Code/Simulations/results/btts_sim_geom.Rda")
beep(2)


#experiment <- simulate_data_bernoulli(variants = 10, 
#means = c(0.1,0.1,0.1,
#          0.35,0.45,0.55,0.65,
#          0.9,0.9,0.9),
#rounds = 200000)$experiment

btts_sim1[[500]]$loss_storage
sim1_loss_storage <- matrix(ncol = 1000, nrow = 500)
for (i in 1:1000) {
  sim1_loss_storage[,i] <- btts_sim1[[i]]$loss_storage
}

sim2_loss_storage <- matrix(ncol = 500, nrow = 500)
for (i in 1:500) {
  sim2_loss_storage[,i] <- btts_sim2[[i]]$loss_storage
}

sim3_loss_storage <- matrix(ncol = 1500, nrow = 500)
for (i in 1:1500) {
  sim3_loss_storage[,i] <- btts_sim3[[i]]$loss_storage
}

plot(1:500, (rowMeans(cbind(sim1_loss_storage, sim2_loss_storage, 
                               sim3_loss_storage))), type = "l")
abline(h=0.05, lty = 2)
plot(1:500, log(rowMeans(cbind(sim1_loss_storage, sim2_loss_storage, 
                               sim3_loss_storage))), type = "l", ylim = c(-8,0),
     ylab = "log(1-est. probability of success)", lwd = 1.5)
abline(h=log(0.05), lty = 2)


sim_geom_loss_storage <- matrix(ncol = 2000, nrow = 500)
for (i in 1:2000) {
  sim_geom_loss_storage[,i] <- btts_sim_geom[[i]]$loss_storage
}

plot(1:500, log(rowMeans(sim_geom_loss_storage)), type = "l", ylim = c(-3,0),
     ylab = "log(1-est. probability of success)", xlab = "horizon", lwd = 1.5)
abline(h=log(0.05), lty = 2)
plot(1:500, (rowMeans(sim_geom_loss_storage)), type = "l")
abline(h=(0.05), lty = 2)
