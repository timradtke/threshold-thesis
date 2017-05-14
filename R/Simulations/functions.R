# simulation_functions

library(plyr); library(dplyr)

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

# Bernoulli Threshold Thompson Sampling



set.seed(512)
set.seed(1024)
set.seed(2028)
set.seed(4056)
experiment <- simulate_data_bernoulli(means = c(0.05,0.3,0.2))$experiment
btts <- run_BTTS_bandit(experiment, 0.225, 0.0125, 2000, seed = 512, a = 1, b = 1)

set.seed(512)
experiment <- simulate_data_bernoulli(variants = 10, means = c(0.1,0.1,0.1,
                                                0.35,0.45,0.55,0.65,
                                                0.9,0.9,0.9),
                                      rounds = 100000)$experiment
btts <- run_BTTS_bandit(experiment, 0.5, 0.1, 500, seed = NULL,
                        a = 1, b = 1)


btts_sim <- run_BTTS_simulation(5000, 500, 0.5)


run_BTTS_simulation <- function(reps, rounds, threshold) {
  successes_storage <- list()
  trials_storage <- list()
  mean_storage <- list()
  set_storage <- list()
  loss_storage <- matrix(ncol = reps, nrow = rounds)
  
  for(i in 1:reps) {
    set.seed(511+i)
    experiment <- simulate_data_bernoulli(variants = 10, 
                                          means = c(0.1,0.1,0.1,
                                                    0.35,0.45,0.55,0.65,
                                                    0.9,0.9,0.9),
                                          rounds = 200000)$experiment
    
    btts <- run_BTTS_bandit(experiment, threshold, 0.1, rounds, seed = NULL,
                            a = 1, b = 1)
    
    successes_storage[[i]] <- btts$successes_storage
    trials_storage[[i]] <- btts$trials_storage
    mean_storage[[i]] <- btts$successes_storage / btts$trials_storage
    set_storage[[i]] <- (btts$successes_storage/btts$trials_storage) > threshold
    
    loss <- function(x) {return(1-(sum(x[1:4]) == 0 & sum(x[7:10]) == 4))}
    
    require(plyr)
    loss_storage[,i] <- as.vector(aaply(set_storage[[i]], 1, loss))
    
    message("Round ", i, " finished.")
  }
  
  return(list(successes_storage = successes_storage,
         trials_storage = trials_storage,
         mean_storage = mean_storage,
         set_storage = set_storage,
         loss_storage = loss_storage))
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




# APT Bandit

# assume the data frame that is passed as argument data has a column
# called 'variant' which contains the name of the groups
# and a column 'obs' containing the values
get_group_means <- function(data, ...) {
  result <- data %>%
    group_by(variant) %>%
    dplyr::summarize(mean = mean(obs)) %>%
    ungroup() %>%
    dplyr::arrange(variant) %>%
    #mutate(gap = abs(mean - tau) + eps,
    #       APT = sqrt(pulls_so_far)) %>% 
    as.data.frame()
  return(result$mean)
}

run_APT_bandit <- function(data, tau, eps, rounds = 500, seed = NA) {
  variants <- sort(unique(data[,"variant"]))
  k <- length(variants)
  storage <- data.frame(variant = variants, round = 1:k,
                        obs_index = rep(NA, k),
                        obs = rep(NA, k))  
  current_row <- 0
  
  # initialization of algorithm: get one observation for each variant
  for(i in 1:k) {
    next_row_info <- get_next_row(data = data, 
                                  variant = storage[i,"variant"],
                                  current_row = current_row)
    storage[i,"obs_index"] <- current_row <- next_row_info$next_row
    storage[i,"obs"] <- data[current_row, "success"]
  } 
  
  variant_summary <- data.frame(
    variant = variants,
    pulls_so_far = rep(1,3),
    mean = get_group_means(storage),
    gap = abs(get_group_means(storage)-tau)+eps,
    quantity = 1*(abs(get_group_means(storage)-tau)+eps),
    stringsAsFactors = FALSE)
  
  summary_over_time <- list(variant_summary)
  if(is.na(seed)) set.seed(seed)
  for(i in 1:(rounds-k)) {
    min_quantity <- min(summary_over_time[[i]]$quantity)
    arg_min_quantity <- which(summary_over_time[[i]]$quantity == min_quantity)
    if(length(arg_min_quantity) > 1) {
      arg_min_quantity <- sample(arg_min_quantity, 1)
    }
    
    next_row_info <- get_next_row(data = data, 
                                  variant = variants[arg_min_quantity],
                                  current_row = current_row)
    if(next_row_info$end == TRUE) {
      break
    }
    
    storage[k+i, "obs_index"] <- current_row <- next_row_info$next_row
    storage[k+i, "round"] <- k+i
    storage[k+i, "variant"] <- variants[arg_min_quantity]
    storage[k+i, "obs"] <- data[current_row, "success"]
    pulls_so_far <- summary_over_time[[i]]$pulls_so_far
    pulls_so_far[arg_min_quantity] <- pulls_so_far[arg_min_quantity] + 1
    
    summary_over_time[[i+1]] <- data.frame(
      variant = variants,
      pulls_so_far = pulls_so_far,
      mean = get_group_means(storage),
      gap = abs(get_group_means(storage)-tau)+eps,
      quantity = sqrt(pulls_so_far)*(abs(get_group_means(storage)-tau)+eps),
      stringsAsFactors = FALSE
    )
  }
  
  # Which variants have mean estimated to be larger than threshold tau?
  return_set <- which(summary_over_time[[length(summary_over_time)]]$mean > tau)
  if(length(return_set) == 0) {
    return_set <- NULL
  } else {
    return_set <- variants[return_set]
  }
  
  return(list(
    return_set = return_set,
    summary = summary_over_time[[length(summary_over_time)]],
    summary_over_time = summary_over_time,
    data_effective = storage,
    data = data
  ))
}

get_loss_at_each_time <- function(summary_over_time, tau, true_set) {
  rounds <- length(summary_over_time)
  compare_return_set <- function(df, tau, true_set) {
    !setequal(df[df$mean>tau,"variant"], true_set)
  }
  loss_over_time <- plyr::ldply(summary_over_time, compare_return_set,
                                tau, true_set)
  names(loss_over_time) <- "loss"
  return(loss_over_time)
}

simulate_bandit <- function(reps, rounds, tau, eps, true_set,
                            true_means, K,
                            seed = NA) {
  # parallelize this!
  loss_df <- data.frame(round = 1:rounds)
  summaries <- list()
  effective_data <- list()
  for(r in 1:reps) {
    data <- simulate_data_bernoulli(variants = K,
                                    means = true_means, 
                                    rounds = 50000,
                                    equal_splits = TRUE,
                                    splits = NULL,
                                    seed = seed+r)$experiment
    current_bandit <- run_APT_bandit(data = data, rounds = rounds, 
                                     tau = tau, eps = eps, seed = seed+r)
    summaries[[r]] <- current_bandit$summary
    effective_data[[r]] <- current_bandit$data_effective
    current_loss <- get_loss_at_each_time(current_bandit$summary_over_time,
                                          tau = tau, true_set = true_set)
    if(dim(current_loss)[1] < rounds) {
      current_loss <- c(current_loss$loss, rep(NA, 
                                          times = rounds-dim(current_loss)[1]))
      loss_df[,r+1] <- current_loss
    } else {
      loss_df[,r+1] <- current_loss$loss
    }
    
    names(loss_df)[r+1] <- paste("sim", r, sep = "")
    message("Simulation ", r, " done!")
  }
  
  return(list(
    summaries = summaries,
    effective_data = effective_data,
    loss = na.omit(loss_df)))
}

sim_res <- simulate_bandit(reps = 10, rounds = 10000, tau = 0.05, eps = 0.005,
                true_set = c("B", "C"), seed = 512,
                true_means = c(0.03, 0.07, 0.09), K = 3)


plot(1:dim(sim_res$loss)[1], rowMeans(sim_res$loss[,-1]), type = "l")
plot(1:dim(sim_res$loss)[1], log(rowMeans(sim_res$loss[,-1])), type = "l")
sim_res$summaries

