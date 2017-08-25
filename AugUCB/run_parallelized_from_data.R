# parallelize bandit simulation

#set.seed(512)
#testt <- data.frame(rep(NA, times = 6000))
#for(i in 1:length(mean_002)) {
#  testt[[i]] <- as.numeric(purrr::rbernoulli(6000, p  = mean_002[i]))
#}

#set.seed(513)
#testt2 <- data.frame(rep(NA, times = 6000))
#for(i in 1:length(mean_002)) {
#  testt2[[i]] <- as.numeric(purrr::rbernoulli(6000, p  = mean_002[i]))
#}

#test_list <- list(testt, testt2)

#para_test_res <- para_bandit_sim_APT(data = test_list, rounds = 5000, 
#                                     tau = tau_002, epsilon = epsilon_002)
#para_test_resKL <- para_bandit_sim_KL(data = test_list, rounds = 5000, 
#                                      tau = tau_002, epsilon = epsilon_002, 
#                                      at_tau = FALSE)
#para_test_resuni <- para_bandit_sim_uniform(data = test_list, rounds = 5000)
#para_test_resPI <- para_bandit_sim_PI(data = test_list, rounds = 5000,
#                                      tau = tau_002, epsilon = epsilon_002,
#                                      alpha = tau_002, beta = 1 - tau_002)
#para_test_resbucb <- para_bandit_sim_bucb(data = test_list, rounds = 5000, 
#                                          rate = "inverse",
#                                          tau = tau_002, epsilon = epsilon_002, 
#                                          alpha = tau_002, beta = 1-tau_002)

#############################################################################

# Get the average performance over all simulations

# expect a list of matrices of size rounds x K as input from the simulations,
# as well as a vector with the true means of length K to compare agains
compare_to_ground_truth <- function(true_means, sim_res, tau, epsilon) {
  message(paste0("Comparing ", length(sim_res), " simulations."))
  true_classification_up <- which(true_means > tau+epsilon)
  true_classification_down <- which(true_means < tau-epsilon)
  
  get_iter_error <- function(x, true_class_up, true_class_down) {
    ifelse(sum(true_class_up %in% which(x >= tau)) == length(true_class_up) &
             sum(true_class_down %in% which(x < tau)) == length(true_class_down),
           0, 1)
  }
  
  get_error <- function(res_df, ...) {
    apply(res_df$mean_storage, 1, get_iter_error, ...)
  }
  
  comp_list <- lapply(sim_res, get_error, true_class_up = true_classification_up, 
                      true_class_down = true_classification_down)
  
  comp_mean <- rowMeans(as.data.frame(comp_list, 
                                      col.names = 1:length(comp_list)))
  
  return(list(mean = comp_mean, full = comp_list))
}

#############################################################################

# Instead of using the hard loss function from Locatelli et al. (2016),
# use a loss that's more similar to simple regret (deviation from true mean)
# for those that are classified wrong
get_simple_regret <- function(true_means, sim_res, tau, epsilon) {
  message(paste0("Comparing ", length(sim_res), " simulations."))
  true_classification_up <- which(true_means > tau+epsilon)
  true_classification_down <- which(true_means < tau-epsilon)
  
  get_iter_error <- function(x, true_class_up, true_class_down) {
    wrong_class_index <- (true_means >= tau & x < tau) | (true_means < tau & x >= tau)
    sum(abs(x-true_means)[wrong_class_index])
  }
  
  get_error <- function(res_df, ...) {
    apply(res_df$mean_storage, 1, get_iter_error, ...)
  }
  
  comp_list <- lapply(sim_res, get_error, 
                      true_class_up = true_classification_up, 
                      true_class_down = true_classification_down)
  
  comp_mean <- rowMeans(as.data.frame(comp_list, 
                                      col.names = 1:length(comp_list)))
  
  return(list(mean = comp_mean, full = comp_list))
}

#testt <- gex_APT
#res <- get_simple_regret(mu_gex, testt, tau_gex, epsilon_gex)

#set.seed(93468734)
#tau_exp2 <- 1
#epsilon_exp2 <- 0
#mu_exp2 <- rexp(20, 1)
#load(paste0(current_path, "exp2_EVT.Rda"))
#res <- get_simple_regret(mu_exp2, exp2_EVT, tau_exp2, epsilon_exp2)
#plot(res$mean)
#############################################################################

# Get the average performance over all simulations

# expect a list of matrices of size rounds x K as input from the simulations,
# as well as a vector with the true means of length K to compare agains
compare_to_cv_data <- function(mean_list, sim_res, thresh, eps) {
  
  get_iter_error <- function(x, true_class_up, true_class_down) {
    ifelse(sum(true_class_up %in% which(x >= thresh)) == length(true_class_up) &&
             sum(true_class_down %in% which(x < thresh)) == length(true_class_down),
           0,1)
  }
  
  iter_error_vectors <- list()
  
  for(i in 1:length(sim_res)) {
    true_classification_up <- which(mean_list[[i]] > thresh+eps)
    true_classification_down <- which(mean_list[[i]] < thresh-eps)
    
    iter_error_vectors[[i]] <- apply(sim_res[[i]]$mean_storage, 1, get_iter_error,
                                     true_class_up = true_classification_up,
                                     true_class_down = true_classification_down)
  }
  
  comp_mean <- rowMeans(as.data.frame(iter_error_vectors, 
                                      col.names = 1:length(iter_error_vectors)))
  
  return(list(mean = comp_mean, full = iter_error_vectors))
}

#################################################################################

# The result is a very long data frame where each row is a (round, arm)
# combination and specifies which arm was wrong in that round; since the
# data frame carries information on all iterations of the simulation,
# some combinations will (and should!) be repeated; this way it's perfect
# for further piping into ggplot to discover which arms are classified
# wrong by the algorithm at what stage of the sampling process

# REQUIRES plyr AND IS SLOW!

get_wrong_arms_per_round <- function(true_means, sim_res, tau, epsilon) {
  n <- length(true_means)
  message(paste0("Comparing ", length(sim_res), " simulations."))
  true_classification_up <- which(true_means > tau+epsilon)
  true_classification_down <- which(true_means < tau-epsilon)
  
  get_wrong_arms_in_iter <- function(x, true_class_up, true_class_down) {
    arms <- sort(c(true_class_up[!(true_class_up %in% which(x >= tau))],
                   true_class_down[!(true_class_down %in% which(x < tau))]))
    ifelse(length(arms) == 0, NA, arms)
  }
  
  get_wrong_arms <- function(res_df, ...) {
    larms <- alply(res_df$mean_storage, .margins = 1, get_wrong_arms_in_iter, ...,
                   .expand = TRUE)
    dfarms <- data.frame(round = rep(0, n), arm = 1:n)
    for(i in 1:length(larms)) {
      dfarms <- rbind(dfarms,
                      data.frame(round = rep(i, length(larms[[i]])),
                                 arm = larms[[i]]))
    }
    return(dfarms)
  }
  
  comp_list <- ldply(sim_res, get_wrong_arms, 
                     true_class_up = true_classification_up,
                     true_class_down = true_classification_down,
                     .progress = "text")
  
  return(comp_list)
}

#################################################################################
# input has to be a list of different data frames
# each data frame is a run of the algorithm
# in the end, return a list of lists
# each object of the first list should be a list that contains 
# the mean_storage data frame and the arm_sequence vector

para_bandit_sim_APT <- function(data, seed = NA, do_verbose = FALSE, ...) {
  require(foreach)
  require(doParallel)
  
  # assume that data is a list of data frames
  reps <- length(data)
  
  gc()
  cl <- makeCluster(max(1,detectCores()-1))
  registerDoParallel(cl)
  res <- foreach(j = 1:reps, #.errorhandling = 'remove',
          .export = c("APT_from_tsdata", "get_next_arm_apt",
                      "get_min"), 
          .verbose = do_verbose, .inorder = TRUE) %dopar% {
            alg_res <- APT_from_tsdata(data = data[[j]], 
                                       seed = 512+j, ...)
            list(mean_storage = alg_res$mean_storage,
                 arm_sequence = alg_res$arm_sequence,
                 input_data = data[[j]])
            }
  stopCluster(cl)
  return(res)
}

#################################################################################
# EVT Algorithm

para_bandit_sim_EVT <- function(data, seed = NA, do_verbose = FALSE, ...) {
  require(foreach)
  require(doParallel)
  
  # assume that data is a list of data frames
  reps <- length(data)
  
  gc()
  cl <- makeCluster(max(1,detectCores()-1))
  registerDoParallel(cl)
  res <- foreach(j = 1:reps, #.errorhandling = 'remove',
                 .export = c("EVT_from_tsdata", "get_next_arm_evt",
                             "get_min"), 
                 .verbose = do_verbose, .inorder = TRUE) %dopar% {
                   alg_res <- EVT_from_tsdata(data = data[[j]], 
                                              seed = 512+j, ...)
                   list(mean_storage = alg_res$mean_storage,
                        arm_sequence = alg_res$arm_sequence,
                        input_data = data[[j]])
                 }
  stopCluster(cl)
  return(res)
}

#################################################################################
# input has to be a list of different data frames
# each data frame is a run of the algorithm
# in the end, return a list of lists
# each object of the first list should be a list that contains 
# the mean_storage data frame and the arm_sequence vector

para_bandit_sim_LR <- function(data, seed = NA, do_verbose = FALSE, ...) {
  require(foreach)
  require(doParallel)
  
  # assume that data is a list of data frames
  reps <- length(data)
  
  gc()
  cl <- makeCluster(max(1,detectCores()-1))
  registerDoParallel(cl)
  res <- foreach(j = 1:reps, #.errorhandling = 'remove',
                 .export = c("LR_bandit_from_tsdata", "get_next_arm_lr",
                             "get_min", "lr_ber"), 
                 .verbose = do_verbose, .inorder = TRUE) %dopar% {
                   alg_res <- LR_bandit_from_tsdata(data = data[[j]], 
                                              seed = 512+j, ...)
                   list(mean_storage = alg_res$mean_storage,
                        arm_sequence = alg_res$arm_sequence,
                        input_data = data[[j]])
                 }
  stopCluster(cl)
  return(res)
}

para_bandit_sim_LR_gaussian <- function(data, seed = NA, do_verbose = FALSE, ...) {
  require(foreach)
  require(doParallel)
  
  # assume that data is a list of data frames
  reps <- length(data)
  
  gc()
  cl <- makeCluster(max(1,detectCores()-1))
  registerDoParallel(cl)
  res <- foreach(j = 1:reps, #.errorhandling = 'remove',
                 .export = c("LR_bandit_from_tsdata_gaussian", "get_min", 
                             "kl_gaussian", "get_next_arm_kl_gaussian"), 
                 .verbose = do_verbose, .inorder = TRUE) %dopar% {
                   alg_res <- LR_bandit_from_tsdata_gaussian(data = data[[j]], 
                                                    seed = 512+j, ...)
                   list(mean_storage = alg_res$mean_storage,
                        arm_sequence = alg_res$arm_sequence,
                        input_data = data[[j]])
                 }
  stopCluster(cl)
  return(res)
}

para_bandit_sim_LR_exponential <- function(data, seed = NA, do_verbose = FALSE, ...) {
  require(foreach)
  require(doParallel)
  
  # assume that data is a list of data frames
  reps <- length(data)
  
  gc()
  cl <- makeCluster(max(1,detectCores()-1))
  registerDoParallel(cl)
  res <- foreach(j = 1:reps, #.errorhandling = 'remove',
                 .export = c("LR_bandit_from_tsdata_exponential", "get_min", 
                             "kl_exponential", "get_next_arm_kl_exponential"), 
                 .verbose = do_verbose, .inorder = TRUE) %dopar% {
                   alg_res <- LR_bandit_from_tsdata_exponential(data = data[[j]], 
                                                             seed = 512+j, ...)
                   list(mean_storage = alg_res$mean_storage,
                        arm_sequence = alg_res$arm_sequence,
                        input_data = data[[j]])
                 }
  stopCluster(cl)
  return(res)
}

para_bandit_sim_LR_poisson <- function(data, seed = NA, do_verbose = FALSE, ...) {
  require(foreach)
  require(doParallel)
  
  # assume that data is a list of data frames
  reps <- length(data)
  
  gc()
  cl <- makeCluster(max(1,detectCores()-1))
  registerDoParallel(cl)
  res <- foreach(j = 1:reps, #.errorhandling = 'remove',
                 .export = c("LR_bandit_from_tsdata_poisson", "get_min", 
                             "kl_poisson", "get_next_arm_kl_poisson"), 
                 .verbose = do_verbose, .inorder = TRUE) %dopar% {
                   alg_res <- LR_bandit_from_tsdata_poisson(data = data[[j]], 
                                                            seed = 512+j, ...)
                   list(mean_storage = alg_res$mean_storage,
                        arm_sequence = alg_res$arm_sequence,
                        input_data = data[[j]])
                 }
  stopCluster(cl)
  return(res)
}

#################################################################################

para_bandit_sim_AugUCB <- function(data, seed = NA, do_verbose = FALSE, ...) {
  require(foreach)
  require(doParallel)
  
  # assume that data is a list of data frames
  reps <- length(data)
  
  gc()
  cl <- makeCluster(max(1,detectCores()-1))
  registerDoParallel(cl)
  res <- foreach(j = 1:reps, #.errorhandling = 'remove',
                 .export = c("AugUCB_from_tsdata", "get_next_arm_augucb",
                             "get_min", "delete_arms"), 
                 .verbose = do_verbose, .inorder = TRUE) %dopar% {
                   alg_res <- AugUCB_from_tsdata(data = data[[j]], 
                                                 seed = 512+j, ...)
                   list(mean_storage = alg_res$mean_storage,
                        arm_sequence = alg_res$arm_sequence,
                        input_data = data[[j]])
                 }
  stopCluster(cl)
  return(res)
}

#################################################################################

para_bandit_sim_KL <- function(data, seed = NA, do_verbose = FALSE, ...) {
  require(foreach)
  require(doParallel)
  
  # assume that data is a list of data frames
  reps <- length(data)
  
  gc()
  cl <- makeCluster(max(1,detectCores()-1))
  registerDoParallel(cl)
  res <- foreach(j = 1:reps, .errorhandling = 'remove',
          .export = c("KL_bandit_from_tsdata", "get_next_arm_kl",
                      "get_min", "kl_ber", "get_next_arm_kl_at_tau"), 
          .verbose = do_verbose, .inorder = TRUE) %dopar% {
                        alg_res <- KL_bandit_from_tsdata(data = data[[j]], 
                                                         seed = 512+j, ...)
                        list(mean_storage = alg_res$mean_storage,
                             arm_sequence = alg_res$arm_sequence,
                             input_data = data[[j]])
          }
  stopCluster(cl)
  return(res)
}

#################################################################################

para_bandit_sim_KLUCB <- function(data, seed = NA, do_verbose = FALSE, ...) {
  require(foreach)
  require(doParallel)
  
  # assume that data is a list of data frames
  reps <- length(data)
  
  gc()
  cl <- makeCluster(max(1,detectCores()-1))
  registerDoParallel(cl)
  res <- foreach(j = 1:reps, .errorhandling = 'remove',
          .export = c("KLUCB_bandit_from_tsdata", "get_next_arm_klucb",
                      "get_min", "kl_ber", "get_klub", "get_kllb"), 
          .verbose = do_verbose, .inorder = TRUE) %dopar% {
            alg_res <- KLUCB_bandit_from_tsdata(data = data[[j]], 
                                                seed = 512+j, ...)
            list(mean_storage = alg_res$mean_storage,
                 arm_sequence = alg_res$arm_sequence,
                 input_data = data[[j]])
          }
  stopCluster(cl)
  return(res)
}

#################################################################################

para_bandit_sim_uniform <- function(data, seed = NA, do_verbose = FALSE, ...) {
  require(foreach)
  require(doParallel)
  
  # assume that data is a list of data frames
  reps <- length(data)
  
  gc()
  cl <- makeCluster(max(1,detectCores()-1))
  registerDoParallel(cl)
  res <- foreach(j = 1:reps, .errorhandling = 'remove',
          .export = c("uniform_bandit_from_tsdata",
                      "get_next_arm_uniform", "get_min"), 
          .verbose = do_verbose, .inorder = TRUE) %dopar% {
            alg_res <- uniform_bandit_from_tsdata(data = data[[j]], 
                                                  #seed = 512+j, 
                                                  ...)
            list(mean_storage = alg_res$mean_storage,
                 arm_sequence = alg_res$arm_sequence,
                 input_data = data[[j]],
                 iter = j)
          }
  stopCluster(cl)
  return(res)
}

#################################################################################

para_bandit_sim_PI <- function(data, seed = NA, do_verbose = FALSE, ...) {
  require(foreach)
  require(doParallel)
  
  # assume that data is a list of data frames
  reps <- length(data)
  
  gc()
  cl <- makeCluster(max(1,detectCores()-1))
  registerDoParallel(cl)
  foreach(j = 1:reps, .errorhandling = 'remove',
          .export = c("PI_bandit_from_tsdata", "get_next_arm_PI",
                      "get_min", "get_posterior_mean"), 
          .verbose = do_verbose, .inorder = TRUE) %dopar% {
            alg_res <- PI_bandit_from_tsdata(data = data[[j]], 
                                             seed = 512+j, ...)
            list(mean_storage = alg_res$mean_storage,
                 arm_sequence = alg_res$arm_sequence,
                 input_data = data[[j]])
          }
}

#############################################################################

para_bandit_sim_TTS <- function(data, seed = NA, do_verbose = FALSE, ...) {
  require(foreach)
  require(doParallel)
  
  # assume that data is a list of data frames
  reps <- length(data)
  
  gc()
  cl <- makeCluster(max(1,detectCores()-1))
  registerDoParallel(cl)
  foreach(j = 1:reps, #.errorhandling = 'remove',
          .export = c("TTS_from_tsdata", "get_next_arm_PI",
                      "draw_arm_tts", "get_posterior_mean"),
          .verbose = do_verbose, .inorder = TRUE) %dopar% {
            alg_res <- TTS_from_tsdata(data = data[[j]], 
                                       seed = 512+j, ...)
            list(mean_storage = alg_res$mean_storage,
                 arm_sequence = alg_res$arm_sequence,
                 input_data = data[[j]])
          }
}

#############################################################################

para_bandit_sim_BETA <- function(data, seed = NA, do_verbose = FALSE, ...) {
  require(foreach)
  require(doParallel)
  
  # assume that data is a list of data frames
  reps <- length(data)
  
  gc()
  cl <- makeCluster(max(1,detectCores()-1))
  registerDoParallel(cl)
  foreach(j = 1:reps, #.errorhandling = 'remove',
          .export = c("BETA_from_tsdata", "get_next_arm_BETA",
                      "get_posterior_mean", "get_min"),
          .verbose = do_verbose, .inorder = TRUE) %dopar% {
            alg_res <- BETA_from_tsdata(data = data[[j]], 
                                        seed = 512+j, ...)
            list(mean_storage = alg_res$mean_storage,
                 arm_sequence = alg_res$arm_sequence,
                 input_data = data[[j]])
          }
}

#############################################################################

para_bandit_sim_bucb <- function(data, seed = NA, do_verbose = FALSE, ...) {
  require(foreach)
  require(doParallel)
  
  # assume that data is a list of data frames
  reps <- length(data)
  
  gc()
  cl <- makeCluster(max(1,detectCores()-1))
  registerDoParallel(cl)
  res <- foreach(j = 1:reps, .errorhandling = 'remove',
          .export = c("BayesUCB_from_tsdata",
                      "get_BayesUCB_metric", "get_posterior_mean", "get_min"),
          .verbose = do_verbose, .inorder = TRUE) %dopar% {
            alg_res <- BayesUCB_from_tsdata(data = data[[j]], 
                                            seed = 512+j, ...)
            list(mean_storage = alg_res$mean_storage,
                 arm_sequence = alg_res$arm_sequence,
                 input_data = data[[j]])
          }
  stopCluster(cl)
  return(res)
}

#############################################################################