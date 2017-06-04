# parallelize bandit simulation

# what is the output that we need?
# For a normal simulation based on simulated data drawn on the go:
# mean estimate at every round for every simulation

#AugUCB(means, K = 4, rounds = 5000, rho = 1/3, tau = 0.5,
#       variances = NA, verbose = FALSE, seed = NA)
# currently returns a list where the object $mean_storage is a rounds x K matrix
# So we need to return a list of these matrices; one matrix for each iteration

para_bandit_sim_AugUCB <- function(reps = 6, seed = NA, do_verbose = FALSE, ...) {
  require(foreach)
  require(doParallel)
  
  gc()
  cl <- makeCluster(max(1,detectCores()-1))
  registerDoParallel(cl)
  foreach(j = 1:reps, .errorhandling = 'remove',
          .export = c("AugUCB", "pull_arm", "get_next_arm_augucb", "delete_arms",
                      "get_min"), .verbose = do_verbose) %dopar% 
    AugUCB(...)$mean_storage
}

means <- c(0.55, 0.38, 0.95)
variances <- c(0.5, 0.9, 1.7)
system.time(
para_res <- para_bandit_sim_AugUCB(reps = 100, means = means, variances = variances, 
                                   K = 3, rounds = 5000)
)
#  user  system elapsed 
# 0.331   0.256  98.026

para_bandit_sim_APT <- function(reps = 6, seed = NA, do_verbose = FALSE, ...) {
  require(foreach)
  require(doParallel)
  
  gc()
  cl <- makeCluster(max(1,detectCores()-1))
  registerDoParallel(cl)
  foreach(j = 1:reps, .errorhandling = 'remove',
          .export = c("APT", "pull_arm", "get_next_arm_apt",
                      "get_min"), .verbose = do_verbose) %dopar% 
    APT(...)$mean_storage
}

system.time(
  para_apt <- para_bandit_sim_APT(reps = 100, means = means, variances = variances, 
                                  K = 3, rounds = 5000)
)
#  user  system elapsed 
# 0.255   0.126  33.540

para_bandit_sim_uniform <- function(reps = 6, seed = NA, do_verbose = FALSE, ...) {
  require(foreach)
  require(doParallel)
  
  gc()
  cl <- makeCluster(max(1,detectCores()-1))
  registerDoParallel(cl)
  foreach(j = 1:reps, .errorhandling = 'remove',
          .export = c("uniform_bandit", "pull_arm", "get_next_arm_uniform",
                      "get_min"), .verbose = do_verbose) %dopar% 
    uniform_bandit(...)$mean_storage
}

system.time(
  para_uniform <- para_bandit_sim_uniform(reps = 100, means = means, 
                                          variances = variances, K = 3, rounds = 5000)
)
#  user  system elapsed 
# 0.241   0.098  28.157

#############################################################################

# Get the average performance over all simulations

# expect a list of matrices of size rounds x K as input from the simulations,
# as well as a vector with the true means of length K to compare agains
compare_to_ground_truth <- function(true_means, sim_res, tau, epsilon) {
  message(paste0("Comparing ", length(sim_res), " simulations."))
  true_classification_up <- which(true_means > tau+epsilon)
  true_classification_down <- which(true_means < tau-epsilon)
  
  get_iter_error <- function(x, true_class_up, true_class_down) {
    ifelse(true_class_up %in% which(x > tau) && 
             true_class_down %in% which(x < tau),
           0,1)
  }
  
  get_error <- function(res_df, ...) {
    apply(res_df, 1, get_iter_error, ...)
  }
  
  comp_list <- lapply(sim_res, get_error, true_class_up = true_classification_up, 
                      true_class_down = true_classification_down)
  
  comp_mean <- rowMeans(as.data.frame(comp_list, 
                                      col.names = 1:length(comp_list)))
  
  return(list(mean = comp_mean, full = comp_list))
}

#compared_res <- compare_to_ground_truth(means, para_res, 0.5, 0.1)
#plot(compared_res$mean, type = "l")

#compared_apt <- compare_to_ground_truth(means, para_apt, 0.5, 0.1)
#plot(compared_apt$mean, type = "l")
#lines(compared_res$mean, col = "blue")

#compared_uniform <- compare_to_ground_truth(means, para_uniform, 0.5, 0.1)
#plot(log(compared_uniform$mean), type = "l", xlim = c(0,1500))
#lines(log(compared_res$mean), col = "blue")
#lines(log(compared_apt$mean), col = "red")
