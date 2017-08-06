########################################################################
# load functions

current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/AugUCB/"
source(paste0(current_path, "run_from_data.R"))
source(paste0(current_path, "run_parallelized_from_data.R"))
current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/synthetic_experiments/gaussian/"

########################################################################
# Create the data

tau_gax1 <- 10
epsilon_gax1 <- 0
mu_gax1 <- c(-5, -5,
             0, 0,
             5, 5,
             13, 15,
             18, 20)
sd_gax1 <- c(5, 10,
             5, 10,
             5, 10,
             5, 10,
             5, 10)

data_list_gax1 <- list()
set.seed(268549023)
for(j in 1:1000) {
  curr_data <- data.frame(rep(NA, times = 2000))
  for(i in 1:length(mu_gax1)) {
    curr_data[[i]] <- as.numeric(rnorm(2000, mu_gax1[i], sd_gax1[i]))
  }
  names(curr_data) <- paste0("V", rep(1:length(mu_gax1)))
  data_list_gax1[[j]] <- curr_data
}
gc()
########################################################################

plot(c(0,2000), c(0, -8), type = "n")
lines(log(gax1_comp_APT), col = "black")
lines(log(gax1_comp_LR), col = "red")

########################################################################
# Gaussian Likelihood Ratio
single_gax1_LR <- LR_bandit_from_tsdata_gaussian(data_list_gax1[[1]], 
                                        rounds = 2000, tau_gax1, epsilon_gax1,
                                        verbose = TRUE, seed = NA)

gax1_LR <- list()
for(j in 1:length(data_list_gax1)) {
  message(j)
  alg_res <- LR_bandit_from_tsdata_gaussian(data_list_gax1[[j]], 
                                 rounds = 2000, tau_gax1, epsilon_gax1,
                                 verbose = TRUE, seed = 512+j)
  gax1_LR[[j]] <- list(mean_storage = alg_res$mean_storage,
                       arm_sequence = alg_res$arm_sequence,
                       input_data = data_list_gax1[[j]])
}

system.time(gax1_LR <- para_bandit_sim_LR_gaussian(data = data_list_gax1, 
                                                   rounds = 2000, 
                                                   tau = tau_gax1, 
                                                   epsilon = epsilon_gax1,
                                                   do_verbose = TRUE))



save(gax1_LR, file = paste0(current_path, "gax1_LR.Rda"))
gax1_comp_LR <- compare_to_ground_truth(mu_gax1, gax1_LR, 
                                        tau_gax1, epsilon_gax1)$mean
save(gax1_comp_LR, file = paste0(current_path, "gax1_comp_LR.Rda"))
gc()

########################################################################

system.time(gax1_APT <- para_bandit_sim_APT(data = data_list_gax1, 
                                            rounds = 2000, 
                                            tau = tau_gax1, 
                                            epsilon = epsilon_gax1,
                                            do_verbose = TRUE))
save(gax1_APT, file = paste0(current_path, "gax1_APT.Rda"))
gax1_comp_APT <- compare_to_ground_truth(mu_gax1, gax1_APT, 
                                        tau_gax1, epsilon_gax1)$mean
save(gax1_comp_APT, file = paste0(current_path, "gax1_comp_APT.Rda"))
gc()
