# run a simulation based on a simple example that can hopefully serve to
# calibrate the algorithms
########################################################################
# load functions

current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/AugUCB/"
source(paste0(current_path, "run_from_data.R"))
source(paste0(current_path, "run_parallelized_from_data.R"))
current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/synthetic_experiments/"

########################################################################
# Create the data

mean_loc2 <- c(0.4-0.2^(1:4),
               0.45, 0.55,
               0.6+0.2^(5-(1:4)))
tau_loc2 <- 0.5
epsilon_loc2 <- 0.1
H_loc2 <- get_complexity(mean_loc2, tau_loc2, epsilon_loc2)

plot(mean_loc2, rep(1,10), main = paste0("Complexity of ", round(H_loc2,2)))
abline(v=tau_loc2)
abline(v=tau_loc2+epsilon_loc2, lty=2)
abline(v=tau_loc2-epsilon_loc2, lty=2)

data_list2 <- list()
set.seed(1024)
for(j in 1:5000) {
  curr_data <- data.frame(rep(NA, times = 4000))
  for(i in 1:length(mean_loc2)) {
    curr_data[[i]] <- as.numeric(purrr::rbernoulli(4000, p  = mean_loc2[i]))
  }
  names(curr_data) <- paste0("V", rep(1:length(mean_loc2)))
  data_list2[[j]] <- curr_data
}

########################################################################
# Plot the results

load(paste0(current_path, "/loc2_comp_BUCB_horizon_4000.Rda"))
load(paste0(current_path, "/loc2_comp_APT_2000.Rda"))
load(paste0(current_path, "/loc2_comp_AugUCB_2000.Rda"))

plot(c(0,2000), c(0, -9), type = "n")
lines(log(loc2_comp_BUCB_horizon_4000), col = "blue")
lines(log(loc2_comp_APT_2000), col = "red")
lines(log(loc2_comp_AugUCB_2000), col = "grey")
#lines(log(loc2_comp_KLUCB100), col = "green")
#lines(log(loc2_comp_KLUCB80), col = "green")
lines(log(loc2_comp_KLUCB80_long), col = "green")
lines(log(loc2_comp_KLUCB15_long), col = "darkgreen")
#lines(log(loc2_comp_KLUCB60), col = "darkgreen")
#lines(log(loc2_comp_KLUCB20), col = "darkgreen")

8000/H_loc2

########################################################################

loc2_BUCB_horizon_4000 <- para_bandit_sim_bucb(data = data_list2, rounds = 4000, 
                                          rate = "inverse_horizon",
                                          tau = tau_loc2, epsilon = epsilon_loc2, 
                                          alpha = tau_loc2, beta = 1-tau_loc2)
save(loc2_BUCB_horizon_4000, 
     file = paste0(current_path, "loc2_BUCB_horizon_4000.Rda"))
loc2_comp_BUCB_horizon_4000 <- compare_to_ground_truth(mean_loc2, 
                                                  loc2_BUCB_horizon_4000,
                                                  tau_loc2,
                                                  epsilon_loc2)$mean
save(loc2_comp_BUCB_horizon_4000, file = paste0(current_path, "loc2_comp_BUCB_horizon_4000.Rda"))

########################################################################
# Standard APT Algorithm
# a seed is built in
system.time(loc2_APT_2000 <- para_bandit_sim_APT(data = data_list2, 
                                                 rounds = 2000, 
                                                 tau = tau_loc2, 
                                                 epsilon = epsilon_loc2))
save(loc2_APT_2000, file = paste0(current_path, "loc2_APT_2000.Rda"))
#load(file = paste0(current_path, "loc2_APT.Rda"))
loc2_comp_APT_2000 <- compare_to_ground_truth(mean_loc2, loc2_APT_2000, tau_loc2, 
                                         epsilon_loc2)$mean
save(loc2_comp_APT_2000, file = paste0(current_path, "loc2_comp_APT_2000.Rda"))
rm(loc2_APT)
gc()

########################################################################
# Standard AugUCB
system.time(loc2_AugUCB_2000 <- para_bandit_sim_AugUCB(data = data_list2, 
                                                  rounds = 2000, 
                                                  tau = tau_loc2))
# user   system  elapsed 
# 13.972   10.042 2284.875
save(loc2_AugUCB_2000, file = paste0(current_path, "loc2_AugUCB_2000.Rda"))
#load(file = paste0(current_path, "loc2_AugUCB.Rda"))
loc2_comp_AugUCB_2000 <- compare_to_ground_truth(mean_loc2, loc2_AugUCB_2000,
                                                 tau_loc2, 
                                                 epsilon_loc2)$mean
save(loc2_comp_AugUCB_2000, file = paste0(current_path,
                                          "loc2_comp_AugUCB_2000.Rda"))
rm(loc2_AugUCB_2000)
gc()

########################################################################
# KL-UCB

error_bound <- function(n, H, K) {
  (K*n+H)*exp(-n/H)
}

plot(1:10000, log(error_bound(1:10000, H_loc2, 10)), type = "l",
     ylim = c(-6,12))
abline(h=log(0.2))

data_list2_small <- data_list2[1:1000]
rm(data_list2)
gc()

system.time(loc2_KLUCB100 <- para_bandit_sim_KLUCB(data = data_list2_small[1:100], 
                                                  rounds = 1000, 
                                                  tau = tau_loc2, 
                                                  epsilon = epsilon_loc2,
                                                  horizon = 10000,
                                                  H = H_loc2))


system.time(loc2_KLUCB80 <- para_bandit_sim_KLUCB(data = data_list2_small, 
                                                   rounds = 1500, 
                                                   tau = tau_loc2, 
                                                   epsilon = epsilon_loc2,
                                                   horizon = 8000,
                                                   H = H_loc2))

system.time(loc2_KLUCB80_long <- para_bandit_sim_KLUCB(data = data_list2, 
                                                  rounds = 1500, 
                                                  tau = tau_loc2, 
                                                  epsilon = epsilon_loc2,
                                                  horizon = 8000,
                                                  H = H_loc2))
save(loc2_KLUCB80_long, file = paste0(current_path, "/loc2_KLUCB80_long.Rda"))
#    user   system  elapsed 
# 16.337   13.441 3564.603
system.time(loc2_KLUCB15_long <- para_bandit_sim_KLUCB(data = data_list2, 
                                                       rounds = 1500, 
                                                       tau = tau_loc2, 
                                                       epsilon = epsilon_loc2,
                                                       horizon = 1500,
                                                       H = H_loc2))
save(loc2_KLUCB15_long, file = paste0(current_path, "/loc2_KLUCB15_long.Rda"))


system.time(loc2_KLUCB60 <- para_bandit_sim_KLUCB(data = data_list2_small[1:100], 
                                                   rounds = 1000, 
                                                   tau = tau_loc2, 
                                                   epsilon = epsilon_loc2,
                                                   horizon = 6000,
                                                   H = H_loc2))

system.time(loc2_KLUCB20 <- para_bandit_sim_KLUCB(data = data_list2_small, 
                                                  rounds = 1500, 
                                                  tau = tau_loc2, 
                                                  epsilon = epsilon_loc2,
                                                  horizon = 2000,
                                                  H = H_loc2))
# user  system elapsed 
#51.048    73.114 34050.023
save(loc2_KLUCB100, file = paste0(current_path, "loc2_KLUCB100.Rda"))
#load(file = paste0(current_path, "loc4_KLUCB40.Rda"))
loc2_comp_KLUCB80 <- compare_to_ground_truth(mean_loc2, loc2_KLUCB80, 
                                              tau_loc2, 
                                              epsilon_loc2)$mean
loc2_comp_KLUCB80_long <- compare_to_ground_truth(mean_loc2, loc2_KLUCB80_long, 
                                             tau_loc2, 
                                             epsilon_loc2)$mean
loc2_comp_KLUCB15_long <- compare_to_ground_truth(mean_loc2, loc2_KLUCB15_long, 
                                                  tau_loc2, 
                                                  epsilon_loc2)$mean
loc2_comp_KLUCB60 <- compare_to_ground_truth(mean_loc2, loc2_KLUCB60, 
                                             tau_loc2, 
                                             epsilon_loc2)$mean
save(loc2_comp_KLUCB60, file = paste0(current_path, "loc2_comp_KLUCB60.Rda"))
loc2_comp_KLUCB20 <- compare_to_ground_truth(mean_loc2, loc2_KLUCB20, 
                                             tau_loc2, 
                                             epsilon_loc2)$mean
save(loc2_comp_KLUCB20, file = paste0(current_path, "loc2_comp_KLUCB20.Rda"))
save(loc2_comp_KLUCB80, file = paste0(current_path, "loc2_comp_KLUCB80.Rda"))
save(loc2_comp_KLUCB80_long, file = paste0(current_path, "loc2_comp_KLUCB80_long.Rda"))
save(loc2_comp_KLUCB15_long, file = paste0(current_path, "loc2_comp_KLUCB15_long.Rda"))
save(loc2_comp_KLUCB100, file = paste0(current_path, "loc2_comp_KLUCB100.Rda"))
#load(file = paste0(current_path, "loc4_comp_KLUCB40.Rda"))
rm(loc2_KLUCB90)
gc()