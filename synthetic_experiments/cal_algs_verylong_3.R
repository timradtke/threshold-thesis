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

mean_loc3 <- c(0.001, 0.005, 0.01, 0.015,
               0.0475, 0.0525,
               0.085, 0.09, 0.095, 0.099)
tau_loc3 <- 0.05
epsilon_loc3 <- 0.005
H_loc3 <- get_complexity(mean_loc3, tau_loc3, epsilon_loc3)

plot(mean_loc3, rep(1,10), main = paste0("Complexity of ", round(H_loc3,2)))
abline(v=tau_loc3)
abline(v=tau_loc3+epsilon_loc3, lty=2)
abline(v=tau_loc3-epsilon_loc3, lty=2)

data_list3 <- list()
set.seed(1024)
for(j in 1:5000) {
  curr_data <- data.frame(rep(NA, times = 7000))
  for(i in 1:length(mean_loc3)) {
    curr_data[[i]] <- as.numeric(purrr::rbernoulli(7000, p  = mean_loc3[i]))
  }
  names(curr_data) <- paste0("V", rep(1:length(mean_loc3)))
  data_list3[[j]] <- curr_data
}

########################################################################
# Plot the results

load(paste0(current_path, "/loc3_comp_BUCB_horizon_7000.Rda"))
load(paste0(current_path, "/loc3_comp_AugUCB_7000.Rda"))
load(paste0(current_path, "/loc3_comp_KLUCB370_5000.Rda"))
load(paste0(current_path, "/loc3_comp_UNIFORM_7000.Rda"))
load(paste0(current_path, "/loc3_comp_APT_7000.Rda"))
load(paste0(current_path, "/loc3_comp_LR_7000.Rda"))

plot(c(0,7000), c(0, -9), type = "n")
lines(log(loc3_comp_BUCB_horizon_7000), col = "blue")
#lines(log(loc3_comp_BUCB_horizon), col = "lightblue")
#lines(log(loc3_comp_APT), col = "pink")
lines(log(loc3_comp_LR_7000), col = "green")
lines(log(loc3_comp_APT_7000), col = "red")
lines(log(loc3_comp_UNIFORM_7000), col = "black")
#lines(log(loc3_comp_AugUCB), col = "grey")
lines(log(loc3_comp_AugUCB_7000), col = "grey")
#lines(log(loc3_comp_KLUCB7), col = "green")
#lines(log(loc3_comp_KLUCB30), col = "green")
#lines(log(loc3_comp_KLUCB370), col = "darkgreen")
lines(log(loc3_comp_KLUCB370_5000), col = "darkgreen")

########################################################################

loc3_BUCB_horizon_7000 <- para_bandit_sim_bucb(data = data_list3, 
                                               rounds = 7000, 
                                               rate = "inverse_horizon",
                                               tau = tau_loc3, 
                                               epsilon = epsilon_loc3, 
                                               alpha = tau_loc3, 
                                               beta = 1-tau_loc3)
save(loc3_BUCB_horizon_7000, 
     file = paste0(current_path, "loc3_BUCB_horizon_7000.Rda"))
loc3_comp_BUCB_horizon_7000 <- compare_to_ground_truth(mean_loc3, 
                                                       loc3_BUCB_horizon_7000,
                                                       tau_loc3,
                                                       epsilon_loc3)$mean
save(loc3_comp_BUCB_horizon_7000, file = paste0(current_path, 
                                                "loc3_comp_BUCB_horizon_7000.Rda"))

########################################################################
# Standard AugUCB
system.time(loc3_AugUCB_7000 <- para_bandit_sim_AugUCB(data = data_list3, 
                                                       rounds = 7000, 
                                                       tau = tau_loc3))
#     user    system   elapsed 
#46.915    50.807 18671.723
save(loc3_AugUCB_7000, file = paste0(current_path, "loc3_AugUCB_7000.Rda"))
loc3_comp_AugUCB_7000 <- compare_to_ground_truth(mean_loc3, loc3_AugUCB_7000, 
                                                 tau_loc3, epsilon_loc3)$mean
save(loc3_comp_AugUCB_7000, file = paste0(current_path, 
                                          "loc3_comp_AugUCB_7000.Rda"))

########################################################################
# Standard APT
system.time(loc3_APT_7000 <- para_bandit_sim_APT(data = data_list3, 
                                                 rounds = 7000, 
                                                 tau = tau_loc3, 
                                                 epsilon = epsilon_loc3))

save(loc3_APT_7000, file = paste0(current_path, "loc3_APT_7000.Rda"))
loc3_comp_APT_7000 <- compare_to_ground_truth(mean_loc3, loc3_APT_7000, 
                                                 tau_loc3, epsilon_loc3)$mean
save(loc3_comp_APT_7000, file = paste0(current_path, 
                                       "loc3_comp_APT_7000.Rda"))

########################################################################
# Standard Likelihood Ratio
system.time(loc3_LR_7000 <- para_bandit_sim_LR(data = data_list3, 
                                               rounds = 7000, 
                                               tau = tau_loc3, 
                                               epsilon = epsilon_loc3))

save(loc3_LR_7000, file = paste0(current_path, "loc3_LR_7000.Rda"))
loc3_comp_LR_7000 <- compare_to_ground_truth(mean_loc3, loc3_LR_7000, 
                                             tau_loc3, epsilon_loc3)$mean
save(loc3_comp_LR_7000, file = paste0(current_path, 
                                       "loc3_comp_LR_7000.Rda"))

########################################################################

system.time(loc3_UNIFORM_7000 <- para_bandit_sim_uniform(data = data_list3, 
                                                         rounds = 7000))
save(loc3_UNIFORM_7000, file = paste0(current_path, "loc3_UNIFORM_7000.Rda"))
loc3_comp_UNIFORM_7000 <- compare_to_ground_truth(mean_loc3, 
                                                  loc3_UNIFORM_7000, 
                                                  tau_loc3, 
                                                  epsilon_loc3)$mean
save(loc3_comp_UNIFORM_7000, file = paste0(current_path, 
                                           "loc3_comp_UNIFORM_7000.Rda"))

########################################################################

system.time(loc3_KLUCB370_5000 <- para_bandit_sim_KLUCB(data = data_list3, 
                                                   rounds = 5000, 
                                                   tau = tau_loc3, 
                                                   epsilon = epsilon_loc3,
                                                   horizon = 370000,
                                                   H = H_loc3))

save(loc3_KLUCB370_5000, file = paste0(current_path, 
                                          "loc3_KLUCB370_5000.Rda"))
loc3_comp_KLUCB370_5000 <- compare_to_ground_truth(mean_loc3, 
                                                   loc3_KLUCB370_5000, 
                                                   tau_loc3, epsilon_loc3)$mean
save(loc3_comp_KLUCB370_5000, file = paste0(current_path, 
                                          "loc3_comp_KLUCB370_5000.Rda"))

########################################################################

error_bound <- function(n, H, K) {
  (K*n+H)*exp(-n/H)
}

plot(1:500000, log(error_bound(1:500000, H_loc3, 10)), type = "l",
     ylim = c(-6,12))
abline(h=log(0.2))

data_list3_small <- data_list3[1:500]
rm(data_list3)
gc()

system.time(loc3_KLUCB30 <- para_bandit_sim_KLUCB(data = data_list3_small, 
                                                   rounds = 2000, 
                                                   tau = tau_loc3, 
                                                   epsilon = epsilon_loc3,
                                                   horizon = round(H_loc3),
                                                   H = H_loc3))

system.time(loc3_KLUCB370 <- para_bandit_sim_KLUCB(data = data_list3_small, 
                                                  rounds = 2000, 
                                                  tau = tau_loc3, 
                                                  epsilon = epsilon_loc3,
                                                  horizon = 370000,
                                                  H = H_loc3))

system.time(loc3_KLUCB7 <- para_bandit_sim_KLUCB(data = data_list3_small, 
                                                   rounds = 2000, 
                                                   tau = tau_loc3, 
                                                   epsilon = epsilon_loc3,
                                                   horizon = 7000,
                                                   H = H_loc3))

save(loc3_KLUCB30, file = paste0(current_path, "/loc3_KLUCB30.Rda"))
save(loc3_KLUCB370, file = paste0(current_path, "/loc3_KLUCB370.Rda"))
save(loc3_KLUCB7, file = paste0(current_path, "/loc3_KLUCB7.Rda"))

loc3_comp_KLUCB30 <- compare_to_ground_truth(mean_loc3, loc3_KLUCB30, 
                                             tau_loc3, 
                                             epsilon_loc3)$mean
loc3_comp_KLUCB370 <- compare_to_ground_truth(mean_loc3, loc3_KLUCB370, 
                                             tau_loc3, 
                                             epsilon_loc3)$mean
loc3_comp_KLUCB7 <- compare_to_ground_truth(mean_loc3, loc3_KLUCB7, 
                                             tau_loc3, 
                                             epsilon_loc3)$mean
save(loc3_comp_KLUCB30, file = paste0(current_path, "/loc3_comp_KLUCB30.Rda"))
save(loc3_comp_KLUCB370, file = paste0(current_path, "/loc3_comp_KLUCB370.Rda"))
save(loc3_comp_KLUCB7, file = paste0(current_path, "/loc3_comp_KLUCB7.Rda"))
