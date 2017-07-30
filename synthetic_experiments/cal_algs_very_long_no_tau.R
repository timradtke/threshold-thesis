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

mean_loc3nt <- c(0.001, 0.005, 0.01, 0.015,
               0.0475, 0.0525,
               0.085, 0.09, 0.095, 0.099)
tau_loc3nt <- 0.05
epsilon_loc3nt <- 0
H_loc3nt <- get_complexity(mean_loc3nt, tau_loc3nt, epsilon_loc3nt)

plot(mean_loc3nt, rep(1,10), main = paste0("Complexity of ", round(H_loc3nt,2)))
abline(v=tau_loc3nt)
abline(v=tau_loc3nt+epsilon_loc3nt, lty=2)
abline(v=tau_loc3nt-epsilon_loc3nt, lty=2)

data_list3 <- list()
set.seed(1024)
for(j in 1:5000) {
  curr_data <- data.frame(rep(NA, times = 7000))
  for(i in 1:length(mean_loc3nt)) {
    curr_data[[i]] <- as.numeric(purrr::rbernoulli(7000, p  = mean_loc3nt[i]))
  }
  names(curr_data) <- paste0("V", rep(1:length(mean_loc3nt)))
  data_list3[[j]] <- curr_data
}

########################################################################
# Plot the results

load(paste0(current_path, "/loc3nt_comp_BUCB_horizon_7000.Rda"))
load(paste0(current_path, "/loc3nt_comp_AugUCB_7000.Rda"))
load(paste0(current_path, "/loc3nt_comp_KLUCB370_5000.Rda"))
load(paste0(current_path, "/loc3nt_comp_UNIFORM_7000.Rda"))
load(paste0(current_path, "/loc3nt_comp_APT_7000.Rda"))
load(paste0(current_path, "/loc3nt_comp_LR_7000.Rda"))

plot(c(0,7000), c(0, -2), type = "n")
lines(log(loc3nt_comp_BUCB_horizon_7000), col = "blue")
#lines(log(loc3nt_comp_BUCB_horizon), col = "lightblue")
#lines(log(loc3nt_comp_APT), col = "pink")
lines(log(loc3nt_comp_LR_7000), col = "green")
lines(log(loc3nt_comp_APT_7000), col = "red")
lines(log(loc3nt_comp_UNIFORM_7000), col = "black")
#lines(log(loc3nt_comp_AugUCB), col = "grey")
lines(log(loc3nt_comp_AugUCB_7000), col = "grey")
#lines(log(loc3nt_comp_KLUCB7), col = "green")
#lines(log(loc3nt_comp_KLUCB30), col = "green")
#lines(log(loc3nt_comp_KLUCB370), col = "darkgreen")
lines(log(loc3nt_comp_KLUCB370_5000), col = "darkgreen")

########################################################################

loc3nt_BUCB_horizon_7000 <- para_bandit_sim_bucb(data = data_list3, 
                                               rounds = 7000, 
                                               rate = "inverse_horizon",
                                               tau = tau_loc3nt, 
                                               epsilon = epsilon_loc3nt, 
                                               alpha = tau_loc3nt, 
                                               beta = 1-tau_loc3nt)
save(loc3nt_BUCB_horizon_7000, 
     file = paste0(current_path, "loc3nt_BUCB_horizon_7000.Rda"))
loc3nt_comp_BUCB_horizon_7000 <- compare_to_ground_truth(mean_loc3nt, 
                                                       loc3nt_BUCB_horizon_7000,
                                                       tau_loc3nt,
                                                       epsilon_loc3nt)$mean
save(loc3nt_comp_BUCB_horizon_7000, file = paste0(current_path, 
                                                "loc3nt_comp_BUCB_horizon_7000.Rda"))

########################################################################
# Standard AugUCB
load(file = paste0(current_path, "loc3_AugUCB_7000.Rda"))
loc3nt_AugUCB_7000 <- loc3_AugUCB_7000
rm(loc3_AugUCB_7000)
gc()
save(loc3nt_AugUCB_7000, file = paste0(current_path, "loc3nt_AugUCB_7000.Rda"))
loc3nt_comp_AugUCB_7000 <- compare_to_ground_truth(mean_loc3nt, loc3nt_AugUCB_7000, 
                                                 tau_loc3nt, epsilon_loc3nt)$mean
save(loc3nt_comp_AugUCB_7000, file = paste0(current_path, 
                                          "loc3nt_comp_AugUCB_7000.Rda"))

########################################################################
# Standard APT
system.time(loc3nt_APT_7000 <- para_bandit_sim_APT(data = data_list3, 
                                                 rounds = 7000, 
                                                 tau = tau_loc3nt, 
                                                 epsilon = epsilon_loc3nt))

save(loc3nt_APT_7000, file = paste0(current_path, "loc3nt_APT_7000.Rda"))
loc3nt_comp_APT_7000 <- compare_to_ground_truth(mean_loc3nt, loc3nt_APT_7000, 
                                              tau_loc3nt, epsilon_loc3nt)$mean
save(loc3nt_comp_APT_7000, file = paste0(current_path, 
                                       "loc3nt_comp_APT_7000.Rda"))

########################################################################
# Standard Likelihood Ratio
system.time(loc3nt_LR_7000 <- para_bandit_sim_LR(data = data_list3, 
                                               rounds = 7000, 
                                               tau = tau_loc3nt, 
                                               epsilon = epsilon_loc3nt))

save(loc3nt_LR_7000, file = paste0(current_path, "loc3nt_LR_7000.Rda"))
loc3nt_comp_LR_7000 <- compare_to_ground_truth(mean_loc3nt, loc3nt_LR_7000, 
                                             tau_loc3nt, epsilon_loc3nt)$mean
save(loc3nt_comp_LR_7000, file = paste0(current_path, 
                                      "loc3nt_comp_LR_7000.Rda"))

########################################################################

load(paste0(current_path, "loc3_UNIFORM_7000.Rda"))
loc3nt_UNIFORM_7000 <- loc3_UNIFORM_7000
rm(loc3_UNIFORM_7000)
gc()
save(loc3nt_UNIFORM_7000, file = paste0(current_path, "loc3nt_UNIFORM_7000.Rda"))
loc3nt_comp_UNIFORM_7000 <- compare_to_ground_truth(mean_loc3nt, 
                                                  loc3nt_UNIFORM_7000, 
                                                  tau_loc3nt, 
                                                  epsilon_loc3nt)$mean
save(loc3nt_comp_UNIFORM_7000, file = paste0(current_path, 
                                           "loc3nt_comp_UNIFORM_7000.Rda"))

########################################################################

system.time(loc3nt_KLUCB370_5000 <- para_bandit_sim_KLUCB(data = data_list3, 
                                                        rounds = 5000, 
                                                        tau = tau_loc3nt, 
                                                        epsilon = epsilon_loc3nt,
                                                        horizon = 370000,
                                                        H = H_loc3nt))

save(loc3nt_KLUCB370_5000, file = paste0(current_path, 
                                       "loc3nt_KLUCB370_5000.Rda"))
loc3nt_comp_KLUCB370_5000 <- compare_to_ground_truth(mean_loc3nt, 
                                                   loc3nt_KLUCB370_5000, 
                                                   tau_loc3nt, epsilon_loc3nt)$mean
save(loc3nt_comp_KLUCB370_5000, file = paste0(current_path, 
                                            "loc3nt_comp_KLUCB370_5000.Rda"))

########################################################################

error_bound <- function(n, H, K) {
  (K*n+H)*exp(-n/H)
}

plot(1:500000, log(error_bound(1:500000, H_loc3nt, 10)), type = "l",
     ylim = c(-6,12))
abline(h=log(0.2))

data_list3_small <- data_list3[1:500]
rm(data_list3)
gc()

system.time(loc3nt_KLUCB30 <- para_bandit_sim_KLUCB(data = data_list3_small, 
                                                  rounds = 2000, 
                                                  tau = tau_loc3nt, 
                                                  epsilon = epsilon_loc3nt,
                                                  horizon = round(H_loc3nt),
                                                  H = H_loc3nt))

system.time(loc3nt_KLUCB370 <- para_bandit_sim_KLUCB(data = data_list3_small, 
                                                   rounds = 2000, 
                                                   tau = tau_loc3nt, 
                                                   epsilon = epsilon_loc3nt,
                                                   horizon = 370000,
                                                   H = H_loc3nt))

system.time(loc3nt_KLUCB7 <- para_bandit_sim_KLUCB(data = data_list3_small, 
                                                 rounds = 2000, 
                                                 tau = tau_loc3nt, 
                                                 epsilon = epsilon_loc3nt,
                                                 horizon = 7000,
                                                 H = H_loc3nt))

save(loc3nt_KLUCB30, file = paste0(current_path, "/loc3nt_KLUCB30.Rda"))
save(loc3nt_KLUCB370, file = paste0(current_path, "/loc3nt_KLUCB370.Rda"))
save(loc3nt_KLUCB7, file = paste0(current_path, "/loc3nt_KLUCB7.Rda"))

loc3nt_comp_KLUCB30 <- compare_to_ground_truth(mean_loc3nt, loc3nt_KLUCB30, 
                                             tau_loc3nt, 
                                             epsilon_loc3nt)$mean
loc3nt_comp_KLUCB370 <- compare_to_ground_truth(mean_loc3nt, loc3nt_KLUCB370, 
                                              tau_loc3nt, 
                                              epsilon_loc3nt)$mean
loc3nt_comp_KLUCB7 <- compare_to_ground_truth(mean_loc3nt, loc3nt_KLUCB7, 
                                            tau_loc3nt, 
                                            epsilon_loc3nt)$mean
save(loc3nt_comp_KLUCB30, file = paste0(current_path, "/loc3nt_comp_KLUCB30.Rda"))
save(loc3nt_comp_KLUCB370, file = paste0(current_path, "/loc3nt_comp_KLUCB370.Rda"))
save(loc3nt_comp_KLUCB7, file = paste0(current_path, "/loc3nt_comp_KLUCB7.Rda"))
