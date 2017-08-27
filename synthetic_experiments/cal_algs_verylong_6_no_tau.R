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

mean_loc6nt <- c(0.001, 0.005, 0.01, 0.015,
                 0.04, 0.06,
                 0.085, 0.09, 0.095, 0.099)
tau_loc6nt <- 0.05
epsilon_loc6nt <- 0
H_loc6nt <- get_complexity(mean_loc6nt, tau_loc6nt, epsilon_loc6nt)

plot(mean_loc6nt, rep(1,10), main = paste0("Complexity of ", round(H_loc6nt,2)))
abline(v=tau_loc6nt)
abline(v=tau_loc6nt+epsilon_loc6nt, lty=2)
abline(v=tau_loc6nt-epsilon_loc6nt, lty=2)

data_list6 <- list()
set.seed(8247502)
for(j in 1:5000) {
  curr_data <- data.frame(rep(NA, times = 7000))
  for(i in 1:length(mean_loc6nt)) {
    curr_data[[i]] <- as.numeric(purrr::rbernoulli(7000, p  = mean_loc6nt[i]))
  }
  names(curr_data) <- paste0("V", rep(1:length(mean_loc6nt)))
  data_list6[[j]] <- curr_data
}

########################################################################
# Empirical Variance Guided Algorithm (Zhong et al., 2017)

loc6nt_EVT <- para_bandit_sim_EVT(data = data_list6, 
                                  rounds = 7000, 
                                  tau = tau_loc6nt, 
                                  epsilon = epsilon_loc6nt)
save(loc6nt_EVT, file = paste0(current_path, "loc6nt_EVT.Rda"))
loc6nt_comp_EVT <- compare_to_ground_truth(mean_loc6nt, loc6nt_EVT, 
                                           tau_loc6nt, epsilon_loc6nt)$mean
save(loc6nt_comp_EVT, file = paste0(current_path, "loc6nt_comp_EVT.Rda"))
rm(loc6nt_EVT)
gc()

########################################################################
# Plot the results

load(paste0(current_path, "/loc6nt_comp_BUCB_horizon_7000.Rda"))
load(paste0(current_path, "/loc6nt_comp_BUCB_horizon_c5_7000.Rda"))
load(paste0(current_path, "/loc6nt_comp_BUCB_horizon_c15_7000.Rda"))
load(paste0(current_path, "/loc6nt_comp_AugUCB_7000.Rda"))
load(paste0(current_path, "/loc6nt_comp_KLUCB370_5000.Rda"))
load(paste0(current_path, "/loc6nt_comp_UNIFORM_7000.Rda"))
load(paste0(current_path, "/loc6nt_comp_APT_7000.Rda"))
load(paste0(current_path, "/loc6nt_comp_LR_7000.Rda"))

plot(c(0,7000), c(0, -5), type = "n")
lines(log(loc6nt_comp_BUCB_horizon_7000), col = "blue")
lines(log(loc6nt_comp_BUCB_horizon_c5_7000), col = "lightblue")
lines(log(loc6nt_comp_BUCB_horizon_c15_7000), col = "darkblue")
#lines(log(loc6nt_comp_BUCB_horizon), col = "lightblue")
#lines(log(loc6nt_comp_APT), col = "pink")
lines(log(loc6nt_comp_LR_7000), col = "green")
lines(log(loc6nt_comp_APT_7000), col = "red")
lines(log(loc6nt_comp_UNIFORM_7000), col = "black")
abline(h=log(0.1))
#lines(log(loc6nt_comp_AugUCB), col = "grey")
lines(log(loc6nt_comp_AugUCB_7000), col = "grey")
#lines(log(loc6nt_comp_KLUCB7), col = "green")
#lines(log(loc6nt_comp_KLUCB30), col = "green")
#lines(log(loc6nt_comp_KLUCB370), col = "darkgreen")
#lines(log(loc6nt_comp_KLUCB370_5000), col = "darkgreen")

########################################################################
# Empirical Variance Guided Algorithm (Zhong et al., 2017)

loc6nt_EVT <- para_bandit_sim_EVT(data = data_list6, 
                                rounds = 7000, 
                                tau = tau_loc6nt, 
                                epsilon = epsilon_loc6nt)
save(loc6nt_EVT, file = paste0(current_path, "loc6nt_EVT.Rda"))
loc6nt_comp_EVT <- compare_to_ground_truth(mean_loc6nt, loc6nt_EVT, 
                                           tau_loc6nt, epsilon_loc6nt)$mean
save(loc6nt_comp_EVT, file = paste0(current_path, "loc6nt_comp_EVT.Rda"))
rm(loc6nt_EVT)
gc()

########################################################################

loc6nt_BUCB_horizon_7000 <- para_bandit_sim_bucb(data = data_list6, 
                                                 rounds = 7000, 
                                                 rate = "inverse_horizon",
                                                 tau = tau_loc6nt, 
                                                 epsilon = epsilon_loc6nt, 
                                                 alpha = tau_loc6nt, 
                                                 beta = 1-tau_loc6nt)
save(loc6nt_BUCB_horizon_7000, 
     file = paste0(current_path, "loc6nt_BUCB_horizon_7000.Rda"))
loc6nt_comp_BUCB_horizon_7000 <- compare_to_ground_truth(mean_loc6nt, 
                                                         loc6nt_BUCB_horizon_7000,
                                                         tau_loc6nt,
                                                         epsilon_loc6nt)$mean
save(loc6nt_comp_BUCB_horizon_7000, file = paste0(current_path, 
                                                  "loc6nt_comp_BUCB_horizon_7000.Rda"))

########################################################################
# Try a different exploration parameter for BUCB

loc6nt_BUCB_horizon_c5_7000 <- para_bandit_sim_bucb(data = data_list6, 
                                                 rounds = 7000, 
                                                 rate = "inverse_horizon_c",
                                                 tau = tau_loc6nt, 
                                                 epsilon = epsilon_loc6nt, 
                                                 alpha = tau_loc6nt, 
                                                 beta = 1-tau_loc6nt,
                                                 const = 5)
save(loc6nt_BUCB_horizon_c5_7000, 
     file = paste0(current_path, "loc6nt_BUCB_horizon_c5_7000.Rda"))
loc6nt_comp_BUCB_horizon_c5_7000 <- compare_to_ground_truth(mean_loc6nt, 
                                                         loc6nt_BUCB_horizon_c5_7000,
                                                         tau_loc6nt,
                                                         epsilon_loc6nt)$mean
save(loc6nt_comp_BUCB_horizon_c5_7000, file = paste0(current_path, 
                                                  "loc6nt_comp_BUCB_horizon_c5_7000.Rda"))
rm(loc6nt_BUCB_horizon_c5_7000)

########################################################################
# Try a different exploration parameter for BUCB

loc6nt_BUCB_horizon_c15_7000 <- para_bandit_sim_bucb(data = data_list6, 
                                                    rounds = 7000, 
                                                    rate = "inverse_horizon_c",
                                                    tau = tau_loc6nt, 
                                                    epsilon = epsilon_loc6nt, 
                                                    alpha = tau_loc6nt, 
                                                    beta = 1-tau_loc6nt,
                                                    const = 1/5)
save(loc6nt_BUCB_horizon_c15_7000, 
     file = paste0(current_path, "loc6nt_BUCB_horizon_c15_7000.Rda"))
loc6nt_comp_BUCB_horizon_c15_7000 <- compare_to_ground_truth(mean_loc6nt, 
                                                            loc6nt_BUCB_horizon_c15_7000,
                                                            tau_loc6nt,
                                                            epsilon_loc6nt)$mean
save(loc6nt_comp_BUCB_horizon_c15_7000, file = paste0(current_path, 
                                                     "loc6nt_comp_BUCB_horizon_c15_7000.Rda"))
rm(loc6nt_BUCB_horizon_c15_7000)

########################################################################
# Standard AugUCB
loc6nt_AugUCB_7000 <- para_bandit_sim_AugUCB(data = data_list6, 
                                             rounds = 7000, 
                                             tau = tau_loc6nt)
save(loc6nt_AugUCB_7000, file = paste0(current_path, "loc6nt_AugUCB_7000.Rda"))
loc6nt_comp_AugUCB_7000 <- compare_to_ground_truth(mean_loc6nt, loc6nt_AugUCB_7000, 
                                                   tau_loc6nt, epsilon_loc6nt)$mean
save(loc6nt_comp_AugUCB_7000, file = paste0(current_path, 
                                            "loc6nt_comp_AugUCB_7000.Rda"))

########################################################################
# Standard APT
system.time(loc6nt_APT_7000 <- para_bandit_sim_APT(data = data_list6, 
                                                   rounds = 7000, 
                                                   tau = tau_loc6nt, 
                                                   epsilon = epsilon_loc6nt))

save(loc6nt_APT_7000, file = paste0(current_path, "loc6nt_APT_7000.Rda"))
loc6nt_comp_APT_7000 <- compare_to_ground_truth(mean_loc6nt, loc6nt_APT_7000, 
                                                tau_loc6nt, epsilon_loc6nt)$mean
save(loc6nt_comp_APT_7000, file = paste0(current_path, 
                                         "loc6nt_comp_APT_7000.Rda"))

########################################################################
# Standard Likelihood Ratio
system.time(loc6nt_LR_7000 <- para_bandit_sim_LR(data = data_list6, 
                                                 rounds = 7000, 
                                                 tau = tau_loc6nt, 
                                                 epsilon = epsilon_loc6nt))

save(loc6nt_LR_7000, file = paste0(current_path, "loc6nt_LR_7000.Rda"))
loc6nt_comp_LR_7000 <- compare_to_ground_truth(mean_loc6nt, loc6nt_LR_7000, 
                                                tau_loc6nt, epsilon_loc6nt)$mean
save(loc6nt_comp_LR_7000, file = paste0(current_path, 
                                        "loc6nt_comp_LR_7000.Rda"))

########################################################################

system.time(loc6nt_UNIFORM_7000 <- para_bandit_sim_uniform(data = data_list6, 
                                                         rounds = 7000))
save(loc6nt_UNIFORM_7000, file = paste0(current_path, "loc6nt_UNIFORM_7000.Rda"))
loc6nt_comp_UNIFORM_7000 <- compare_to_ground_truth(mean_loc6nt, 
                                                    loc6nt_UNIFORM_7000, 
                                                    tau_loc6nt, 
                                                    epsilon_loc6nt)$mean
save(loc6nt_comp_UNIFORM_7000, file = paste0(current_path, 
                                             "loc6nt_comp_UNIFORM_7000.Rda"))

########################################################################

system.time(loc6nt_KLUCB370_5000 <- para_bandit_sim_KLUCB(data = data_list6, 
                                                          rounds = 5000, 
                                                          tau = tau_loc6nt, 
                                                          epsilon = epsilon_loc6nt,
                                                          horizon = 370000,
                                                          H = H_loc6nt))

save(loc6nt_KLUCB370_5000, file = paste0(current_path, 
                                         "loc6nt_KLUCB370_5000.Rda"))
loc6nt_comp_KLUCB370_5000 <- compare_to_ground_truth(mean_loc6nt, 
                                                     loc6nt_KLUCB370_5000, 
                                                     tau_loc6nt, epsilon_loc6nt)$mean
save(loc6nt_comp_KLUCB370_5000, file = paste0(current_path, 
                                              "loc6nt_comp_KLUCB370_5000.Rda"))

########################################################################

error_bound <- function(n, H, K) {
  (K*n+H)*exp(-n/H)
}

plot(1:500000, log(error_bound(1:500000, H_loc6nt, 10)), type = "l",
     ylim = c(-6,12))
abline(h=log(0.2))

data_list6_small <- data_list6[1:500]
rm(data_list6)
gc()

system.time(loc6nt_KLUCB30 <- para_bandit_sim_KLUCB(data = data_list6_small, 
                                                    rounds = 2000, 
                                                    tau = tau_loc6nt, 
                                                    epsilon = epsilon_loc6nt,
                                                    horizon = round(H_loc6nt),
                                                    H = H_loc6nt))

system.time(loc6nt_KLUCB370 <- para_bandit_sim_KLUCB(data = data_list6_small, 
                                                     rounds = 2000, 
                                                     tau = tau_loc6nt, 
                                                     epsilon = epsilon_loc6nt,
                                                     horizon = 370000,
                                                     H = H_loc6nt))

system.time(loc6nt_KLUCB7 <- para_bandit_sim_KLUCB(data = data_list6_small, 
                                                   rounds = 2000, 
                                                   tau = tau_loc6nt, 
                                                   epsilon = epsilon_loc6nt,
                                                   horizon = 7000,
                                                   H = H_loc6nt))

save(loc6nt_KLUCB30, file = paste0(current_path, "/loc6nt_KLUCB30.Rda"))
save(loc6nt_KLUCB370, file = paste0(current_path, "/loc6nt_KLUCB370.Rda"))
save(loc6nt_KLUCB7, file = paste0(current_path, "/loc6nt_KLUCB7.Rda"))

loc6nt_comp_KLUCB30 <- compare_to_ground_truth(mean_loc6nt, loc6nt_KLUCB30, 
                                               tau_loc6nt, 
                                               epsilon_loc6nt)$mean
loc6nt_comp_KLUCB370 <- compare_to_ground_truth(mean_loc6nt, loc6nt_KLUCB370, 
                                                tau_loc6nt, 
                                                epsilon_loc6nt)$mean
loc6nt_comp_KLUCB7 <- compare_to_ground_truth(mean_loc6nt, loc6nt_KLUCB7, 
                                              tau_loc6nt, 
                                              epsilon_loc6nt)$mean
save(loc6nt_comp_KLUCB30, file = paste0(current_path, "/loc6nt_comp_KLUCB30.Rda"))
save(loc6nt_comp_KLUCB370, file = paste0(current_path, "/loc6nt_comp_KLUCB370.Rda"))
save(loc6nt_comp_KLUCB7, file = paste0(current_path, "/loc6nt_comp_KLUCB7.Rda"))
