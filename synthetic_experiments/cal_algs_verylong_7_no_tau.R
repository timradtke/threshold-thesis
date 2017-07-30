########################################################################
# load functions

current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/AugUCB/"
source(paste0(current_path, "run_from_data.R"))
source(paste0(current_path, "run_parallelized_from_data.R"))
current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/synthetic_experiments/"

########################################################################
# Create the data

# do 2000 rounds
mean_loc7nt <- c(10^-4, 10^-3.5, 10^-3.25, 10^-3, 
                 10^-2.75, 10^-2.5, 10^-2.25, 10^-1.85, 10^-1.75,
                 10^-1.5, 10^-1)
tau_loc7nt <- 10^-2
epsilon_loc7nt <- 0
H_loc7nt <- get_complexity(mean_loc7nt, tau_loc7nt, epsilon_loc7nt)

plot(mean_loc7nt, rep(1,11), main = paste0("Complexity of ", round(H_loc7nt,2)))
abline(v=tau_loc7nt)
abline(v=tau_loc7nt+epsilon_loc7nt, lty=2)
abline(v=tau_loc7nt-epsilon_loc7nt, lty=2)

plot(log10(mean_loc7nt), rep(1,11), main = paste0("Complexity of ", round(H_loc7nt,2)))
abline(v=log10(tau_loc7nt))

data_list7_10000 <- list()
set.seed(386303469)
for(j in 1:5000) {
  curr_data <- data.frame(rep(NA, times = 10000))
  for(i in 1:length(mean_loc7nt)) {
    curr_data[[i]] <- as.numeric(purrr::rbernoulli(10000, p  = mean_loc7nt[i]))
  }
  names(curr_data) <- paste0("V", rep(1:length(mean_loc7nt)))
  data_list7_10000[[j]] <- curr_data
}
gc()
########################################################################

system.time(loc7nt_KLUCB90_long <- para_bandit_sim_KLUCB(data = data_list7_10000, 
                                                       rounds = 10000, 
                                                       tau = tau_loc7nt, 
                                                       epsilon = epsilon_loc7nt,
                                                       horizon = 90000,
                                                       H = H_loc7nt))
save(loc7nt_KLUCB90_long, file = paste0(current_path, "loc7nt_KLUCB90_long.Rda"))
loc7nt_comp_KLUCB90_long <- compare_to_ground_truth(mean_loc7nt, loc7nt_KLUCB90_long, 
                                                  tau_loc7nt, 
                                                  epsilon_loc7nt)$mean
save(loc7nt_comp_KLUCB90_long, file = paste0(current_path, 
                                           "loc7nt_comp_KLUCB90_long.Rda"))

########################################################################
# Standard Likelihood Ratio
system.time(loc7nt_LR_10000 <- para_bandit_sim_LR(data = data_list7_10000, 
                                                 rounds = 10000, 
                                                 tau = tau_loc7nt, 
                                                 epsilon = epsilon_loc7nt))

save(loc7nt_LR_10000, file = paste0(current_path, "loc7nt_LR_10000.Rda"))
loc7nt_comp_LR_10000 <- compare_to_ground_truth(mean_loc7nt, loc7nt_LR_10000, 
                                               tau_loc7nt, epsilon_loc7nt)$mean
save(loc7nt_comp_LR_10000, file = paste0(current_path, 
                                        "loc7nt_comp_LR_10000.Rda"))
rm(loc7nt_LR_10000, loc7nt_comp_LR_10000)
gc()

########################################################################
# Standard Uniform
system.time(loc7nt_UNIFORM_10000 <- para_bandit_sim_uniform(data = data_list7_10000, 
                                                         rounds = 10000))
save(loc7nt_UNIFORM_10000, file = paste0(current_path, "loc7nt_UNIFORM_10000.Rda"))
loc7nt_comp_UNIFORM_10000 <- compare_to_ground_truth(mean_loc7nt, loc7nt_UNIFORM_10000, 
                                                    tau_loc7nt, 
                                                    epsilon_loc7nt)$mean
save(loc7nt_comp_UNIFORM_10000, file = paste0(current_path,
                                           "loc7nt_comp_UNIFORM_10000.Rda"))
rm(loc7nt_UNIFORM_10000, loc7nt_comp_UNIFORM_10000)
gc()
########################################################################
# Standard APT Algorithm
# a seed is built in
system.time(loc7nt_APT_10000 <- para_bandit_sim_APT(data = data_list7_10000, 
                                                 rounds = 10000, 
                                                 tau = tau_loc7nt, 
                                                 epsilon = epsilon_loc7nt))
save(loc7nt_APT_10000, file = paste0(current_path, "loc7nt_APT_10000.Rda"))
loc7nt_comp_APT_10000 <- compare_to_ground_truth(mean_loc7nt, loc7nt_APT_10000, 
                                              tau_loc7nt, 
                                              epsilon_loc7nt)$mean
save(loc7nt_comp_APT_10000, file = paste0(current_path, 
                                       "loc7nt_comp_APT_10000.Rda"))

########################################################################

loc7nt_BUCB_horizon_10000 <- para_bandit_sim_bucb(data = data_list7_10000, 
                                               rounds = 10000, 
                                               rate = "inverse_horizon",
                                               tau = tau_loc7nt, 
                                               epsilon = epsilon_loc7nt, 
                                               alpha = tau_loc7nt, 
                                               beta = 1-tau_loc7nt)
save(loc7nt_BUCB_horizon_10000, 
     file = paste0(current_path, "loc7nt_BUCB_horizon_10000.Rda"))
loc7nt_comp_BUCB_horizon_10000 <- compare_to_ground_truth(mean_loc7nt, 
                                                       loc7nt_BUCB_horizon_10000,
                                                       tau_loc7nt,
                                                       epsilon_loc7nt)$mean
save(loc7nt_comp_BUCB_horizon_10000,
     file = paste0(current_path, "loc7nt_comp_BUCB_horizon_10000.Rda"))
rm(loc7nt_BUCB_horizon_10000)
gc()

########################################################################
# Standard AugUCB
system.time(loc7nt_AugUCB <- para_bandit_sim_AugUCB(data = data_list7, 
                                                  rounds = 2000, 
                                                  tau = tau_loc7nt))
save(loc7nt_AugUCB, file = paste0(current_path, "loc7nt_AugUCB.Rda"))
load(file = paste0(current_path, "loc7nt_AugUCB.Rda"))
loc7nt_comp_AugUCB <- compare_to_ground_truth(mean_loc7nt, loc7nt_AugUCB, tau_loc7nt, 
                                            epsilon_loc7nt)$mean
save(loc7nt_comp_AugUCB, file = paste0(current_path, "loc7nt_comp_AugUCB.Rda"))

# now with 10000

system.time(loc7nt_AugUCB_10000 <- para_bandit_sim_AugUCB(data = data_list7_10000, 
                                                       rounds = 10000, 
                                                       tau = tau_loc7nt))
#user    system   elapsed 
#60.824    67.047 24179.679 
save(loc7nt_AugUCB_10000, file = paste0(current_path, "loc7nt_AugUCB_10000.Rda"))
#load(file = paste0(current_path, "loc7nt_AugUCB_10000.Rda"))
loc7nt_comp_AugUCB_10000 <- compare_to_ground_truth(mean_loc7nt, loc7nt_AugUCB_10000, 
                                                 tau_loc7nt, epsilon_loc7nt)$mean
save(loc7nt_comp_AugUCB_10000, file = paste0(current_path, 
                                          "loc7nt_comp_AugUCB_10000.Rda"))

########################################################################

load(file = paste0(current_path, "loc7nt_comp_UNIFORM_10000.Rda"))
load(file = paste0(current_path, "loc7nt_comp_APT_10000.Rda"))
load(file = paste0(current_path, "loc7nt_comp_LR_10000.Rda"))
load(file = paste0(current_path, "loc7nt_comp_AugUCB_10000.Rda"))

plot(c(0,10000), c(0, -5), type = "n")
lines(log(loc7nt_comp_BUCB_horizon_10000), col = "red")
lines(log(loc7nt_comp_UNIFORM_10000), col = "darkblue")
lines(log(loc7nt_comp_APT_10000), col = "green")
lines(log(loc7nt_comp_LR_10000), col = "darkred")
