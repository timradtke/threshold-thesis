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
  curr_data <- data.frame(rep(NA, times = 1000))
  for(i in 1:length(mean_loc2)) {
    curr_data[[i]] <- as.numeric(purrr::rbernoulli(1000, p  = mean_loc2[i]))
  }
  names(curr_data) <- paste0("V", rep(1:length(mean_loc2)))
  data_list2[[j]] <- curr_data
}

########################################################################
# KL-UCB Algorithm
system.time(loc2_KLUCB60 <- para_bandit_sim_KLUCB(data = data_list2[1:100], 
                                                  rounds = 1000, 
                                                  tau = tau_loc2, 
                                                  epsilon = epsilon_loc2,
                                                  horizon = 6000,
                                                  H = H_loc2))
# user  system elapsed 
#1.019   1.369 553.760
save(loc2_KLUCB60, file = paste0(current_path, "loc2_KLUCB60.Rda"))
#load(file = paste0(current_path, "loc4_KLUCB40.Rda"))
loc2_comp_KLUCB60 <- compare_to_ground_truth(mean_loc2, loc2_KLUCB60, tau_loc2, 
                                             epsilon_loc2)$mean
save(loc2_comp_KLUCB60, file = paste0(current_path, "loc2_comp_KLUCB60.Rda"))
#load(file = paste0(current_path, "loc4_comp_KLUCB40.Rda"))
rm(loc4_APT)
gc()

########################################################################
# Standard APT Algorithm
# a seed is built in
system.time(loc2_APT <- para_bandit_sim_APT(data = data_list2, rounds = 1000, 
                                           tau = tau_loc2, epsilon = epsilon_loc2))
save(loc2_APT, file = paste0(current_path, "loc2_APT.Rda"))
#load(file = paste0(current_path, "loc2_APT.Rda"))
loc2_comp_APT <- compare_to_ground_truth(mean_loc2, loc2_APT, tau_loc2, 
                                        epsilon_loc2)$mean
save(loc2_comp_APT, file = paste0(current_path, "loc2_comp_APT.Rda"))
rm(loc2_APT)
gc()

########################################################################
# Standard AugUCB
system.time(loc2_AugUCB <- para_bandit_sim_AugUCB(data = data_list2, 
                                                  rounds = 1000, 
                                                  tau = tau_loc2))
# user   system  elapsed 
# 13.972   10.042 2284.875
save(loc2_AugUCB, file = paste0(current_path, "loc2_AugUCB.Rda"))
#load(file = paste0(current_path, "loc2_AugUCB.Rda"))
loc2_comp_AugUCB <- compare_to_ground_truth(mean_loc2, loc2_AugUCB, tau_loc2, 
                                           epsilon_loc2)$mean
save(loc2_comp_AugUCB, file = paste0(current_path, "loc2_comp_AugUCB.Rda"))
rm(loc2_AugUCB)
gc()
########################################################################
# Standard Uniform

system.time(loc2_UNIFORM <- para_bandit_sim_uniform(data = data_list2, 
                                                    rounds = 1000))
# user  system elapsed 
# 6.268   3.400 211.002 
save(loc2_UNIFORM, file = paste0(current_path, "loc2_UNIFORM.Rda"))
#load(file = paste0(current_path, "loc2_UNIFORM.Rda"))
loc2_comp_UNIFORM <- compare_to_ground_truth(mean_loc2, loc2_UNIFORM, tau_loc2, 
                                            epsilon_loc2)$mean
save(loc2_comp_UNIFORM, file = paste0(current_path, "loc2_comp_UNIFORM.Rda"))
rm(loc2_UNIFORM)
gc()
########################################################################

loc2_BUCB_horizon <- para_bandit_sim_bucb(data = data_list2, rounds = 1000, 
                                         rate = "inverse_horizon",
                                         tau = tau_loc2, epsilon = epsilon_loc2, 
                                         alpha = tau_loc2, beta = 1-tau_loc2)
save(loc2_BUCB_horizon, 
     file = paste0(current_path, "loc2_BUCB_horizon.Rda"))
loc2_comp_BUCB_horizon <- compare_to_ground_truth(mean_loc2, 
                                                 loc2_BUCB_horizon,
                                                 tau_loc2,
                                                 epsilon_loc2)$mean
save(loc2_comp_BUCB_horizon, file = paste0(current_path, "loc2_comp_BUCB_horizon.Rda"))

########################################################################

system.time(loc2_KL_tau_horizon <- para_bandit_sim_KL(data = data_list2, 
                                                     rounds = 1000, 
                                                     tau = tau_loc2, 
                                                     epsilon = epsilon_loc2, 
                                                     at_tau = TRUE, 
                                                     horizon = 1000^2))
save(loc2_KL_tau_horizon, file = paste0(current_path, "loc2_KL_tau_horizon.Rda"))
loc2_comp_KL_tau_horizon <- compare_to_ground_truth(mean_loc2, 
                                                   loc2_KL_tau_horizon, 
                                                   tau_loc2, epsilon_loc2)$mean
save(loc2_comp_KL_tau_horizon, file = paste0(current_path,
                                            "loc2_comp_KL_tau_horizon.Rda"))
rm(loc2_KL_tau_horizon)
gc()

########################################################################

# Do KL by comparing against tau directly

system.time(loc2_KL <- para_bandit_sim_KL(data = data_list2, rounds = 1000, 
                                         tau = tau_loc2, epsilon = epsilon_loc2, 
                                         at_tau = TRUE, horizon = 1000))
# user  system elapsed 
#6.932   4.007 585.704 
save(loc2_KL, file = paste0(current_path, "loc2_KL.Rda"))
#load(file = paste0(current_path, "loc2_KL.Rda"))
loc2_comp_KL <- compare_to_ground_truth(mean_loc2, loc2_KL, tau_loc2, 
                                       epsilon_loc2)$mean
save(loc2_comp_KL, file = paste0(current_path, "loc2_comp_KL.Rda"))
rm(loc2_KL)
gc()
########################################################################
# Do KL by comparing against tau and epsilon

system.time(loc2_KL_not_tau_horizon1000 <- para_bandit_sim_KL(data = data_list2, rounds = 1000, 
                                                             tau = tau_loc2, epsilon = epsilon_loc2, 
                                                             at_tau = FALSE, horizon = 1000))
save(loc2_KL_not_tau_horizon1000, file = paste0(current_path, "loc2_KL_not_tau_horizon1000.Rda"))
#load(file = paste0(current_path, "loc2_KL_not_tau_horizon1000.Rda"))
loc2_comp_KL_not_tau_horizon1000 <- compare_to_ground_truth(mean_loc2, 
                                                           loc2_KL_not_tau_horizon1000, 
                                                           tau_loc2, epsilon_loc2)$mean
save(loc2_comp_KL_not_tau_horizon1000, file = paste0(current_path,
                                                    "loc2_comp_KL_not_tau_horizon1000.Rda"))
rm(loc2_KL_not_tau_horizon1000)
gc()

########################################################################
# Probability that arm is above or below tau+-epsilon

loc2_PI <- para_bandit_sim_PI(data = data_list2, rounds = 1000,
                             tau = tau_loc2, epsilon = epsilon_loc2,
                             alpha = tau_loc2, beta = 1 - tau_loc2)
save(loc2_PI, file = paste0(current_path, "loc2_PI.Rda"))
#load(file = paste0(current_path, "loc2_PI.Rda"))
loc2_comp_PI <- compare_to_ground_truth(mean_loc2, loc2_PI, 
                                       tau_loc2, epsilon_loc2)$mean
save(loc2_comp_PI, file = paste0(current_path, "loc2_comp_PI.Rda"))
rm(loc2_PI)
gc()
########################################################################

loc2_TTS <- para_bandit_sim_TTS(data = data_list2, rounds = 1000,
                               tau = tau_loc2, epsilon = epsilon_loc2,
                               alpha = tau_loc2, beta = 1 - tau_loc2)
save(loc2_TTS, file = paste0(current_path, "loc2_TTS.Rda"))
#load(file = paste0(current_path, "loc2_TTS.Rda"))
loc2_comp_TTS <- compare_to_ground_truth(mean_loc2, loc2_TTS, 
                                        tau_loc2, epsilon_loc2)$mean
save(loc2_comp_TTS, file = paste0(current_path, "loc2_comp_TTS.Rda"))
rm(loc2_TTS)
gc()

########################################################################

loc2_BUCB <- para_bandit_sim_bucb(data = data_list2, rounds = 1000, 
                                 rate = "inverse",
                                 tau = tau_loc2, epsilon = epsilon_loc2, 
                                 alpha = tau_loc2, beta = 1-tau_loc2)
save(loc2_BUCB, file = paste0(current_path, "loc2_BUCB.Rda"))
#load(file = paste0(current_path, "loc2_BUCB.Rda"))
loc2_comp_BUCB <- compare_to_ground_truth(mean_loc2, loc2_BUCB, 
                                         tau_loc2, epsilon_loc2)$mean
save(loc2_comp_BUCB, file = paste0(current_path, "loc2_comp_BUCB.Rda"))
rm(loc2_BUCB)
gc()
########################################################################

loc2_BUCB_squared <- para_bandit_sim_bucb(data = data_list2, rounds = 1000, 
                                         rate = "inverse_squared",
                                         tau = tau_loc2, epsilon = epsilon_loc2, 
                                         alpha = tau_loc2, beta = 1-tau_loc2)
save(loc2_BUCB_squared, file = paste0(current_path, "loc2_BUCB_squared.Rda"))
#load(file = paste0(current_path, "loc2_BUCB_squared.Rda"))
loc2_comp_BUCB_squared <- compare_to_ground_truth(mean_loc2, loc2_BUCB_squared, 
                                                 tau_loc2, epsilon_loc2)$mean
save(loc2_comp_BUCB_squared, file = paste0(current_path, 
                                          "loc2_comp_BUCB_squared.Rda"))
rm(loc2_BUCB_squared)
gc()
########################################################################

loc2_BUCB_power5 <- para_bandit_sim_bucb(data = data_list2, rounds = 1000, 
                                        rate = "inverse_power5",
                                        tau = tau_loc2, epsilon = epsilon_loc2, 
                                        alpha = tau_loc2, beta = 1-tau_loc2)
save(loc2_BUCB_power5, file = paste0(current_path, "loc2_BUCB_power5.Rda"))
load(file = paste0(current_path, "loc2_BUCB_power5.Rda"))
loc2_comp_BUCB_power5 <- compare_to_ground_truth(mean_loc2, loc2_BUCB_power5, 
                                                tau_loc2, epsilon_loc2)$mean
save(loc2_comp_BUCB_power5, file = paste0(current_path, 
                                         "loc2_comp_BUCB_power5.Rda"))
rm(loc2_BUCB_power5)
gc()
########################################################################

load(paste0(current_path, "loc2_comp_APT.Rda"))
load(paste0(current_path, "loc2_comp_AugUCB.Rda"))
load(paste0(current_path, "loc2_comp_UNIFORM.Rda"))
load(paste0(current_path, "loc2_comp_PI.Rda"))
load(paste0(current_path, "loc2_comp_TTS.Rda"))
load(paste0(current_path, "loc2_comp_BUCB.Rda"))
load(paste0(current_path, "loc2_comp_BUCB_horizon.Rda"))
load(paste0(current_path, "loc2_comp_BUCB_squared.Rda"))
#load(paste0(current_path, "loc2_comp_BUCB_power5.Rda"))
load(paste0(current_path, "loc2_comp_KL.Rda"))
load(paste0(current_path, "loc2_comp_KL_tau_horizon.Rda"))
load(paste0(current_path, "loc2_comp_KL_horizon1.Rda"))
load(paste0(current_path, "loc2_comp_KL_not_tau_horizon1000.Rda"))

plot(c(0,1000), c(0, -7), type = "n")
lines(log(loc2_comp_APT), col = "red")
lines(log(loc2_comp_AugUCB), col = "grey")
lines(log(loc2_comp_UNIFORM), col = "black")
lines(log(loc2_comp_PI), col = "darkgreen")
lines(log(loc2_comp_TTS), col = "lightgreen")
lines(log(loc2_comp_BUCB), col = "orange")
lines(log(loc2_comp_BUCB_horizon), col = "darkorange")
lines(log(loc2_comp_BUCB_squared), col = "blue")
lines(log(loc2_comp_KLUCB10), col = "blue")
lines(log(loc2_comp_KLUCB185), col = "darkblue")
lines(log(loc2_comp_KLUCB60), col = "lightblue")
#lines(log(loc2_comp_BUCB_power5), col = "darkred")
lines(log(loc2_comp_KL), col = "blue", lty = 2)
lines(log(loc2_comp_KL_tau_horizon), col = "black", lty = 2)
#lines(log(loc2_comp_KL_horizon1), col = "black", lty = 2)
lines(log(loc2_comp_KL_not_tau_horizon1000), col = "black", lty = 2)

########################################################################

library(plyr)
loc2_UNIFORM_as <- colMeans(ldply(loc2_UNIFORM, function(x) table(x$arm_sequence)))
loc2_APT_as <- colMeans(ldply(loc2_APT, function(x) table(x$arm_sequence)))
loc2_PI_as <- colMeans(ldply(loc2_PI, function(x) table(x$arm_sequence)))
loc2_TTS_as <- colMeans(ldply(loc2_TTS, function(x) table(x$arm_sequence)))
loc2_BUCB_as <- colMeans(ldply(loc2_BUCB, function(x) table(x$arm_sequence)))
loc2_KL_as <- colMeans(ldply(loc2_KL, function(x) table(x$arm_sequence)))
loc2_KL_not_tau_as <- colMeans(ldply(loc2_KL_not_tau, function(x) table(x$arm_sequence)))
loc2_KL_horizon_as <- colMeans(ldply(loc2_KL_horizon, function(x) table(x$arm_sequence)))

round(data.frame(loc2_UNIFORM_as, loc2_APT_as, loc2_PI_as, loc2_TTS_as, loc2_BUCB_as, 
                 loc2_KL_as, loc2_KL_not_tau_as, loc2_KL_horizon_as))
