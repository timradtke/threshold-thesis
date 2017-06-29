# run a simulation based on a simple example that can hopefully serve to
# calibrate the algorithms
########################################################################
# load functions

current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/AugUCB/"
source(paste0(current_path, "run_from_data.R"))
source(paste0(current_path, "run_parallelized_from_data.R"))

########################################################################
# Create the data

mean_loc2 <- c(0.2+(0:3)*0.05,
              0.45, 0.55,
              0.65+(0:3)*0.05)
tau_loc2 <- 0.5
epsilon_loc2 <- 0.1

data_list2 <- list()
for(j in 1:5000) {
  set.seed(1024+j)
  curr_data <- data.frame(rep(NA, times = 1000))
  for(i in 1:length(mean_loc2)) {
    curr_data[[i]] <- as.numeric(purrr::rbernoulli(1000, p  = mean_loc2[i]))
  }
  names(curr_data) <- paste0("V", rep(1:length(mean_loc2)))
  data_list2[[j]] <- curr_data
}

########################################################################
# Standard APT Algorithm
# a seed is built in
system.time(loc2_APT <- para_bandit_sim_APT(data = data_list2, rounds = 1000, 
                                           tau = tau_loc2, epsilon = epsilon_loc2))
# user  system elapsed 
#6.508   3.544 300.383
save(loc2_APT, file = paste0(current_path, "loc2_APT.Rda"))
loc2_APT[[1]]
loc2_APT[[5000]]

loc2_comp_APT <- compare_to_ground_truth(mean_loc2, loc2_APT, tau_loc2, 
                                        epsilon_loc2)$mean
save(loc2_comp_APT, file = paste0(current_path, "loc2_comp_APT.Rda"))
########################################################################
# Standard Uniform

system.time(loc2_UNIFORM <- para_bandit_sim_uniform(data = data_list2, 
                                                    rounds = 1000))
# user  system elapsed 
# 6.268   3.400 211.002 
save(loc2_UNIFORM, file = paste0(current_path, "loc2_UNIFORM.Rda"))
loc2_comp_UNIFORM <- compare_to_ground_truth(mean_loc2, loc2_UNIFORM, tau_loc2, 
                                            epsilon_loc2)$mean
save(loc2_comp_UNIFORM, file = paste0(current_path, "loc2_comp_UNIFORM.Rda"))

########################################################################
# Do KL by comparing against tau directly

system.time(loc2_KL <- para_bandit_sim_KL(data = data_list2, rounds = 1000, 
                                         tau = tau_loc2, epsilon = epsilon_loc2, 
                                         at_tau = TRUE))
# user  system elapsed 
#6.824   4.335 391.240
save(loc2_KL, file = paste0(current_path, "loc2_KL.Rda"))
loc2_comp_KL <- compare_to_ground_truth(mean_loc2, loc2_KL, tau_loc2, 
                                       epsilon_loc2)$mean
save(loc2_comp_KL, file = paste0(current_path, "loc2_comp_KL.Rda"))

########################################################################
# Do KL by comparing against tau and epsilon

system.time(loc2_KL_not_tau <- para_bandit_sim_KL(data = data_list2, rounds = 1000, 
                                                 tau = tau_loc2, epsilon = epsilon_loc2, 
                                                 at_tau = FALSE))
save(loc2_KL_not_tau, file = paste0(current_path, "loc2_KL_not_tau.Rda"))
loc2_comp_KL_not_tau <- compare_to_ground_truth(mean_loc2, loc2_KL_not_tau, 
                                               tau_loc2, epsilon_loc2)$mean
save(loc2_comp_KL_not_tau, file = paste0(current_path,
                                        "loc2_comp_KL_not_tau.Rda"))

########################################################################
# Do KL by comparing against tau and epsilon with 0/1 adjustment

system.time(loc2_KL_horizon <- para_bandit_sim_KL(data = data_list2, rounds = 1000, 
                                                 tau = tau_loc2, epsilon = epsilon_loc2, 
                                                 at_tau = FALSE, horizon = 10000))
save(loc2_KL_horizon, file = paste0(current_path, "loc2_KL_horizon.Rda"))
loc2_comp_KL_horizon <- compare_to_ground_truth(mean_loc2, loc2_KL_horizon, 
                                               tau_loc2, epsilon_loc2)$mean
save(loc2_comp_KL_horizon, file = paste0(current_path,
                                        "loc2_comp_KL_horizon.Rda"))


########################################################################
# Do KL by comparing against tau and epsilon with 0/1 adjustment

system.time(loc2_KL_horizon2 <- para_bandit_sim_KL(data = data_list2, rounds = 1000, 
                                                  tau = tau_loc2, epsilon = epsilon_loc2, 
                                                  at_tau = FALSE, horizon = 1000^2))
save(loc2_KL_horizon2, file = paste0(current_path, "loc2_KL_horizon2.Rda"))
loc2_comp_KL_horizon2 <- compare_to_ground_truth(mean_loc2, loc2_KL_horizon2, 
                                                tau_loc2, epsilon_loc2)$mean
save(loc2_comp_KL_horizon2, file = paste0(current_path,
                                         "loc2_comp_KL_horizon2.Rda"))

########################################################################
# Do KL by comparing against tau and epsilon with 0/1 adjustment

system.time(loc2_KL_horizon1 <- para_bandit_sim_KL(data = data_list2, rounds = 1000, 
                                                  tau = tau_loc2, epsilon = epsilon_loc2, 
                                                  at_tau = FALSE, horizon = 1000))
save(loc2_KL_horizon1, file = paste0(current_path, "loc2_KL_horizon1.Rda"))
loc2_comp_KL_horizon1 <- compare_to_ground_truth(mean_loc2, loc2_KL_horizon1, 
                                                tau_loc2, epsilon_loc2)$mean
save(loc2_comp_KL_horizon1, file = paste0(current_path,
                                         "loc2_comp_KL_horizon1.Rda"))

########################################################################
# Probability that arm is above or below tau+-epsilon

loc2_PI <- para_bandit_sim_PI(data = data_list2, rounds = 1000,
                             tau = tau_loc2, epsilon = epsilon_loc2,
                             alpha = tau_loc2, beta = 1 - tau_loc2)
save(loc2_PI, file = paste0(current_path, "loc2_PI.Rda"))
loc2_comp_PI <- compare_to_ground_truth(mean_loc2, loc2_PI, 
                                       tau_loc2, epsilon_loc2)$mean
save(loc2_comp_PI, file = paste0(current_path, "loc2_comp_PI.Rda"))

########################################################################

loc2_TTS <- para_bandit_sim_TTS(data = data_list2, rounds = 1000,
                               tau = tau_loc2, epsilon = epsilon_loc2,
                               alpha = tau_loc2, beta = 1 - tau_loc2)
save(loc2_TTS, file = paste0(current_path, "loc2_TTS.Rda"))
loc2_comp_TTS <- compare_to_ground_truth(mean_loc2, loc2_TTS, 
                                        tau_loc2, epsilon_loc2)$mean
save(loc2_comp_TTS, file = paste0(current_path, "loc2_comp_TTS.Rda"))

########################################################################

#system.time(loc2_BETA <- para_bandit_sim_BETA(data = data_list2, rounds = 1000,
#                                             tau = tau_loc2, epsilon = epsilon_loc2,
#                                             alpha = tau_loc2, beta = 1 - tau_loc2))
#save(loc2_BETA, file = paste0(current_path, "loc2_BETA.Rda"))
#loc2_comp_BETA <- compare_to_ground_truth(mean_loc2, loc2_BETA, 
#                                         tau_loc2, epsilon_loc2)$mean
#save(loc2_comp_BETA, file = paste0(current_path, "loc2_comp_BETA.Rda"))

########################################################################

loc2_BUCB <- para_bandit_sim_bucb(data = data_list2, rounds = 1000, 
                                 rate = "inverse",
                                 tau = tau_loc2, epsilon = epsilon_loc2, 
                                 alpha = tau_loc2, beta = 1-tau_loc2)
save(loc2_BUCB, file = paste0(current_path, "loc2_BUCB.Rda"))
loc2_comp_BUCB <- compare_to_ground_truth(mean_loc2, loc2_BUCB, 
                                         tau_loc2, epsilon_loc2)$mean
save(loc2_comp_BUCB, file = paste0(current_path, "loc2_comp_BUCB.Rda"))

########################################################################

plot(c(0,600), c(0, -10), type = "n")
lines(log(loc2_comp_APT), col = "red")
lines(log(loc2_comp_UNIFORM), col = "black")
lines(log(loc2_comp_PI), col = "darkgreen")
lines(log(loc2_comp_TTS), col = "lightgreen")
lines(log(loc2_comp_BUCB), col = "orange")
#lines(log(loc2_comp_KL), col = "lightblue")
lines(log(loc2_comp_KL_not_tau), col = "blue")
lines(log(loc2_comp_KL_horizon), col = "darkblue")
lines(log(loc2_comp_KL_horizon2), col = "darkred")
lines(log(loc2_comp_KL_horizon1), col = "red")

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
