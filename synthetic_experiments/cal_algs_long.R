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

mean_loc <- c(0.1, 0.1, 0.1, 
              0.35, 0.45, 0.55, 0.65,
              0.9, 0.9, 0.9)
tau_loc <- 0.5
epsilon_loc <- 0.1

data_list <- list()
set.seed(1024)
for(j in 1:5000) {
  curr_data <- data.frame(rep(NA, times = 1000))
  for(i in 1:length(mean_loc)) {
    curr_data[[i]] <- as.numeric(purrr::rbernoulli(1000, p  = mean_loc[i]))
  }
  names(curr_data) <- paste0("V", rep(1:length(mean_loc)))
  data_list[[j]] <- curr_data
}

########################################################################
# Standard APT Algorithm
# a seed is built in
system.time(loc_APT <- para_bandit_sim_APT(data = data_list, rounds = 1000, 
                                           tau = tau_loc, epsilon = epsilon_loc))
# user  system elapsed 
#6.508   3.544 300.383
save(loc_APT, file = paste0(current_path, "loc_APT.Rda"))
load(file = paste0(current_path, "loc_APT.Rda"))
loc_comp_APT <- compare_to_ground_truth(mean_loc, loc_APT, tau_loc, 
                                        epsilon_loc)$mean
save(loc_comp_APT, file = paste0(current_path, "loc_comp_APT.Rda"))
rm(loc_APT)
gc()

########################################################################
# Standard AugUCB
system.time(loc_AugUCB <- para_bandit_sim_AugUCB(data = data_list, rounds = 1000, 
                                                 tau = tau_loc))
# user   system  elapsed 
# 13.972   10.042 2284.875
save(loc_AugUCB, file = paste0(current_path, "loc_AugUCB.Rda"))
load(file = paste0(current_path, "loc_AugUCB.Rda"))
loc_comp_AugUCB <- compare_to_ground_truth(mean_loc, loc_AugUCB, tau_loc, 
                                           epsilon_loc)$mean
save(loc_comp_AugUCB, file = paste0(current_path, "loc_comp_AugUCB.Rda"))
rm(loc_AugUCB)
gc()
########################################################################
# Standard Uniform

system.time(loc_UNIFORM <- para_bandit_sim_uniform(data = data_list, 
                                                   rounds = 1000))
# user  system elapsed 
# 6.268   3.400 211.002 
save(loc_UNIFORM, file = paste0(current_path, "loc_UNIFORM.Rda"))
load(file = paste0(current_path, "loc_UNIFORM.Rda"))
loc_comp_UNIFORM <- compare_to_ground_truth(mean_loc, loc_UNIFORM, tau_loc, 
                                            epsilon_loc)$mean
save(loc_comp_UNIFORM, file = paste0(current_path, "loc_comp_UNIFORM.Rda"))
rm(loc_UNIFORM)
gc()
########################################################################

loc_BUCB_horizon <- para_bandit_sim_bucb(data = data_list, rounds = 1000, 
                                          rate = "inverse_horizon",
                                          tau = tau_loc, epsilon = epsilon_loc, 
                                          alpha = tau_loc, beta = 1-tau_loc)
save(loc_BUCB_horizon, 
     file = paste0(current_path, "loc_BUCB_horizon.Rda"))
loc_comp_BUCB_horizon <- compare_to_ground_truth(mean_loc, 
                                                 loc_BUCB_horizon,
                                                 tau_loc,
                                                 epsilon_loc)$mean
save(loc_comp_BUCB_horizon, file = paste0(current_path, "loc_comp_BUCB_horizon.Rda"))

########################################################################

system.time(loc_KL_tau_horizon <- para_bandit_sim_KL(data = data_list, 
                                                      rounds = 1000, 
                                                      tau = tau_loc, 
                                                      epsilon = epsilon_loc, 
                                                      at_tau = TRUE, 
                                                      horizon = 1000^2))
save(loc_KL_tau_horizon, file = paste0(current_path, "loc_KL_tau_horizon.Rda"))
loc_comp_KL_tau_horizon <- compare_to_ground_truth(mean_loc, 
                                                    loc_KL_tau_horizon, 
                                                    tau_loc, epsilon_loc)$mean
save(loc_comp_KL_tau_horizon, file = paste0(current_path,
                                             "loc_comp_KL_tau_horizon.Rda"))
rm(loc_KL_tau_horizon)
gc()

########################################################################

# Do KL by comparing against tau directly

system.time(loc_KL <- para_bandit_sim_KL(data = data_list, rounds = 1000, 
                                         tau = tau_loc, epsilon = epsilon_loc, 
                                         at_tau = TRUE, horizon = 1000))
# user  system elapsed 
#6.932   4.007 585.704 
save(loc_KL, file = paste0(current_path, "loc_KL.Rda"))
load(file = paste0(current_path, "loc_KL.Rda"))
loc_comp_KL <- compare_to_ground_truth(mean_loc, loc_KL, tau_loc, 
                                       epsilon_loc)$mean
save(loc_comp_KL, file = paste0(current_path, "loc_comp_KL.Rda"))
rm(loc_KL)
gc()
########################################################################
# Do KL by comparing against tau and epsilon

system.time(loc_KL_not_tau_horizon1000 <- para_bandit_sim_KL(data = data_list, rounds = 1000, 
                                                 tau = tau_loc, epsilon = epsilon_loc, 
                                                 at_tau = FALSE, horizon = 1000))
save(loc_KL_not_tau_horizon1000, file = paste0(current_path, "loc_KL_not_tau_horizon1000.Rda"))
load(file = paste0(current_path, "loc_KL_not_tau_horizon1000.Rda"))
loc_comp_KL_not_tau_horizon1000 <- compare_to_ground_truth(mean_loc, 
                                                           loc_KL_not_tau_horizon1000, 
                                               tau_loc, epsilon_loc)$mean
save(loc_comp_KL_not_tau_horizon1000, file = paste0(current_path,
                                        "loc_comp_KL_not_tau_horizon1000.Rda"))
rm(loc_KL_not_tau_horizon1000)
gc()

########################################################################
# Do KL by comparing against tau and epsilon with 0/1 adjustment

system.time(loc_KL_horizon <- para_bandit_sim_KL(data = data_list, rounds = 1000, 
                                                 tau = tau_loc, epsilon = epsilon_loc, 
                                                 at_tau = FALSE, horizon = 10000))
save(loc_KL_horizon, file = paste0(current_path, "loc_KL_horizon.Rda"))
loc_comp_KL_horizon <- compare_to_ground_truth(mean_loc, loc_KL_horizon, 
                                               tau_loc, epsilon_loc)$mean
save(loc_comp_KL_horizon, file = paste0(current_path,
                                        "loc_comp_KL_horizon.Rda"))


########################################################################
# Do KL by comparing against tau and epsilon with 0/1 adjustment

system.time(loc_KL_horizon2 <- para_bandit_sim_KL(data = data_list, rounds = 1000, 
                                                  tau = tau_loc, epsilon = epsilon_loc, 
                                                  at_tau = FALSE, horizon = 1000^2))
save(loc_KL_horizon2, file = paste0(current_path, "loc_KL_horizon2.Rda"))
loc_comp_KL_horizon2 <- compare_to_ground_truth(mean_loc, loc_KL_horizon2, 
                                                tau_loc, epsilon_loc)$mean
save(loc_comp_KL_horizon2, file = paste0(current_path,
                                         "loc_comp_KL_horizon2.Rda"))

########################################################################
# Do KL by comparing against tau and epsilon with 0/1 adjustment

system.time(loc_KL_horizon1 <- para_bandit_sim_KL(data = data_list, rounds = 1000, 
                                                  tau = tau_loc, epsilon = epsilon_loc, 
                                                  at_tau = FALSE, horizon = 1000))
save(loc_KL_horizon1, file = paste0(current_path, "loc_KL_horizon1.Rda"))
load(file = paste0(current_path, "loc_KL_horizon1.Rda"))
loc_comp_KL_horizon1 <- compare_to_ground_truth(mean_loc, loc_KL_horizon1, 
                                                tau_loc, epsilon_loc)$mean
save(loc_comp_KL_horizon1, file = paste0(current_path,
                                         "loc_comp_KL_horizon1.Rda"))

########################################################################
# Probability that arm is above or below tau+-epsilon

loc_PI <- para_bandit_sim_PI(data = data_list, rounds = 1000,
                             tau = tau_loc, epsilon = epsilon_loc,
                             alpha = tau_loc, beta = 1 - tau_loc)
save(loc_PI, file = paste0(current_path, "loc_PI.Rda"))
load(file = paste0(current_path, "loc_PI.Rda"))
loc_comp_PI <- compare_to_ground_truth(mean_loc, loc_PI, 
                                       tau_loc, epsilon_loc)$mean
save(loc_comp_PI, file = paste0(current_path, "loc_comp_PI.Rda"))
rm(loc_PI)
gc()
########################################################################

loc_TTS <- para_bandit_sim_TTS(data = data_list, rounds = 1000,
                               tau = tau_loc, epsilon = epsilon_loc,
                               alpha = tau_loc, beta = 1 - tau_loc)
save(loc_TTS, file = paste0(current_path, "loc_TTS.Rda"))
load(file = paste0(current_path, "loc_TTS.Rda"))
loc_comp_TTS <- compare_to_ground_truth(mean_loc, loc_TTS, 
                                        tau_loc, epsilon_loc)$mean
save(loc_comp_TTS, file = paste0(current_path, "loc_comp_TTS.Rda"))
rm(loc_TTS)
gc()
########################################################################

#system.time(loc_BETA <- para_bandit_sim_BETA(data = data_list, rounds = 1000,
#                                             tau = tau_loc, epsilon = epsilon_loc,
#                                             alpha = tau_loc, beta = 1 - tau_loc))
#save(loc_BETA, file = paste0(current_path, "loc_BETA.Rda"))
#loc_comp_BETA <- compare_to_ground_truth(mean_loc, loc_BETA, 
#                                         tau_loc, epsilon_loc)$mean
#save(loc_comp_BETA, file = paste0(current_path, "loc_comp_BETA.Rda"))

########################################################################

loc_BUCB <- para_bandit_sim_bucb(data = data_list, rounds = 1000, 
                                 rate = "inverse",
                                 tau = tau_loc, epsilon = epsilon_loc, 
                                 alpha = tau_loc, beta = 1-tau_loc)
save(loc_BUCB, file = paste0(current_path, "loc_BUCB.Rda"))
load(file = paste0(current_path, "loc_BUCB.Rda"))
loc_comp_BUCB <- compare_to_ground_truth(mean_loc, loc_BUCB, 
                                         tau_loc, epsilon_loc)$mean
save(loc_comp_BUCB, file = paste0(current_path, "loc_comp_BUCB.Rda"))
rm(loc_BUCB)
gc()
########################################################################

loc_BUCB_squared <- para_bandit_sim_bucb(data = data_list, rounds = 1000, 
                                         rate = "inverse_squared",
                                         tau = tau_loc, epsilon = epsilon_loc, 
                                         alpha = tau_loc, beta = 1-tau_loc)
save(loc_BUCB_squared, file = paste0(current_path, "loc_BUCB_squared.Rda"))
load(file = paste0(current_path, "loc_BUCB_squared.Rda"))
loc_comp_BUCB_squared <- compare_to_ground_truth(mean_loc, loc_BUCB_squared, 
                                                 tau_loc, epsilon_loc)$mean
save(loc_comp_BUCB_squared, file = paste0(current_path, 
                                          "loc_comp_BUCB_squared.Rda"))
rm(loc_BUCB_squared)
gc()
########################################################################

loc_BUCB_power5 <- para_bandit_sim_bucb(data = data_list, rounds = 1000, 
                                        rate = "inverse_power5",
                                        tau = tau_loc, epsilon = epsilon_loc, 
                                        alpha = tau_loc, beta = 1-tau_loc)
save(loc_BUCB_power5, file = paste0(current_path, "loc_BUCB_power5.Rda"))
load(file = paste0(current_path, "loc_BUCB_power5.Rda"))
loc_comp_BUCB_power5 <- compare_to_ground_truth(mean_loc, loc_BUCB_power5, 
                                                tau_loc, epsilon_loc)$mean
save(loc_comp_BUCB_power5, file = paste0(current_path, 
                                         "loc_comp_BUCB_power5.Rda"))
rm(loc_BUCB_power5)
gc()
########################################################################

load(paste0(current_path, "loc_comp_APT.Rda"))
load(paste0(current_path, "loc_comp_AugUCB.Rda"))
load(paste0(current_path, "loc_comp_UNIFORM.Rda"))
load(paste0(current_path, "loc_comp_PI.Rda"))
load(paste0(current_path, "loc_comp_TTS.Rda"))
load(paste0(current_path, "loc_comp_BUCB.Rda"))
load(paste0(current_path, "loc_comp_BUCB_horizon.Rda"))
load(paste0(current_path, "loc_comp_BUCB_squared.Rda"))
load(paste0(current_path, "loc_comp_BUCB_power5.Rda"))
load(paste0(current_path, "loc_comp_KL.Rda"))
load(paste0(current_path, "loc_comp_KL_tau_horizon.Rda"))
load(paste0(current_path, "loc_comp_KL_horizon1.Rda"))
load(paste0(current_path, "loc_comp_KL_not_tau_horizon1000.Rda"))

plot(c(0,600), c(0, -9), type = "n")
lines(log(loc_comp_APT), col = "red")
lines(log(loc_comp_AugUCB), col = "grey")
lines(log(loc_comp_UNIFORM), col = "black")
lines(log(loc_comp_PI), col = "darkgreen")
lines(log(loc_comp_TTS), col = "lightgreen")
lines(log(loc_comp_BUCB), col = "orange")
lines(log(loc_comp_BUCB_horizon), col = "darkorange")
lines(log(loc_comp_BUCB_squared), col = "blue")
lines(log(loc_comp_BUCB_power5), col = "darkred")
lines(log(loc_comp_KL), col = "blue", lty = 2)
lines(log(loc_comp_KL_tau_horizon), col = "black", lty = 2)
lines(log(loc_comp_KL_horizon1), col = "black", lty = 2)
lines(log(loc_comp_KL_not_tau_horizon1000), col = "red", lty = 2)
#lines(log(loc_comp_KL), col = "lightblue")
#lines(log(loc_comp_KL_not_tau), col = "blue")
#lines(log(loc_comp_KL_horizon), col = "darkblue")
#lines(log(loc_comp_KL_horizon2), col = "darkred")
#lines(log(loc_comp_KL_horizon1), col = "red")

########################################################################

library(plyr)
loc_UNIFORM_as <- colMeans(ldply(loc_UNIFORM, function(x) table(x$arm_sequence)))
loc_APT_as <- colMeans(ldply(loc_APT, function(x) table(x$arm_sequence)))
loc_PI_as <- colMeans(ldply(loc_PI, function(x) table(x$arm_sequence)))
loc_TTS_as <- colMeans(ldply(loc_TTS, function(x) table(x$arm_sequence)))
loc_BUCB_as <- colMeans(ldply(loc_BUCB, function(x) table(x$arm_sequence)))
loc_KL_as <- colMeans(ldply(loc_KL, function(x) table(x$arm_sequence)))
loc_KL_not_tau_as <- colMeans(ldply(loc_KL_not_tau, function(x) table(x$arm_sequence)))
loc_KL_horizon_as <- colMeans(ldply(loc_KL_horizon, function(x) table(x$arm_sequence)))

round(data.frame(loc_UNIFORM_as, loc_APT_as, loc_PI_as, loc_TTS_as, loc_BUCB_as, 
                 loc_KL_as, loc_KL_not_tau_as, loc_KL_horizon_as))
