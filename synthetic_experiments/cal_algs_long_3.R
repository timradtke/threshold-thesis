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

plot(mean_loc3, rep(1,10))
abline(v=tau_loc3)
abline(v=tau_loc3+epsilon_loc3, lty=2)
abline(v=tau_loc3-epsilon_loc3, lty=2)

data_list3 <- list()
set.seed(1024)
for(j in 1:5000) {
  curr_data <- data.frame(rep(NA, times = 2500))
  for(i in 1:length(mean_loc3)) {
    curr_data[[i]] <- as.numeric(purrr::rbernoulli(2500, p  = mean_loc3[i]))
  }
  names(curr_data) <- paste0("V", rep(1:length(mean_loc3)))
  data_list3[[j]] <- curr_data
}

########################################################################
# Standard APT Algorithm
# a seed is built in
system.time(loc3_APT <- para_bandit_sim_APT(data = data_list3, rounds = 2500, 
                                            tau = tau_loc3, 
                                            epsilon = epsilon_loc3))
# user  system elapsed 
#12.497    9.065 1227.363
save(loc3_APT, file = paste0(current_path, "loc3_APT.Rda"))
#load(file = paste0(current_path, "loc3_APT.Rda"))
loc3_comp_APT <- compare_to_ground_truth(mean_loc3, loc3_APT, tau_loc3, 
                                         epsilon_loc3)$mean
save(loc3_comp_APT, file = paste0(current_path, "loc3_comp_APT.Rda"))
rm(loc3_APT)
gc()

########################################################################

loc3_BUCB <- para_bandit_sim_bucb(data = data_list3, rounds = 2500, 
                                  rate = "inverse",
                                  tau = tau_loc3, epsilon = epsilon_loc3, 
                                  alpha = tau_loc3, beta = 1-tau_loc3)
save(loc3_BUCB, file = paste0(current_path, "loc3_BUCB.Rda"))
#load(file = paste0(current_path, "loc3_BUCB.Rda"))
loc3_comp_BUCB <- compare_to_ground_truth(mean_loc3, loc3_BUCB, 
                                          tau_loc3, epsilon_loc3)$mean
save(loc3_comp_BUCB, file = paste0(current_path, "loc3_comp_BUCB.Rda"))
rm(loc3_BUCB)
gc()

########################################################################

loc3_BUCB_horizon <- para_bandit_sim_bucb(data = data_list3, rounds = 2500, 
                                          rate = "inverse_horizon",
                                          tau = tau_loc3, epsilon = epsilon_loc3, 
                                          alpha = tau_loc3, beta = 1-tau_loc3)
save(loc3_BUCB_horizon, 
     file = paste0(current_path, "loc3_BUCB_horizon.Rda"))
loc3_comp_BUCB_horizon <- compare_to_ground_truth(mean_loc3, 
                                                  loc3_BUCB_horizon,
                                                  tau_loc3,
                                                  epsilon_loc3)$mean
save(loc3_comp_BUCB_horizon, file = paste0(current_path, "loc3_comp_BUCB_horizon.Rda"))
rm(loc3_BUCB_horizon)
gc()

########################################################################

loc3_TTS <- para_bandit_sim_TTS(data = data_list3, rounds = 2500,
                                tau = tau_loc3, epsilon = epsilon_loc3,
                                alpha = tau_loc3, beta = 1 - tau_loc3)
save(loc3_TTS, file = paste0(current_path, "loc3_TTS.Rda"))
loc3_comp_TTS <- compare_to_ground_truth(mean_loc3, loc3_TTS, 
                                         tau_loc3, epsilon_loc3)$mean
save(loc3_comp_TTS, file = paste0(current_path, "loc3_comp_TTS.Rda"))

########################################################################
# Do KL by comparing against tau with 0/1 adjustment

system.time(loc3_KL_tau_horizon <- para_bandit_sim_KL(data = data_list3, 
                                                      rounds = 2500, 
                                                      tau = tau_loc3, 
                                                      epsilon = epsilon_loc3, 
                                                      at_tau = TRUE, 
                                                      horizon = 2500^2))
save(loc3_KL_tau_horizon, file = paste0(current_path, "loc3_KL_tau_horizon.Rda"))
loc3_comp_KL_tau_horizon <- compare_to_ground_truth(mean_loc3, 
                                                    loc3_KL_tau_horizon, 
                                                    tau_loc3, epsilon_loc3)$mean
save(loc3_comp_KL_tau_horizon, file = paste0(current_path,
                                            "loc3_comp_KL_tau_horizon.Rda"))
rm(loc3_KL_tau_horizon)
gc()

########################################################################

system.time(loc3_UNIFORM <- para_bandit_sim_uniform(data = data_list3, 
                                                    rounds = 2500))
#   user  system elapsed 
# 11.406   7.720 781.403 
save(loc3_UNIFORM, file = paste0(current_path, "loc3_UNIFORM.Rda"))
load(file = paste0(current_path, "loc3_UNIFORM.Rda"))
loc3_comp_UNIFORM <- compare_to_ground_truth(mean_loc3, loc3_UNIFORM, tau_loc3, 
                                             epsilon_loc3)$mean
save(loc3_comp_UNIFORM, file = paste0(current_path, "loc3_comp_UNIFORM.Rda"))
rm(loc3_UNIFORM)
gc()

########################################################################
# Probability that arm is above or below tau+-epsilon

loc3_PI <- para_bandit_sim_PI(data = data_list3, rounds = 2500,
                              tau = tau_loc3, epsilon = epsilon_loc3,
                              alpha = tau_loc3, beta = 1 - tau_loc3)
save(loc3_PI, file = paste0(current_path, "loc3_PI.Rda"))
loc3_comp_PI <- compare_to_ground_truth(mean_loc3, loc3_PI, 
                                        tau_loc3, epsilon_loc3)$mean
save(loc3_comp_PI, file = paste0(current_path, "loc3_comp_PI.Rda"))

########################################################################
# Standard AugUCB
system.time(loc3_AugUCB <- para_bandit_sim_AugUCB(data = data_list3, rounds = 2500, 
                                                  tau = tau_loc3))
save(loc3_AugUCB, file = paste0(current_path, "loc3_AugUCB.Rda"))
load(file = paste0(current_path, "loc3_AugUCB.Rda"))
loc3_comp_AugUCB <- compare_to_ground_truth(mean_loc3, loc3_AugUCB, tau_loc3, 
                                            epsilon_loc3)$mean
save(loc3_comp_AugUCB, file = paste0(current_path, "loc3_comp_AugUCB.Rda"))
rm(loc3_AugUCB)
gc()

########################################################################
# Do KL by comparing against tau and epsilon with 0/1 adjustment

system.time(loc3_KL_horizon <- para_bandit_sim_KL(data = data_list3, rounds = 2500, 
                                                  tau = tau_loc3, epsilon = epsilon_loc3, 
                                                  at_tau = FALSE, horizon = 2500))
save(loc3_KL_horizon, file = paste0(current_path, "loc3_KL_horizon.Rda"))
load(file = paste0(current_path, "loc3_KL_horizon.Rda"))
loc3_comp_KL_horizon <- compare_to_ground_truth(mean_loc3, loc3_KL_horizon, 
                                                tau_loc3, epsilon_loc3)$mean
save(loc3_comp_KL_horizon, file = paste0(current_path,
                                         "loc3_comp_KL_horizon.Rda"))
rm(loc3_KL_horizon)
gc()

########################################################################
########################################################################
########################################################################
########################################################################
########################################################################












########################################################################
# Standard Uniform

system.time(loc3_UNIFORM <- para_bandit_sim_uniform(data = data_list3, 
                                                    rounds = 1000))
# user  system elapsed 
# 6.268   3.400 211.002 
save(loc3_UNIFORM, file = paste0(current_path, "loc3_UNIFORM.Rda"))
load(file = paste0(current_path, "loc3_UNIFORM.Rda"))
loc3_comp_UNIFORM <- compare_to_ground_truth(mean_loc3, loc3_UNIFORM, tau_loc3, 
                                             epsilon_loc3)$mean
save(loc3_comp_UNIFORM, file = paste0(current_path, "loc3_comp_UNIFORM.Rda"))
rm(loc3_UNIFORM)
gc()
########################################################################
# Do KL by comparing against tau directly

system.time(loc3_KL <- para_bandit_sim_KL(data = data_list3, rounds = 1000, 
                                          tau = tau_loc3, epsilon = epsilon_loc3, 
                                          at_tau = TRUE, horizon = 1000))
# user  system elapsed 
#6.824   4.335 391.240
save(loc3_KL, file = paste0(current_path, "loc3_KL.Rda"))
load(file = paste0(current_path, "loc3_KL.Rda"))
loc3_comp_KL <- compare_to_ground_truth(mean_loc3, loc3_KL, tau_loc3, 
                                        epsilon_loc3)$mean
save(loc3_comp_KL, file = paste0(current_path, "loc3_comp_KL.Rda"))
rm(loc3_KL)
gc()
########################################################################
# Do KL by comparing against tau and epsilon

system.time(loc3_KL_not_tau <- para_bandit_sim_KL(data = data_list3, rounds = 1000, 
                                                  tau = tau_loc3, epsilon = epsilon_loc3, 
                                                  at_tau = FALSE))
save(loc3_KL_not_tau, file = paste0(current_path, "loc3_KL_not_tau.Rda"))
load(file = paste0(current_path, "loc3_KL_not_tau.Rda"))
loc3_comp_KL_not_tau <- compare_to_ground_truth(mean_loc3, loc3_KL_not_tau, 
                                                tau_loc3, epsilon_loc3)$mean
save(loc3_comp_KL_not_tau, file = paste0(current_path,
                                         "loc3_comp_KL_not_tau.Rda"))
rm(loc3_KL_not_tau)
gc()

########################################################################
# Do KL by comparing against tau and epsilon with 0/1 adjustment

system.time(loc3_KL_horizon2 <- para_bandit_sim_KL(data = data_list3, rounds = 1000, 
                                                   tau = tau_loc3, epsilon = epsilon_loc3, 
                                                   at_tau = FALSE, horizon = 1000^2))
save(loc3_KL_horizon2, file = paste0(current_path, "loc3_KL_horizon2.Rda"))
loc3_comp_KL_horizon2 <- compare_to_ground_truth(mean_loc3, loc3_KL_horizon2, 
                                                 tau_loc3, epsilon_loc3)$mean
save(loc3_comp_KL_horizon2, file = paste0(current_path,
                                          "loc3_comp_KL_horizon2.Rda"))

########################################################################
# Do KL by comparing against tau and epsilon with 0/1 adjustment

system.time(loc3_KL_horizon1 <- para_bandit_sim_KL(data = data_list3, rounds = 1000, 
                                                   tau = tau_loc3, epsilon = epsilon_loc3, 
                                                   at_tau = FALSE, horizon = 1000))
save(loc3_KL_horizon1, file = paste0(current_path, "loc3_KL_horizon1.Rda"))
loc3_comp_KL_horizon1 <- compare_to_ground_truth(mean_loc3, loc3_KL_horizon1, 
                                                 tau_loc3, epsilon_loc3)$mean
save(loc3_comp_KL_horizon1, file = paste0(current_path,
                                          "loc3_comp_KL_horizon1.Rda"))

########################################################################
# Probability that arm is above or below tau+-epsilon

loc3_PI <- para_bandit_sim_PI(data = data_list3, rounds = 1000,
                              tau = tau_loc3, epsilon = epsilon_loc3,
                              alpha = tau_loc3, beta = 1 - tau_loc3)
save(loc3_PI, file = paste0(current_path, "loc3_PI.Rda"))
loc3_comp_PI <- compare_to_ground_truth(mean_loc3, loc3_PI, 
                                        tau_loc3, epsilon_loc3)$mean
save(loc3_comp_PI, file = paste0(current_path, "loc3_comp_PI.Rda"))

########################################################################
# Probability that arm is above or below tau+-epsilon

loc3_PI_tadj <- para_bandit_sim_PI(data = data_list3, rounds = 1000,
                                   tau = tau_loc3, epsilon = epsilon_loc3,
                                   alpha = tau_loc3, beta = 1 - tau_loc3, tadj = TRUE)
save(loc3_PI_tadj, file = paste0(current_path, "loc3_PI_tadj.Rda"))
loc3_comp_PI_tadj <- compare_to_ground_truth(mean_loc3, loc3_PI_tadj, 
                                             tau_loc3, epsilon_loc3)$mean
save(loc3_comp_PI_tadj, file = paste0(current_path, "loc3_comp_PI_tadj.Rda"))

########################################################################

loc3_TTS <- para_bandit_sim_TTS(data = data_list3, rounds = 1000,
                                tau = tau_loc3, epsilon = epsilon_loc3,
                                alpha = tau_loc3, beta = 1 - tau_loc3)
save(loc3_TTS, file = paste0(current_path, "loc3_TTS.Rda"))
loc3_comp_TTS <- compare_to_ground_truth(mean_loc3, loc3_TTS, 
                                         tau_loc3, epsilon_loc3)$mean
save(loc3_comp_TTS, file = paste0(current_path, "loc3_comp_TTS.Rda"))

########################################################################

#system.time(loc3_BETA <- para_bandit_sim_BETA(data = data_list3, rounds = 1000,
#                                              tau = tau_loc3, epsilon = epsilon_loc3,
#                                              alpha = tau_loc3, beta = 1 - tau_loc3))
#save(loc3_BETA, file = paste0(current_path, "loc3_BETA.Rda"))
#loc3_comp_BETA <- compare_to_ground_truth(mean_loc3, loc3_BETA, 
#                                          tau_loc3, epsilon_loc3)$mean
#save(loc3_comp_BETA, file = paste0(current_path, "loc3_comp_BETA.Rda"))


########################################################################

loc3_BUCB_squared <- para_bandit_sim_bucb(data = data_list3, rounds = 1000, 
                                          rate = "inverse_squared",
                                          tau = tau_loc3, epsilon = epsilon_loc3, 
                                          alpha = tau_loc3, beta = 1-tau_loc3)
save(loc3_BUCB_squared, file = paste0(current_path, "loc3_BUCB_squared.Rda"))
load(file = paste0(current_path, "loc3_BUCB_squared.Rda"))
loc3_comp_BUCB_squared <- compare_to_ground_truth(mean_loc3, loc3_BUCB_squared, 
                                                  tau_loc3, epsilon_loc3)$mean
save(loc3_comp_BUCB_squared, file = paste0(current_path, "loc3_comp_BUCB_squared.Rda"))
rm(loc3_BUCB_squared)
gc()
########################################################################

loc3_BUCB_cubic <- para_bandit_sim_bucb(data = data_list3, rounds = 1000, 
                                        rate = "inverse_cubic",
                                        tau = tau_loc3, epsilon = epsilon_loc3, 
                                        alpha = tau_loc3, beta = 1-tau_loc3)
save(loc3_BUCB_cubic, file = paste0(current_path, "loc3_BUCB_cubic.Rda"))
loc3_comp_BUCB_cubic <- compare_to_ground_truth(mean_loc3, loc3_BUCB_cubic, 
                                                tau_loc3, epsilon_loc3)$mean
save(loc3_comp_BUCB_cubic, file = paste0(current_path, "loc3_comp_BUCB_cubic.Rda"))

########################################################################

loc3_BUCB_power5 <- para_bandit_sim_bucb(data = data_list3, rounds = 1000, 
                                         rate = "inverse_power5",
                                         tau = tau_loc3, epsilon = epsilon_loc3, 
                                         alpha = tau_loc3, beta = 1-tau_loc3)
save(loc3_BUCB_power5, file = paste0(current_path, "loc3_BUCB_power5.Rda"))
loc3_comp_BUCB_power5 <- compare_to_ground_truth(mean_loc3, loc3_BUCB_power5, 
                                                 tau_loc3, epsilon_loc3)$mean
save(loc3_comp_BUCB_power5, file = paste0(current_path, "loc3_comp_BUCB_power5.Rda"))


########################################################################

loc3_BUCB_squared_no_eps <- para_bandit_sim_bucb(data = data_list3, rounds = 1000, 
                                                 rate = "inverse_squared",
                                                 tau = tau_loc3, epsilon = epsilon_loc3, 
                                                 alpha = tau_loc3, beta = 1-tau_loc3,
                                                 with_epsilon = FALSE)
save(loc3_BUCB_squared_no_eps, 
     file = paste0(current_path, "loc3_BUCB_squared_no_eps.Rda"))
loc3_comp_BUCB_squared_no_eps <- compare_to_ground_truth(mean_loc3, 
                                                         loc3_BUCB_squared_no_eps,
                                                         tau_loc3,
                                                         epsilon_loc3)$mean
save(loc3_comp_BUCB_squared_no_eps, file = paste0(current_path, "loc3_comp_BUCB_squared_no_eps.Rda"))

########################################################################

loc3_BUCB_horizon <- para_bandit_sim_bucb(data = data_list3, rounds = 1000, 
                                          rate = "inverse_horizon",
                                          tau = tau_loc3, epsilon = epsilon_loc3, 
                                          alpha = tau_loc3, beta = 1-tau_loc3,
                                          with_epsilon = FALSE, const = 5)
save(loc3_BUCB_horizon, 
     file = paste0(current_path, "loc3_BUCB_horizon.Rda"))
loc3_comp_BUCB_horizon <- compare_to_ground_truth(mean_loc3, 
                                                  loc3_BUCB_horizon,
                                                  tau_loc3,
                                                  epsilon_loc3)$mean
save(loc3_comp_BUCB_horizon, file = paste0(current_path, "loc3_comp_BUCB_horizon.Rda"))

########################################################################

loc3_BUCB_horizon5 <- para_bandit_sim_bucb(data = data_list3, rounds = 1000, 
                                           rate = "inverse_horizon_c",
                                           tau = tau_loc3, epsilon = epsilon_loc3, 
                                           alpha = tau_loc3, beta = 1-tau_loc3,
                                           with_epsilon = FALSE, const = 5)
save(loc3_BUCB_horizon5, 
     file = paste0(current_path, "loc3_BUCB_horizon5.Rda"))
load(file = paste0(current_path, "loc3_BUCB_horizon5.Rda"))
loc3_comp_BUCB_horizon5 <- compare_to_ground_truth(mean_loc3, 
                                                   loc3_BUCB_horizon5,
                                                   tau_loc3,
                                                   epsilon_loc3)$mean
save(loc3_comp_BUCB_horizon5, file = paste0(current_path, "loc3_comp_BUCB_horizon5.Rda"))
rm(loc3_BUCB_horizon5)
gc()
########################################################################

load(paste0(current_path, "loc3_comp_APT.Rda"))
load(paste0(current_path, "loc3_comp_UNIFORM.Rda"))
load(paste0(current_path, "loc3_comp_BUCB.Rda"))
load(paste0(current_path, "loc3_comp_BUCB_horizon.Rda"))
load(paste0(current_path, "loc3_comp_TTS.Rda"))
load(paste0(current_path, "loc3_comp_PI.Rda"))
load(paste0(current_path, "loc3_comp_AugUCB.Rda"))
load(paste0(current_path, "loc3_comp_KL_tau_horizon.Rda"))
load(paste0(current_path, "loc3_comp_KL_horizon.Rda"))

plot(c(0,2500), c(0, -6), type = "n")
lines(log(loc3_comp_APT), col = "blue")
lines(log(loc3_comp_BUCB), col = "green")
lines(log(loc3_comp_BUCB_horizon), col = "pink")
lines(log(loc3_comp_AugUCB), col = "grey")
lines(log(loc3_comp_UNIFORM), col = "black")
lines(log(loc3_comp_TTS), col = "lightgreen")
lines(log(loc3_comp_PI), col = "darkgreen")
lines(log(loc3_comp_KL_tau_horizon), col = "darkred")
lines(log(loc3_comp_KL_horizon), col = "red")

########################################################################

library(plyr)
loc3_UNIFORM_as <- colMeans(ldply(loc3_UNIFORM, function(x) table(x$arm_sequence)))
loc3_APT_as <- colMeans(ldply(loc3_APT, function(x) table(x$arm_sequence)))
loc3_PI_as <- colMeans(ldply(loc3_PI, function(x) table(x$arm_sequence)))
loc3_PI_tadj_as <- colMeans(ldply(loc3_PI_tadj, function(x) table(x$arm_sequence)))
loc3_TTS_as <- colMeans(ldply(loc3_TTS, function(x) table(x$arm_sequence)))
loc3_BUCB_as <- colMeans(ldply(loc3_BUCB, function(x) table(x$arm_sequence)))
loc3_BUCB_squared_as <- colMeans(ldply(loc3_BUCB_squared, function(x) table(x$arm_sequence)))
loc3_BUCB_cubic_as <- colMeans(ldply(loc3_BUCB_cubic, function(x) table(x$arm_sequence)))
loc3_KL_as <- colMeans(ldply(loc3_KL, function(x) table(x$arm_sequence)))
loc3_KL_not_tau_as <- colMeans(ldply(loc3_KL_not_tau, function(x) table(x$arm_sequence)))
loc3_KL_horizon1_as <- colMeans(ldply(loc3_KL_horizon1, function(x) table(x$arm_sequence)))

round(data.frame(loc3_UNIFORM_as, loc3_APT_as, 
                 loc3_BUCB_as, loc3_BUCB_squared_as, loc3_BUCB_cubic_as,
                 loc3_PI_as, loc3_PI_tadj_as, loc3_TTS_as,
                 loc3_KL_as, loc3_KL_horizon1_as))


loc3_UNIFORM_as_1 <- colMeans(ldply(loc3_UNIFORM, function(x) table(x$arm_sequence[1:500])))
loc3_APT_as_1 <- colMeans(ldply(loc3_APT, function(x) table(x$arm_sequence[1:500])))
loc3_PI_as_1 <- colMeans(ldply(loc3_PI, function(x) table(x$arm_sequence[1:500])))
loc3_PI_tadj_as_1 <- colMeans(ldply(loc3_PI_tadj, function(x) table(x$arm_sequence[1:500])))
loc3_TTS_as_1 <- colMeans(ldply(loc3_TTS, function(x) table(x$arm_sequence[1:500])))
loc3_BUCB_as_1 <- colMeans(ldply(loc3_BUCB, function(x) table(x$arm_sequence[1:500])))
loc3_KL_as_1 <- colMeans(ldply(loc3_KL, function(x) table(x$arm_sequence[1:500])))
loc3_KL_not_tau_as_1 <- colMeans(ldply(loc3_KL_not_tau, function(x) table(x$arm_sequence[1:500])))
loc3_KL_horizon1_as_1 <- colMeans(ldply(loc3_KL_horizon1, function(x) table(x$arm_sequence[1:500])))

round(data.frame(loc3_UNIFORM_as_1, loc3_APT_as_1, loc3_BUCB_as_1,
                 loc3_PI_as_1, loc3_PI_tadj_as_1, loc3_TTS_as_1,
                 loc3_KL_as_1, loc3_KL_horizon1_as_1))/500


loc3_UNIFORM_as_2 <- colMeans(ldply(loc3_UNIFORM, function(x) table(x$arm_sequence[501:1000])))
loc3_APT_as_2 <- colMeans(ldply(loc3_APT, function(x) table(x$arm_sequence[501:1000])))
loc3_PI_as_2 <- colMeans(ldply(loc3_PI, function(x) table(x$arm_sequence[501:1000])))
loc3_PI_tadj_as_2 <- colMeans(ldply(loc3_PI_tadj, function(x) table(x$arm_sequence[501:1000])))
loc3_TTS_as_2 <- colMeans(ldply(loc3_TTS, function(x) table(x$arm_sequence[501:1000])))
loc3_BUCB_as_2 <- colMeans(ldply(loc3_BUCB, function(x) table(x$arm_sequence[501:1000])))
loc3_KL_as_2 <- colMeans(ldply(loc3_KL, function(x) table(x$arm_sequence[501:1000])))
loc3_KL_not_tau_as_2 <- colMeans(ldply(loc3_KL_not_tau, function(x) table(x$arm_sequence[501:1000])))
loc3_KL_horizon1_as_2 <- colMeans(ldply(loc3_KL_horizon1, function(x) table(x$arm_sequence[501:1000])))

round(data.frame(loc3_UNIFORM_as, loc3_APT_as, loc3_BUCB_as,
                 loc3_PI_as, loc3_PI_tadj_as, loc3_TTS_as,
                 loc3_KL_as, loc3_KL_horizon1_as))

########################################################################

plot(loc3_BUCB_squared[[1]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc3_BUCB_squared[[2]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc3_BUCB_squared[[3]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc3_BUCB_squared[[4]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc3_BUCB_squared[[5]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc3_BUCB_squared[[6]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc3_BUCB_squared[[7]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc3_BUCB_squared[[8]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc3_BUCB_squared[[9]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc3_BUCB_squared[[10]]$arm_sequence, pch = 19, cex = 0.3)

plot(loc3_APT[[1]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc3_APT[[2]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc3_APT[[3]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc3_APT[[4]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc3_APT[[5]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc3_APT[[6]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc3_APT[[7]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc3_APT[[8]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc3_APT[[9]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc3_APT[[10]]$arm_sequence, pch = 19, cex = 0.3)

arm_seq_res_APT <- data.frame(t(laply(loc3_APT, function(x) x$arm_sequence)))
library(tidyr)
library(dplyr)
#arm_seq_table <- arm_seq_res %>% tbl_df() %>% mutate(index = 1:1000) %>%
#  gather(key = iter, value = arm, -index) %>%
#  group_by(index, arm) %>% summarize(count = n())

library(ggplot2)
arm_seq_res_APT %>% tbl_df() %>% mutate(index = 1:1000) %>%
  gather(key = iter, value = arm, -index) %>%
  filter(index > 10) %>%
  ggplot(aes(index, arm)) + geom_count() +
  scale_size_area(max_size = 10)
#ggplot(arm_seq_table[11:9910,], aes(index, arm)) + geom_tile(aes(fill = count))

arm_seq_res_BUCB <- data.frame(t(laply(loc3_BUCB, function(x) x$arm_sequence)))
arm_seq_res_BUCB %>% tbl_df() %>% mutate(index = 1:1000) %>%
  gather(key = iter, value = arm, -index) %>%
  filter(index > 10) %>%
  ggplot(aes(index, arm)) + geom_count() +
  scale_size_area(max_size = 10)

arm_seq_res_BUCB_squared <- data.frame(t(laply(loc3_BUCB_squared, function(x) x$arm_sequence)))
arm_seq_res_BUCB_squared %>% tbl_df() %>% mutate(index = 1:1000) %>%
  gather(key = iter, value = arm, -index) %>%
  filter(index > 10) %>%
  ggplot(aes(index, arm)) + geom_count() +
  scale_size_area(max_size = 10)
