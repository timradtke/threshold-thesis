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

mean_loc5 <- c(0.0005, 0.001,
               0.05, 0.055,
               0.065, 0.075,
               0.085, 0.09,
               0.1, 0.1)
tau_loc5 <- 0.07
epsilon_loc5 <- 0.01

plot(mean_loc5, rep(1,10))
abline(v=tau_loc5)
abline(v=tau_loc5+epsilon_loc5, lty=2)
abline(v=tau_loc5-epsilon_loc5, lty=2)

data_list5 <- list()
set.seed(256)
for(j in 1:5000) {
  curr_data <- data.frame(rep(NA, times = 2500))
  for(i in 1:length(mean_loc5)) {
    curr_data[[i]] <- as.numeric(purrr::rbernoulli(2500, p  = mean_loc5[i]))
  }
  names(curr_data) <- paste0("V", rep(1:length(mean_loc5)))
  data_list5[[j]] <- curr_data
}

########################################################################
# Standard APT Algorithm
# a seed is built in
system.time(loc5_APT <- para_bandit_sim_APT(data = data_list5, rounds = 2500, 
                                            tau = tau_loc5, epsilon = epsilon_loc5))
# user  system elapsed 
#6.508   3.544 300.383
save(loc5_APT, file = paste0(current_path, "loc5_APT.Rda"))
load(file = paste0(current_path, "loc5_APT.Rda"))
loc5_comp_APT <- compare_to_ground_truth(mean_loc5, loc5_APT, tau_loc5, 
                                         epsilon_loc5)$mean
save(loc5_comp_APT, file = paste0(current_path, "loc5_comp_APT.Rda"))
load(file = paste0(current_path, "loc5_comp_APT.Rda"))
rm(loc5_APT)
gc()

########################################################################

loc5_BUCB <- para_bandit_sim_bucb(data = data_list5, rounds = 2000, 
                                  rate = "inverse",
                                  tau = tau_loc5, epsilon = epsilon_loc5, 
                                  alpha = tau_loc5, beta = 1-tau_loc5)
save(loc5_BUCB, file = paste0(current_path, "loc5_BUCB.Rda"))
load(file = paste0(current_path, "loc5_BUCB.Rda"))
loc5_comp_BUCB <- compare_to_ground_truth(mean_loc5, loc5_BUCB, 
                                          tau_loc5, epsilon_loc5)$mean
save(loc5_comp_BUCB, file = paste0(current_path, "loc5_comp_BUCB.Rda"))
load(file = paste0(current_path, "loc5_comp_BUCB.Rda"))

########################################################################

loc5_BUCB_horizon <- para_bandit_sim_bucb(data = data_list5, rounds = 2500, 
                                          rate = "inverse_horizon",
                                          tau = tau_loc5, epsilon = epsilon_loc5, 
                                          alpha = tau_loc5, beta = 1-tau_loc5)
save(loc5_BUCB_horizon, 
     file = paste0(current_path, "loc5_BUCB_horizon.Rda"))
loc5_comp_BUCB_horizon <- compare_to_ground_truth(mean_loc5, 
                                                  loc5_BUCB_horizon,
                                                  tau_loc5,
                                                  epsilon_loc5)$mean
save(loc5_comp_BUCB_horizon, file = paste0(current_path, "loc5_comp_BUCB_horizon.Rda"))
load(file = paste0(current_path, "loc5_comp_BUCB_horizon.Rda"))
########################################################################

loc5_BUCB_squared <- para_bandit_sim_bucb(data = data_list5, rounds = 2500, 
                                          rate = "inverse_squared",
                                          tau = tau_loc5, epsilon = epsilon_loc5, 
                                          alpha = tau_loc5, beta = 1-tau_loc5)
save(loc5_BUCB_squared, file = paste0(current_path, "loc5_BUCB_squared.Rda"))
load(file = paste0(current_path, "loc5_BUCB_squared.Rda"))
loc5_comp_BUCB_squared <- compare_to_ground_truth(mean_loc5, loc5_BUCB_squared, 
                                                  tau_loc5, epsilon_loc5)$mean
save(loc5_comp_BUCB_squared,
     file = paste0(current_path, "loc5_comp_BUCB_squared.Rda"))
load(file = paste0(current_path, "loc5_comp_BUCB_squared.Rda"))
rm(loc5_BUCB_squared)
gc()

########################################################################

loc5_TTS <- para_bandit_sim_TTS(data = data_list5, rounds = 2500,
                                tau = tau_loc5, epsilon = epsilon_loc5,
                                alpha = tau_loc5, beta = 1 - tau_loc5)
save(loc5_TTS, file = paste0(current_path, "loc5_TTS.Rda"))
#load(file = paste0(current_path, "loc5_TTS.Rda"))
loc5_comp_TTS <- compare_to_ground_truth(mean_loc5, loc5_TTS, 
                                         tau_loc5, epsilon_loc5)$mean
save(loc5_comp_TTS, file = paste0(current_path, "loc5_comp_TTS.Rda"))
rm(loc5_TTS)
gc()

########################################################################
# Standard AugUCB
system.time(loc5_AugUCB <- para_bandit_sim_AugUCB(data = data_list5, rounds = 2500, 
                                                  tau = tau_loc5))
save(loc5_AugUCB, file = paste0(current_path, "loc5_AugUCB.Rda"))
loc5_comp_AugUCB <- compare_to_ground_truth(mean_loc5, loc5_AugUCB, tau_loc5, 
                                            epsilon_loc5)$mean
save(loc5_comp_AugUCB, file = paste0(current_path, "loc5_comp_AugUCB.Rda"))

########################################################################

# Standard Uniform
system.time(loc5_UNIFORM <- para_bandit_sim_uniform(data = data_list5, 
                                                    rounds = 2500))
save(loc5_UNIFORM, file = paste0(current_path, "loc5_UNIFORM.Rda"))
#load(file = paste0(current_path, "loc5_UNIFORM.Rda"))
loc5_comp_UNIFORM <- compare_to_ground_truth(mean_loc5, loc5_UNIFORM, tau_loc5, 
                                             epsilon_loc5)$mean
save(loc5_comp_UNIFORM, file = paste0(current_path, "loc5_comp_UNIFORM.Rda"))
rm(loc5_UNIFORM)
gc()

########################################################################
# Do KL by comparing against tau directly

system.time(loc5_KL <- para_bandit_sim_KL(data = data_list5, rounds = 1000, 
                                          tau = tau_loc5, epsilon = epsilon_loc5, 
                                          at_tau = TRUE, horizon = 1000))
# user  system elapsed 
#6.824   4.335 391.240
save(loc5_KL, file = paste0(current_path, "loc5_KL.Rda"))
loc5_comp_KL <- compare_to_ground_truth(mean_loc5, loc5_KL, tau_loc5, 
                                        epsilon_loc5)$mean
save(loc5_comp_KL, file = paste0(current_path, "loc5_comp_KL.Rda"))

########################################################################
# Do KL by comparing against tau and epsilon

system.time(loc5_KL_not_tau <- para_bandit_sim_KL(data = data_list5, rounds = 1000, 
                                                  tau = tau_loc5, epsilon = epsilon_loc5, 
                                                  at_tau = FALSE))
save(loc5_KL_not_tau, file = paste0(current_path, "loc5_KL_not_tau.Rda"))
loc5_comp_KL_not_tau <- compare_to_ground_truth(mean_loc5, loc5_KL_not_tau, 
                                                tau_loc5, epsilon_loc5)$mean
save(loc5_comp_KL_not_tau, file = paste0(current_path,
                                         "loc5_comp_KL_not_tau.Rda"))

########################################################################
# Do KL by comparing against tau and epsilon with 0/1 adjustment

system.time(loc5_KL_horizon <- para_bandit_sim_KL(data = data_list5, rounds = 1000, 
                                                  tau = tau_loc5, epsilon = epsilon_loc5, 
                                                  at_tau = FALSE, horizon = 10000))
save(loc5_KL_horizon, file = paste0(current_path, "loc5_KL_horizon.Rda"))
loc5_comp_KL_horizon <- compare_to_ground_truth(mean_loc5, loc5_KL_horizon, 
                                                tau_loc5, epsilon_loc5)$mean
save(loc5_comp_KL_horizon, file = paste0(current_path,
                                         "loc5_comp_KL_horizon.Rda"))


########################################################################
# Do KL by comparing against tau and epsilon with 0/1 adjustment

system.time(loc5_KL_horizon2 <- para_bandit_sim_KL(data = data_list5, rounds = 1000, 
                                                   tau = tau_loc5, epsilon = epsilon_loc5, 
                                                   at_tau = FALSE, horizon = 1000^2))
save(loc5_KL_horizon2, file = paste0(current_path, "loc5_KL_horizon2.Rda"))
loc5_comp_KL_horizon2 <- compare_to_ground_truth(mean_loc5, loc5_KL_horizon2, 
                                                 tau_loc5, epsilon_loc5)$mean
save(loc5_comp_KL_horizon2, file = paste0(current_path,
                                          "loc5_comp_KL_horizon2.Rda"))

########################################################################
# Do KL by comparing against tau and epsilon with 0/1 adjustment

system.time(loc5_KL_horizon1 <- para_bandit_sim_KL(data = data_list5, rounds = 1000, 
                                                   tau = tau_loc5, epsilon = epsilon_loc5, 
                                                   at_tau = FALSE, horizon = 1000))
save(loc5_KL_horizon1, file = paste0(current_path, "loc5_KL_horizon1.Rda"))
load(file = paste0(current_path, "loc5_KL_horizon1.Rda"))
loc5_comp_KL_horizon1 <- compare_to_ground_truth(mean_loc5, loc5_KL_horizon1, 
                                                 tau_loc5, epsilon_loc5)$mean
save(loc5_comp_KL_horizon1, file = paste0(current_path,
                                          "loc5_comp_KL_horizon1.Rda"))
rm(loc5_KL_horizon1)
gc()
########################################################################
# Probability that arm is above or below tau+-epsilon

loc5_PI <- para_bandit_sim_PI(data = data_list5, rounds = 1000,
                              tau = tau_loc5, epsilon = epsilon_loc5,
                              alpha = tau_loc5, beta = 1 - tau_loc5)
save(loc5_PI, file = paste0(current_path, "loc5_PI.Rda"))
loc5_comp_PI <- compare_to_ground_truth(mean_loc5, loc5_PI, 
                                        tau_loc5, epsilon_loc5)$mean
save(loc5_comp_PI, file = paste0(current_path, "loc5_comp_PI.Rda"))

########################################################################
# Probability that arm is above or below tau+-epsilon

loc5_PI_tadj <- para_bandit_sim_PI(data = data_list5, rounds = 1000,
                                   tau = tau_loc5, epsilon = epsilon_loc5,
                                   alpha = tau_loc5, beta = 1 - tau_loc5, tadj = TRUE)
save(loc5_PI_tadj, file = paste0(current_path, "loc5_PI_tadj.Rda"))
loc5_comp_PI_tadj <- compare_to_ground_truth(mean_loc5, loc5_PI_tadj, 
                                             tau_loc5, epsilon_loc5)$mean
save(loc5_comp_PI_tadj, file = paste0(current_path, "loc5_comp_PI_tadj.Rda"))

########################################################################

loc5_TTS <- para_bandit_sim_TTS(data = data_list5, rounds = 1000,
                                tau = tau_loc5, epsilon = epsilon_loc5,
                                alpha = tau_loc5, beta = 1 - tau_loc5)
save(loc5_TTS, file = paste0(current_path, "loc5_TTS.Rda"))
load(file = paste0(current_path, "loc5_TTS.Rda"))
loc5_comp_TTS <- compare_to_ground_truth(mean_loc5, loc5_TTS, 
                                         tau_loc5, epsilon_loc5)$mean
save(loc5_comp_TTS, file = paste0(current_path, "loc5_comp_TTS.Rda"))
rm(loc5_TTS)
gc()
########################################################################

#system.time(loc5_BETA <- para_bandit_sim_BETA(data = data_list5, rounds = 1000,
#                                              tau = tau_loc5, epsilon = epsilon_loc5,
#                                              alpha = tau_loc5, beta = 1 - tau_loc5))
#save(loc5_BETA, file = paste0(current_path, "loc5_BETA.Rda"))
#loc5_comp_BETA <- compare_to_ground_truth(mean_loc5, loc5_BETA, 
#                                          tau_loc5, epsilon_loc5)$mean
#save(loc5_comp_BETA, file = paste0(current_path, "loc5_comp_BETA.Rda"))

########################################################################

loc5_BUCB_squared_long <- para_bandit_sim_bucb(data = data_list5_long, rounds = 2000, 
                                               rate = "inverse_squared",
                                               tau = tau_loc5, epsilon = epsilon_loc5, 
                                               alpha = tau_loc5, beta = 1-tau_loc5)
save(loc5_BUCB_squared_long, file = paste0(current_path, "loc5_BUCB_squared_long.Rda"))
load(file = paste0(current_path, "loc5_BUCB_squared_long.Rda"))
loc5_comp_BUCB_squared_long <- compare_to_ground_truth(mean_loc5, loc5_BUCB_squared_long, 
                                                       tau_loc5, epsilon_loc5)$mean
save(loc5_comp_BUCB_squared_long, file = paste0(current_path, "loc5_comp_BUCB_squared_long.Rda"))
rm(loc5_BUCB_squared_long)
gc()
########################################################################

loc5_BUCB_cubic <- para_bandit_sim_bucb(data = data_list5, rounds = 1000, 
                                        rate = "inverse_cubic",
                                        tau = tau_loc5, epsilon = epsilon_loc5, 
                                        alpha = tau_loc5, beta = 1-tau_loc5)
save(loc5_BUCB_cubic, file = paste0(current_path, "loc5_BUCB_cubic.Rda"))
loc5_comp_BUCB_cubic <- compare_to_ground_truth(mean_loc5, loc5_BUCB_cubic, 
                                                tau_loc5, epsilon_loc5)$mean
save(loc5_comp_BUCB_cubic, file = paste0(current_path, "loc5_comp_BUCB_cubic.Rda"))

########################################################################

loc5_BUCB_power5 <- para_bandit_sim_bucb(data = data_list5, rounds = 1000, 
                                         rate = "inverse_power5",
                                         tau = tau_loc5, epsilon = epsilon_loc5, 
                                         alpha = tau_loc5, beta = 1-tau_loc5)
save(loc5_BUCB_power5, file = paste0(current_path, "loc5_BUCB_power5.Rda"))
load(file = paste0(current_path, "loc5_BUCB_power5.Rda"))
loc5_comp_BUCB_power5 <- compare_to_ground_truth(mean_loc5, loc5_BUCB_power5, 
                                                 tau_loc5, epsilon_loc5)$mean
save(loc5_comp_BUCB_power5, file = paste0(current_path, "loc5_comp_BUCB_power5.Rda"))
rm(loc5_BUCB_power5)
gc()
########################################################################

loc5_BUCB_squared_no_eps <- para_bandit_sim_bucb(data = data_list5, rounds = 1000, 
                                                 rate = "inverse_squared",
                                                 tau = tau_loc5, epsilon = epsilon_loc5, 
                                                 alpha = tau_loc5, beta = 1-tau_loc5,
                                                 with_epsilon = FALSE)
save(loc5_BUCB_squared_no_eps, 
     file = paste0(current_path, "loc5_BUCB_squared_no_eps.Rda"))
loc5_comp_BUCB_squared_no_eps <- compare_to_ground_truth(mean_loc5, 
                                                         loc5_BUCB_squared_no_eps,
                                                         tau_loc5,
                                                         epsilon_loc5)$mean
save(loc5_comp_BUCB_squared_no_eps, file = paste0(current_path, "loc5_comp_BUCB_squared_no_eps.Rda"))

########################################################################

loc5_BUCB_horizon <- para_bandit_sim_bucb(data = data_list5, rounds = 1000, 
                                          rate = "inverse_horizon",
                                          tau = tau_loc5, epsilon = epsilon_loc5, 
                                          alpha = tau_loc5, beta = 1-tau_loc5,
                                          with_epsilon = FALSE)
save(loc5_BUCB_horizon, 
     file = paste0(current_path, "loc5_BUCB_horizon.Rda"))
loc5_comp_BUCB_horizon <- compare_to_ground_truth(mean_loc5, 
                                                  loc5_BUCB_horizon,
                                                  tau_loc5,
                                                  epsilon_loc5)$mean
save(loc5_comp_BUCB_horizon, file = paste0(current_path, "loc5_comp_BUCB_horizon.Rda"))

########################################################################

loc5_BUCB_horizon5 <- para_bandit_sim_bucb(data = data_list5, rounds = 1000, 
                                           rate = "inverse_horizon_c",
                                           tau = tau_loc5, epsilon = epsilon_loc5, 
                                           alpha = tau_loc5, beta = 1-tau_loc5,
                                           with_epsilon = FALSE, const = 5)
save(loc5_BUCB_horizon5, 
     file = paste0(current_path, "loc5_BUCB_horizon5.Rda"))
loc5_comp_BUCB_horizon5 <- compare_to_ground_truth(mean_loc5, 
                                                   loc5_BUCB_horizon5,
                                                   tau_loc5,
                                                   epsilon_loc5)$mean
save(loc5_comp_BUCB_horizon5, file = paste0(current_path, "loc5_comp_BUCB_horizon5.Rda"))

########################################################################

load(paste0(current_path, "loc5_comp_APT.Rda"))
load(paste0(current_path, "loc5_comp_AugUCB.Rda"))
load(paste0(current_path, "loc5_comp_TTS.Rda"))
load(paste0(current_path, "loc5_comp_BUCB.Rda"))
load(paste0(current_path, "loc5_comp_BUCB_horizon.Rda"))
load(paste0(current_path, "loc5_comp_BUCB_squared.Rda"))

plot(c(0,3000), c(0,-2), type = "n")
lines(log(loc5_comp_APT), col = "blue")
lines(log(loc5_comp_TTS), col = "lightblue")
lines(log(loc5_comp_BUCB), col = "green")
lines(log(loc5_comp_BUCB_horizon), col = "pink")
lines(log(loc5_comp_BUCB_squared), col = "red")
lines(log(loc5_comp_AugUCB), col = "grey")

########################################################################

library(plyr)
loc5_UNIFORM_as <- colMeans(ldply(loc5_UNIFORM, function(x) table(x$arm_sequence)))
loc5_APT_as <- colMeans(ldply(loc5_APT, function(x) table(x$arm_sequence)))
loc5_PI_as <- colMeans(ldply(loc5_PI, function(x) table(x$arm_sequence)))
loc5_PI_tadj_as <- colMeans(ldply(loc5_PI_tadj, function(x) table(x$arm_sequence)))
loc5_TTS_as <- colMeans(ldply(loc5_TTS, function(x) table(x$arm_sequence)))
loc5_BUCB_as <- colMeans(ldply(loc5_BUCB, function(x) table(x$arm_sequence)))
loc5_BUCB_squared_as <- colMeans(ldply(loc5_BUCB_squared, function(x) table(x$arm_sequence)))
loc5_BUCB_cubic_as <- colMeans(ldply(loc5_BUCB_cubic, function(x) table(x$arm_sequence)))
loc5_KL_as <- colMeans(ldply(loc5_KL, function(x) table(x$arm_sequence)))
loc5_KL_not_tau_as <- colMeans(ldply(loc5_KL_not_tau, function(x) table(x$arm_sequence)))
loc5_KL_horizon1_as <- colMeans(ldply(loc5_KL_horizon1, function(x) table(x$arm_sequence)))

round(data.frame(loc5_UNIFORM_as, loc5_APT_as, 
                 loc5_BUCB_as, loc5_BUCB_squared_as,
                 loc5_TTS_as#, loc5_KL_as, loc5_KL_horizon1_as
))


loc5_UNIFORM_as_1 <- colMeans(ldply(loc5_UNIFORM, function(x) table(x$arm_sequence[1:500])))
loc5_APT_as_1 <- colMeans(ldply(loc5_APT, function(x) table(x$arm_sequence[1:500])))
loc5_PI_as_1 <- colMeans(ldply(loc5_PI, function(x) table(x$arm_sequence[1:500])))
loc5_PI_tadj_as_1 <- colMeans(ldply(loc5_PI_tadj, function(x) table(x$arm_sequence[1:500])))
loc5_TTS_as_1 <- colMeans(ldply(loc5_TTS, function(x) table(x$arm_sequence[1:500])))
loc5_BUCB_as_1 <- colMeans(ldply(loc5_BUCB, function(x) table(x$arm_sequence[1:500])))
loc5_KL_as_1 <- colMeans(ldply(loc5_KL, function(x) table(x$arm_sequence[1:500])))
loc5_KL_not_tau_as_1 <- colMeans(ldply(loc5_KL_not_tau, function(x) table(x$arm_sequence[1:500])))
loc5_KL_horizon1_as_1 <- colMeans(ldply(loc5_KL_horizon1, function(x) table(x$arm_sequence[1:500])))

round(data.frame(loc5_UNIFORM_as_1, loc5_APT_as_1, loc5_BUCB_as_1,
                 loc5_PI_as_1, loc5_PI_tadj_as_1, loc5_TTS_as_1,
                 loc5_KL_as_1, loc5_KL_horizon1_as_1))/500


loc5_UNIFORM_as_2 <- colMeans(ldply(loc5_UNIFORM, function(x) table(x$arm_sequence[501:1000])))
loc5_APT_as_2 <- colMeans(ldply(loc5_APT, function(x) table(x$arm_sequence[501:1000])))
loc5_PI_as_2 <- colMeans(ldply(loc5_PI, function(x) table(x$arm_sequence[501:1000])))
loc5_PI_tadj_as_2 <- colMeans(ldply(loc5_PI_tadj, function(x) table(x$arm_sequence[501:1000])))
loc5_TTS_as_2 <- colMeans(ldply(loc5_TTS, function(x) table(x$arm_sequence[501:1000])))
loc5_BUCB_as_2 <- colMeans(ldply(loc5_BUCB, function(x) table(x$arm_sequence[501:1000])))
loc5_KL_as_2 <- colMeans(ldply(loc5_KL, function(x) table(x$arm_sequence[501:1000])))
loc5_KL_not_tau_as_2 <- colMeans(ldply(loc5_KL_not_tau, function(x) table(x$arm_sequence[501:1000])))
loc5_KL_horizon1_as_2 <- colMeans(ldply(loc5_KL_horizon1, function(x) table(x$arm_sequence[501:1000])))

round(data.frame(loc5_UNIFORM_as, loc5_APT_as, loc5_BUCB_as,
                 loc5_PI_as, loc5_PI_tadj_as, loc5_TTS_as,
                 loc5_KL_as, loc5_KL_horizon1_as))

########################################################################

plot(loc5_BUCB_squared[[1]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc5_BUCB_squared[[2]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc5_BUCB_squared[[3]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc5_BUCB_squared[[4]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc5_BUCB_squared[[5]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc5_BUCB_squared[[6]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc5_BUCB_squared[[7]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc5_BUCB_squared[[8]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc5_BUCB_squared[[9]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc5_BUCB_squared[[10]]$arm_sequence, pch = 19, cex = 0.3)

plot(loc5_APT[[1]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc5_APT[[2]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc5_APT[[3]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc5_APT[[4]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc5_APT[[5]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc5_APT[[6]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc5_APT[[7]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc5_APT[[8]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc5_APT[[9]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc5_APT[[10]]$arm_sequence, pch = 19, cex = 0.3)

arm_seq_res_APT <- data.frame(t(laply(loc5_APT, function(x) x$arm_sequence)))
library(tidyr)
library(dplyr)
#arm_seq_table <- arm_seq_res %>% tbl_df() %>% mutate(index = 1:1000) %>%
#  gather(key = iter, value = arm, -index) %>%
#  group_by(index, arm) %>% summarize(count = n())

library(ggplot2)
arm_seq_res_APT %>% tbl_df() %>% mutate(index = 1:2000) %>%
  gather(key = iter, value = arm, -index) %>%
  filter(index > 11) %>%
  ggplot(aes(index, arm)) + geom_count() +
  scale_size_area(max_size = 10)
#ggplot(arm_seq_table[11:9910,], aes(index, arm)) + geom_tile(aes(fill = count))

arm_seq_res_BUCB_squared <- data.frame(t(laply(loc5_BUCB_squared, function(x) x$arm_sequence)))
arm_seq_res_BUCB_squared %>% tbl_df() %>% mutate(index = 1:2000) %>%
  gather(key = iter, value = arm, -index) %>%
  filter(index > 11) %>%
  ggplot(aes(index, arm)) + geom_count() +
  scale_size_area(max_size = 10)

arm_seq_res_TTS <- data.frame(t(laply(loc5_TTS, function(x) x$arm_sequence)))
arm_seq_res_TTS %>% tbl_df() %>% mutate(index = 1:2000) %>%
  gather(key = iter, value = arm, -index) %>%
  filter(index > 11) %>%
  ggplot(aes(index, arm)) + geom_count() +
  scale_size_area(max_size = 10)
