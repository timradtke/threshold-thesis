########################################################################
# load functions

current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/AugUCB/"
source(paste0(current_path, "run_from_data.R"))
source(paste0(current_path, "run_parallelized_from_data.R"))
current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/synthetic_experiments/"

error_bound <- function(n, H, K) {
  (K*n+H)*exp(-n/H)
}

plot(1:360000, log(error_bound(1:360000, H_loc4, 11)), type = "l",
     ylim = c(-6,12), main = round(H_loc4))
abline(h=log(0.2))

########################################################################
# Create the data

# do 2000 rounds
mean_loc4 <- c(0.0005, 0.001, 0.001, 0.005, 0.009,
               0.02, 0.0275, 0.0325, 0.04,
               0.08, 0.09)
tau_loc4 <- 0.03
epsilon_loc4 <- 0.005
H_loc4 <- get_complexity(mean_loc4, tau_loc4, epsilon_loc4)

plot(mean_loc4, rep(1,11), main = paste0("Complexity of ", round(H_loc4,2)))
abline(v=tau_loc4)
abline(v=tau_loc4+epsilon_loc4, lty=2)
abline(v=tau_loc4-epsilon_loc4, lty=2)

data_list4 <- list()
set.seed(1024)
for(j in 1:5000) {
  curr_data <- data.frame(rep(NA, times = 6000))
  for(i in 1:length(mean_loc4)) {
    curr_data[[i]] <- as.numeric(purrr::rbernoulli(6000, p  = mean_loc4[i]))
  }
  names(curr_data) <- paste0("V", rep(1:length(mean_loc4)))
  data_list4[[j]] <- curr_data
}

data_list4_8000 <- list()
set.seed(1024)
for(j in 1:5000) {
  curr_data <- data.frame(rep(NA, times = 8000))
  for(i in 1:length(mean_loc4)) {
    curr_data[[i]] <- as.numeric(purrr::rbernoulli(8000, p  = mean_loc4[i]))
  }
  names(curr_data) <- paste0("V", rep(1:length(mean_loc4)))
  data_list4_8000[[j]] <- curr_data
}

########################################################################

system.time(loc4_KLUCB90_long <- para_bandit_sim_KLUCB(data = data_list4_8000, 
                                                    rounds = 8000, 
                                                    tau = tau_loc4, 
                                                    epsilon = epsilon_loc4,
                                                    horizon = 90000,
                                                    H = H_loc4))
save(loc4_KLUCB90_long, file = paste0(current_path, "loc4_KLUCB90_long.Rda"))
loc4_comp_KLUCB90_long <- compare_to_ground_truth(mean_loc4, loc4_KLUCB90_long, 
                                               tau_loc4, 
                                               epsilon_loc4)$mean
save(loc4_comp_KLUCB90_long, file = paste0(current_path, 
                                        "loc4_comp_KLUCB90_long.Rda"))

########################################################################
data_list4_8000_small <- data_list4_8000[1:500]
rm(data_list4)
gc()

system.time(loc4_KLUCB24100 <- para_bandit_sim_KLUCB(data = data_list4_8000_small, 
                                                     rounds = 2000, 
                                                     tau = tau_loc4, 
                                                     epsilon = epsilon_loc4,
                                                     horizon = 241000,
                                                     H = H_loc4))
save(loc4_KLUCB24100, file = paste0(current_path, "loc4_KLUCB24100.Rda"))
loc4_comp_KLUCB24100 <- compare_to_ground_truth(mean_loc4, loc4_KLUCB24100, 
                                                tau_loc4, 
                                                epsilon_loc4)$mean
save(loc4_comp_KLUCB24100, file = paste0(current_path, 
                                         "loc4_comp_KLUCB24100.Rda"))

system.time(loc4_KLUCB12 <- para_bandit_sim_KLUCB(data = data_list4_8000_small, 
                                                     rounds = 2000, 
                                                     tau = tau_loc4, 
                                                     epsilon = epsilon_loc4,
                                                     horizon = 12000,
                                                     H = H_loc4))
save(loc4_KLUCB12, file = paste0(current_path, "loc4_KLUCB12.Rda"))
loc4_comp_KLUCB12 <- compare_to_ground_truth(mean_loc4, loc4_KLUCB12, 
                                                tau_loc4, 
                                                epsilon_loc4)$mean
save(loc4_comp_KLUCB12, file = paste0(current_path, 
                                         "loc4_comp_KLUCB12.Rda"))

system.time(loc4_KLUCB70 <- para_bandit_sim_KLUCB(data = data_list4_8000_small, 
                                                  rounds = 2000, 
                                                  tau = tau_loc4, 
                                                  epsilon = epsilon_loc4,
                                                  horizon = 70000,
                                                  H = H_loc4))
save(loc4_KLUCB70, file = paste0(current_path, "loc4_KLUCB70.Rda"))
loc4_comp_KLUCB70 <- compare_to_ground_truth(mean_loc4, loc4_KLUCB70, 
                                             tau_loc4, 
                                             epsilon_loc4)$mean
save(loc4_comp_KLUCB70, file = paste0(current_path, 
                                      "loc4_comp_KLUCB70.Rda"))

system.time(loc4_KLUCB90_2 <- para_bandit_sim_KLUCB(data = data_list4_8000_small, 
                                                   rounds = 2000, 
                                                   tau = tau_loc4, 
                                                   epsilon = epsilon_loc4,
                                                   horizon = 90000,
                                                   H = H_loc4))
save(loc4_KLUCB90_2, file = paste0(current_path, "loc4_KLUCB90_2.Rda"))
loc4_comp_KLUCB90_2 <- compare_to_ground_truth(mean_loc4, loc4_KLUCB90_2, 
                                              tau_loc4, 
                                              epsilon_loc4)$mean
save(loc4_comp_KLUCB90_2, file = paste0(current_path, 
                                       "loc4_comp_KLUCB90_2.Rda"))

system.time(loc4_KLUCB110 <- para_bandit_sim_KLUCB(data = data_list4_8000_small, 
                                                  rounds = 2000, 
                                                  tau = tau_loc4, 
                                                  epsilon = epsilon_loc4,
                                                  horizon = 110000,
                                                  H = H_loc4))
save(loc4_KLUCB110, file = paste0(current_path, "loc4_KLUCB110.Rda"))
loc4_comp_KLUCB110 <- compare_to_ground_truth(mean_loc4, loc4_KLUCB110, 
                                             tau_loc4, 
                                             epsilon_loc4)$mean
save(loc4_comp_KLUCB110, file = paste0(current_path, 
                                      "loc4_comp_KLUCB110.Rda"))



########################################################################
# Standard Uniform
system.time(loc4_UNIFORM <- para_bandit_sim_uniform(data = data_list4, 
                                                    rounds = 6000))
# user  system elapsed 
# 6.268   3.400 211.002 
save(loc4_UNIFORM, file = paste0(current_path, "loc4_UNIFORM.Rda"))
load(file = paste0(current_path, "loc4_UNIFORM.Rda"))
loc4_comp_UNIFORM <- compare_to_ground_truth(mean_loc4, loc4_UNIFORM, tau_loc4, 
                                             epsilon_loc4)$mean
save(loc4_comp_UNIFORM, file = paste0(current_path, "loc4_comp_UNIFORM.Rda"))
rm(loc4_UNIFORM)
gc()
########################################################################
# Standard APT Algorithm
# a seed is built in
system.time(loc4_APT <- para_bandit_sim_APT(data = data_list4, rounds = 6000, 
                                            tau = tau_loc4, epsilon = epsilon_loc4))
# user  system elapsed 
#6.508   3.544 300.383
save(loc4_APT, file = paste0(current_path, "loc4_APT.Rda"))
load(file = paste0(current_path, "loc4_APT.Rda"))
loc4_comp_APT <- compare_to_ground_truth(mean_loc4, loc4_APT, tau_loc4, 
                                         epsilon_loc4)$mean
save(loc4_comp_APT, file = paste0(current_path, "loc4_comp_APT.Rda"))
load(file = paste0(current_path, "loc4_comp_APT.Rda"))
rm(loc4_APT)
gc()

########################################################################
# KL-UCB Algorithm with tuned exploration rate through horizon=40000
system.time(loc4_KLUCB180 <- para_bandit_sim_KLUCB(data = data_list4[1:100], 
                                                  rounds = 1000, 
                                                  tau = tau_loc4, 
                                                  epsilon = epsilon_loc4,
                                                  horizon = 180000,
                                                  H = H_loc4))
# user  system elapsed 
#1.019   1.369 553.760
save(loc4_KLUCB180, file = paste0(current_path, "loc4_KLUCB180.Rda"))
#load(file = paste0(current_path, "loc4_KLUCB40.Rda"))
loc4_comp_KLUCB180 <- compare_to_ground_truth(mean_loc4, loc4_KLUCB180, tau_loc4, 
                                         epsilon_loc4)$mean
save(loc4_comp_KLUCB180, file = paste0(current_path, "loc4_comp_KLUCB180.Rda"))
#load(file = paste0(current_path, "loc4_comp_KLUCB40.Rda"))
rm(loc4_KLUCB180)
gc()

data_list4_small <- data_list4[1:1000]
rm(data_list4)
gc()
system.time(loc4_KLUCB90 <- para_bandit_sim_KLUCB(data = data_list4_small, 
                                                   rounds = 6000, 
                                                   tau = tau_loc4, 
                                                   epsilon = epsilon_loc4,
                                                   horizon = 90000,
                                                   H = H_loc4))
# user  system elapsed 
#51.048    73.114 34050.023
save(loc4_KLUCB90, file = paste0(current_path, "loc4_KLUCB90.Rda"))
#load(file = paste0(current_path, "loc4_KLUCB40.Rda"))
loc4_comp_KLUCB90 <- compare_to_ground_truth(mean_loc4, loc4_KLUCB90, tau_loc4, 
                                              epsilon_loc4)$mean
save(loc4_comp_KLUCB90, file = paste0(current_path, "loc4_comp_KLUCB90.Rda"))
#load(file = paste0(current_path, "loc4_comp_KLUCB40.Rda"))
rm(loc4_KLUCB90)
gc()

system.time(loc4_KLUCB90 <- para_bandit_sim_KLUCB(data = data_list4_small, 
                                                  rounds = 2000, 
                                                  tau = tau_loc4, 
                                                  epsilon = epsilon_loc4,
                                                  horizon = 90000,
                                                  H = H_loc4))
# user  system elapsed 
#51.048    73.114 34050.023
save(loc4_KLUCB90, file = paste0(current_path, "loc4_KLUCB90.Rda"))
#load(file = paste0(current_path, "loc4_KLUCB40.Rda"))
loc4_comp_KLUCB90 <- compare_to_ground_truth(mean_loc4, loc4_KLUCB90, tau_loc4, 
                                             epsilon_loc4)$mean
save(loc4_comp_KLUCB90, file = paste0(current_path, "loc4_comp_KLUCB90.Rda"))
#load(file = paste0(current_path, "loc4_comp_KLUCB40.Rda"))
rm(loc4_KLUCB90)
gc()

########################################################################

loc4_BUCB_horizon_long <- para_bandit_sim_bucb(data = data_list4, rounds = 6000, 
                                          rate = "inverse_horizon",
                                          tau = tau_loc4, epsilon = epsilon_loc4, 
                                          alpha = tau_loc4, beta = 1-tau_loc4)
save(loc4_BUCB_horizon_long, 
     file = paste0(current_path, "loc4_BUCB_horizon_long.Rda"))
loc4_comp_BUCB_horizon_long <- compare_to_ground_truth(mean_loc4, 
                                                  loc4_BUCB_horizon_long,
                                                  tau_loc4,
                                                  epsilon_loc4)$mean
save(loc4_comp_BUCB_horizon_long,
     file = paste0(current_path, "loc4_comp_BUCB_horizon_long.Rda"))
rm(loc4_BUCB_horizon_long)
gc()


# now with 8000

loc4_BUCB_horizon_8000 <- para_bandit_sim_bucb(data = data_list4_8000, 
                                               rounds = 8000, 
                                               rate = "inverse_horizon",
                                               tau = tau_loc4, 
                                               epsilon = epsilon_loc4, 
                                               alpha = tau_loc4, 
                                               beta = 1-tau_loc4)
save(loc4_BUCB_horizon_8000, 
     file = paste0(current_path, "loc4_BUCB_horizon_8000.Rda"))
loc4_comp_BUCB_horizon_8000 <- compare_to_ground_truth(mean_loc4, 
                                                       loc4_BUCB_horizon_8000,
                                                       tau_loc4,
                                                       epsilon_loc4)$mean
save(loc4_comp_BUCB_horizon_8000,
     file = paste0(current_path, "loc4_comp_BUCB_horizon_8000.Rda"))
rm(loc4_BUCB_horizon_8000)
gc()

########################################################################

loc4_BUCB_horizon_uninfprior_long <- para_bandit_sim_bucb(data = data_list4, rounds = 6000, 
                                               rate = "inverse_horizon",
                                               tau = tau_loc4, epsilon = epsilon_loc4, 
                                               alpha = 1, beta = 1)
save(loc4_BUCB_horizon_uninfprior_long, 
     file = paste0(current_path, "loc4_BUCB_horizon_uninfprior_long.Rda"))
loc4_comp_BUCB_horizon_uninfprior_long <- compare_to_ground_truth(mean_loc4, 
                                                       loc4_BUCB_horizon_uninfprior_long,
                                                       tau_loc4,
                                                       epsilon_loc4)$mean
save(loc4_comp_BUCB_horizon_uninfprior_long,
     file = paste0(current_path, "loc4_comp_BUCB_horizon_uninfprior_long.Rda"))
rm(loc4_BUCB_horizon_uninfprior_long)
gc()

########################################################################

loc4_BUCB_linearhorizon_long <- para_bandit_sim_bucb(data = data_list4, rounds = 6000, 
                                               rate = "inverse_horizon_linear",
                                               tau = tau_loc4, epsilon = epsilon_loc4, 
                                               alpha = tau_loc4, beta = 1-tau_loc4)
save(loc4_BUCB_linearhorizon_long, 
     file = paste0(current_path, "loc4_BUCB_linearhorizon_long.Rda"))
loc4_comp_BUCB_linearhorizon_long <- compare_to_ground_truth(mean_loc4, 
                                                       loc4_BUCB_linearhorizon_long,
                                                       tau_loc4,
                                                       epsilon_loc4)$mean
save(loc4_comp_BUCB_linearhorizon_long,
     file = paste0(current_path, "loc4_comp_BUCB_linearhorizon_long.Rda"))
rm(loc4_BUCB_linearhorizon_long)
gc()

########################################################################

H <- sum(1/(abs(mean_loc4-tau_loc4)+epsilon_loc4)^2)
loc4_BUCB_linearhorizonC_long <- para_bandit_sim_bucb(data = data_list4, rounds = 6000, 
                                                     rate = "inverse_horizon_linear_c",
                                                     tau = tau_loc4, epsilon = epsilon_loc4, 
                                                     alpha = tau_loc4, beta = 1-tau_loc4,
                                                     const = H)
save(loc4_BUCB_linearhorizonC_long, 
     file = paste0(current_path, "loc4_BUCB_linearhorizonC_long.Rda"))
loc4_comp_BUCB_linearhorizonC_long <- compare_to_ground_truth(mean_loc4, 
                                                             loc4_BUCB_linearhorizonC_long,
                                                             tau_loc4,
                                                             epsilon_loc4)$mean
save(loc4_comp_BUCB_linearhorizonC_long,
     file = paste0(current_path, "loc4_comp_BUCB_linearhorizonC_long.Rda"))
rm(loc4_BUCB_linearhorizonC_long)
gc()

########################################################################

H <- sum(1/(abs(mean_loc4-tau_loc4)+epsilon_loc4)^2)
loc4_BUCB_linearhorizonC_uninf_long <- para_bandit_sim_bucb(data = data_list4, rounds = 6000, 
                                                      rate = "inverse_horizon_linear_c",
                                                      tau = tau_loc4, epsilon = epsilon_loc4, 
                                                      alpha = 1, beta = 1,
                                                      const = H)
save(loc4_BUCB_linearhorizonC_uninf_long, 
     file = paste0(current_path, "loc4_BUCB_linearhorizonC_uninf_long.Rda"))
loc4_comp_BUCB_linearhorizonC_uninf_long <- compare_to_ground_truth(mean_loc4, 
                                                              loc4_BUCB_linearhorizonC_uninf_long,
                                                              tau_loc4,
                                                              epsilon_loc4)$mean
save(loc4_comp_BUCB_linearhorizonC_uninf_long,
     file = paste0(current_path, "loc4_comp_BUCB_linearhorizonC_uninf_long.Rda"))
rm(loc4_BUCB_linearhorizonC_uninf_long)
gc()

########################################################################

loc4_BUCB_horizon5 <- para_bandit_sim_bucb(data = data_list4, rounds = 2000, 
                                           rate = "inverse_horizon_c",
                                           tau = tau_loc4, epsilon = epsilon_loc4, 
                                           alpha = tau_loc4, beta = 1-tau_loc4,
                                           const = 5)
save(loc4_BUCB_horizon5, 
     file = paste0(current_path, "loc4_BUCB_horizon5.Rda"))
loc4_comp_BUCB_horizon5 <- compare_to_ground_truth(mean_loc4, 
                                                   loc4_BUCB_horizon5,
                                                   tau_loc4,
                                                   epsilon_loc4)$mean
save(loc4_comp_BUCB_horizon5, file = paste0(current_path, "loc4_comp_BUCB_horizon5.Rda"))
rm(loc4_BUCB_horizon5)
gc()

########################################################################

loc4_BUCB <- para_bandit_sim_bucb(data = data_list4, rounds = 2000, 
                                  rate = "inverse",
                                  tau = tau_loc4, epsilon = epsilon_loc4, 
                                  alpha = tau_loc4, beta = 1-tau_loc4)
save(loc4_BUCB, file = paste0(current_path, "loc4_BUCB.Rda"))
load(file = paste0(current_path, "loc4_BUCB.Rda"))
loc4_comp_BUCB <- compare_to_ground_truth(mean_loc4, loc4_BUCB, 
                                          tau_loc4, epsilon_loc4)$mean
save(loc4_comp_BUCB, file = paste0(current_path, "loc4_comp_BUCB.Rda"))
load(file = paste0(current_path, "loc4_comp_BUCB.Rda"))
rm(loc4_BUCB)
gc()

########################################################################

loc4_BUCB_squared_long <- para_bandit_sim_bucb(data = data_list4, rounds = 6000, 
                                          rate = "inverse_squared",
                                          tau = tau_loc4, epsilon = epsilon_loc4, 
                                          alpha = tau_loc4, beta = 1-tau_loc4)
save(loc4_BUCB_squared_long, file = paste0(current_path, 
                                           "loc4_BUCB_squared_long.Rda"))
#load(file = paste0(current_path, "loc4_BUCB_squared_long.Rda"))
loc4_comp_BUCB_squared_long <- compare_to_ground_truth(mean_loc4, loc4_BUCB_squared_long, 
                                                  tau_loc4, epsilon_loc4)$mean
save(loc4_comp_BUCB_squared_long, file = paste0(current_path, "loc4_comp_BUCB_squared_long.Rda"))
rm(loc4_BUCB_squared_long)
gc()

########################################################################

loc4_BUCB_sqrt <- para_bandit_sim_bucb(data = data_list4, rounds = 2000, 
                                       rate = "inverse_sqrt",
                                       tau = tau_loc4, epsilon = epsilon_loc4, 
                                       alpha = tau_loc4, beta = 1-tau_loc4)
save(loc4_BUCB_sqrt, file = paste0(current_path, "loc4_BUCB_sqrt.Rda"))
#load(file = paste0(current_path, "loc4_BUCB_sqrt.Rda"))
loc4_comp_BUCB_sqrt <- compare_to_ground_truth(mean_loc4, loc4_BUCB_sqrt, 
                                               tau_loc4, epsilon_loc4)$mean
save(loc4_comp_BUCB_sqrt, file = paste0(current_path, "loc4_comp_BUCB_sqrt.Rda"))
rm(loc4_BUCB_sqrt)
gc()

########################################################################
# Do KL by comparing against tau with 0/1 adjustment

system.time(loc4_KL_tau_horizon <- para_bandit_sim_KL(data = data_list4, 
                                                      rounds = 2000, 
                                                      tau = tau_loc4, 
                                                      epsilon = epsilon_loc4, 
                                                      at_tau = TRUE, 
                                                      horizon = 2000^2))
save(loc4_KL_tau_horizon, file = paste0(current_path, "loc4_KL_tau_horizon.Rda"))
loc4_comp_KL_tau_horizon <- compare_to_ground_truth(mean_loc4, 
                                                    loc4_KL_tau_horizon, 
                                                    tau_loc4, epsilon_loc4)$mean
save(loc4_comp_KL_tau_horizon, file = paste0(current_path,
                                             "loc4_comp_KL_tau_horizon.Rda"))
rm(loc4_KL_tau_horizon)
gc()

########################################################################
# Do KL by comparing against tau and epsilon with 0/1 adjustment

system.time(loc4_KL_horizon1 <- para_bandit_sim_KL(data = data_list4, rounds = 2000, 
                                                   tau = tau_loc4, epsilon = epsilon_loc4, 
                                                   at_tau = FALSE, horizon = 2000))
save(loc4_KL_horizon1, file = paste0(current_path, "loc4_KL_horizon1.Rda"))
load(file = paste0(current_path, "loc4_KL_horizon1.Rda"))
loc4_comp_KL_horizon1 <- compare_to_ground_truth(mean_loc4, loc4_KL_horizon1, 
                                                 tau_loc4, epsilon_loc4)$mean
save(loc4_comp_KL_horizon1, file = paste0(current_path,
                                          "loc4_comp_KL_horizon1.Rda"))

########################################################################

loc4_TTS <- para_bandit_sim_TTS(data = data_list4, rounds = 2000,
                                tau = tau_loc4, epsilon = epsilon_loc4,
                                alpha = tau_loc4, beta = 1 - tau_loc4)
save(loc4_TTS, file = paste0(current_path, "loc4_TTS.Rda"))
load(file = paste0(current_path, "loc4_TTS.Rda"))
loc4_comp_TTS <- compare_to_ground_truth(mean_loc4, loc4_TTS, 
                                         tau_loc4, epsilon_loc4)$mean
save(loc4_comp_TTS, file = paste0(current_path, "loc4_comp_TTS.Rda"))
rm(loc4_TTS)
gc()

########################################################################
# Standard AugUCB
system.time(loc4_AugUCB <- para_bandit_sim_AugUCB(data = data_list4, 
                                                  rounds = 2000, 
                                                  tau = tau_loc4))
save(loc4_AugUCB, file = paste0(current_path, "loc4_AugUCB.Rda"))
load(file = paste0(current_path, "loc4_AugUCB.Rda"))
loc4_comp_AugUCB <- compare_to_ground_truth(mean_loc4, loc4_AugUCB, tau_loc4, 
                                            epsilon_loc4)$mean
save(loc4_comp_AugUCB, file = paste0(current_path, "loc4_comp_AugUCB.Rda"))

# now with 8000

system.time(loc4_AugUCB_8000 <- para_bandit_sim_AugUCB(data = data_list4_8000, 
                                                       rounds = 8000, 
                                                       tau = tau_loc4))
#user    system   elapsed 
#60.824    67.047 24179.679 
save(loc4_AugUCB_8000, file = paste0(current_path, "loc4_AugUCB_8000.Rda"))
#load(file = paste0(current_path, "loc4_AugUCB_8000.Rda"))
loc4_comp_AugUCB_8000 <- compare_to_ground_truth(mean_loc4, loc4_AugUCB_8000, 
                                                 tau_loc4, epsilon_loc4)$mean
save(loc4_comp_AugUCB_8000, file = paste0(current_path, 
                                          "loc4_comp_AugUCB_8000.Rda"))
########################################################################
# Do KL by comparing against tau directly

system.time(loc4_KL <- para_bandit_sim_KL(data = data_list4, rounds = 2000, 
                                          tau = tau_loc4, epsilon = epsilon_loc4, 
                                          at_tau = TRUE, horizon = 2000))
# user  system elapsed 
#6.824   4.335 391.240
save(loc4_KL, file = paste0(current_path, "loc4_KL.Rda"))
loc4_comp_KL <- compare_to_ground_truth(mean_loc4, loc4_KL, tau_loc4, 
                                        epsilon_loc4)$mean
save(loc4_comp_KL, file = paste0(current_path, "loc4_comp_KL.Rda"))

########################################################################
# Do KL by comparing against tau and epsilon

system.time(loc4_KL_not_tau <- para_bandit_sim_KL(data = data_list4, rounds = 2000, 
                                                  tau = tau_loc4, epsilon = epsilon_loc4, 
                                                  at_tau = FALSE))
save(loc4_KL_not_tau, file = paste0(current_path, "loc4_KL_not_tau.Rda"))
loc4_comp_KL_not_tau <- compare_to_ground_truth(mean_loc4, loc4_KL_not_tau, 
                                                tau_loc4, epsilon_loc4)$mean
save(loc4_comp_KL_not_tau, file = paste0(current_path,
                                         "loc4_comp_KL_not_tau.Rda"))

########################################################################
# Do KL by comparing against tau and epsilon with 0/1 adjustment

system.time(loc4_KL_horizon <- para_bandit_sim_KL(data = data_list4, rounds = 2000, 
                                                  tau = tau_loc4, epsilon = epsilon_loc4, 
                                                  at_tau = FALSE, horizon = 20000))
save(loc4_KL_horizon, file = paste0(current_path, "loc4_KL_horizon.Rda"))
loc4_comp_KL_horizon <- compare_to_ground_truth(mean_loc4, loc4_KL_horizon, 
                                                tau_loc4, epsilon_loc4)$mean
save(loc4_comp_KL_horizon, file = paste0(current_path,
                                         "loc4_comp_KL_horizon.Rda"))


########################################################################
# Do KL by comparing against tau and epsilon with 0/1 adjustment

system.time(loc4_KL_horizon2 <- para_bandit_sim_KL(data = data_list4, rounds = 2000, 
                                                   tau = tau_loc4, epsilon = epsilon_loc4, 
                                                   at_tau = FALSE, horizon = 2000^2))
save(loc4_KL_horizon2, file = paste0(current_path, "loc4_KL_horizon2.Rda"))
loc4_comp_KL_horizon2 <- compare_to_ground_truth(mean_loc4, loc4_KL_horizon2, 
                                                 tau_loc4, epsilon_loc4)$mean
save(loc4_comp_KL_horizon2, file = paste0(current_path,
                                          "loc4_comp_KL_horizon2.Rda"))

########################################################################
# Probability that arm is above or below tau+-epsilon

loc4_PI <- para_bandit_sim_PI(data = data_list4, rounds = 2000,
                              tau = tau_loc4, epsilon = epsilon_loc4,
                              alpha = tau_loc4, beta = 1 - tau_loc4)
save(loc4_PI, file = paste0(current_path, "loc4_PI.Rda"))
load(file = paste0(current_path, "loc4_PI.Rda"))
loc4_comp_PI <- compare_to_ground_truth(mean_loc4, loc4_PI, 
                                        tau_loc4, epsilon_loc4)$mean
save(loc4_comp_PI, file = paste0(current_path, "loc4_comp_PI.Rda"))
rm(loc4_PI)
gc()

########################################################################
# Probability that arm is above or below tau+-epsilon

loc4_PI_tadj <- para_bandit_sim_PI(data = data_list4, rounds = 2000,
                                   tau = tau_loc4, epsilon = epsilon_loc4,
                                   alpha = tau_loc4, beta = 1 - tau_loc4, tadj = TRUE)
save(loc4_PI_tadj, file = paste0(current_path, "loc4_PI_tadj.Rda"))
load(file = paste0(current_path, "loc4_PI_tadj.Rda"))
loc4_comp_PI_tadj <- compare_to_ground_truth(mean_loc4, loc4_PI_tadj, 
                                             tau_loc4, epsilon_loc4)$mean
save(loc4_comp_PI_tadj, file = paste0(current_path, "loc4_comp_PI_tadj.Rda"))
rm(loc4_PI_tadj)
gc()
########################################################################

#system.time(loc4_BETA <- para_bandit_sim_BETA(data = data_list4, rounds = 1000,
#                                              tau = tau_loc4, epsilon = epsilon_loc4,
#                                              alpha = tau_loc4, beta = 1 - tau_loc4))
#save(loc4_BETA, file = paste0(current_path, "loc4_BETA.Rda"))
#loc4_comp_BETA <- compare_to_ground_truth(mean_loc4, loc4_BETA, 
#                                          tau_loc4, epsilon_loc4)$mean
#save(loc4_comp_BETA, file = paste0(current_path, "loc4_comp_BETA.Rda"))

########################################################################

loc4_BUCB_cubic <- para_bandit_sim_bucb(data = data_list4, rounds = 2000, 
                                        rate = "inverse_cubic",
                                        tau = tau_loc4, epsilon = epsilon_loc4, 
                                        alpha = tau_loc4, beta = 1-tau_loc4)
save(loc4_BUCB_cubic, file = paste0(current_path, "loc4_BUCB_cubic.Rda"))
loc4_comp_BUCB_cubic <- compare_to_ground_truth(mean_loc4, loc4_BUCB_cubic, 
                                                tau_loc4, epsilon_loc4)$mean
save(loc4_comp_BUCB_cubic, file = paste0(current_path, "loc4_comp_BUCB_cubic.Rda"))

########################################################################

loc4_BUCB_power5 <- para_bandit_sim_bucb(data = data_list4, rounds = 2000, 
                                         rate = "inverse_power5",
                                         tau = tau_loc4, epsilon = epsilon_loc4, 
                                         alpha = tau_loc4, beta = 1-tau_loc4)
save(loc4_BUCB_power5, file = paste0(current_path, "loc4_BUCB_power5.Rda"))
loc4_comp_BUCB_power5 <- compare_to_ground_truth(mean_loc4, loc4_BUCB_power5, 
                                                 tau_loc4, epsilon_loc4)$mean
save(loc4_comp_BUCB_power5, file = paste0(current_path, "loc4_comp_BUCB_power5.Rda"))

########################################################################

loc4_BUCB_squared_no_eps <- para_bandit_sim_bucb(data = data_list4, rounds = 2000, 
                                                 rate = "inverse_squared",
                                                 tau = tau_loc4, epsilon = epsilon_loc4, 
                                                 alpha = tau_loc4, beta = 1-tau_loc4,
                                                 with_epsilon = FALSE)
save(loc4_BUCB_squared_no_eps, 
     file = paste0(current_path, "loc4_BUCB_squared_no_eps.Rda"))
loc4_comp_BUCB_squared_no_eps <- compare_to_ground_truth(mean_loc4, 
                                                         loc4_BUCB_squared_no_eps,
                                                         tau_loc4,
                                                         epsilon_loc4)$mean
save(loc4_comp_BUCB_squared_no_eps, file = paste0(current_path, "loc4_comp_BUCB_squared_no_eps.Rda"))

########################################################################

loc4_BUCB_horizon5 <- para_bandit_sim_bucb(data = data_list4, rounds = 2000, 
                                           rate = "inverse_horizon_c",
                                           tau = tau_loc4, epsilon = epsilon_loc4, 
                                           alpha = tau_loc4, beta = 1-tau_loc4,
                                           with_epsilon = FALSE, const = 5)
save(loc4_BUCB_horizon5, 
     file = paste0(current_path, "loc4_BUCB_horizon5.Rda"))
loc4_comp_BUCB_horizon5 <- compare_to_ground_truth(mean_loc4, 
                                                   loc4_BUCB_horizon5,
                                                   tau_loc4,
                                                   epsilon_loc4)$mean
save(loc4_comp_BUCB_horizon5, file = paste0(current_path, "loc4_comp_BUCB_horizon5.Rda"))

########################################################################

load(file = paste0(current_path, "loc4_comp_UNIFORM.Rda"))
load(file = paste0(current_path, "loc4_comp_APT.Rda"))
load(file = paste0(current_path, "loc4_comp_AugUCB.Rda"))
load(file = paste0(current_path, "loc4_comp_AugUCB_8000.Rda"))
load(file = paste0(current_path, "loc4_comp_KL_tau_horizon.Rda"))
load(file = paste0(current_path, "loc4_comp_TTS.Rda"))
load(file = paste0(current_path, "loc4_comp_BUCB_sqrt.Rda"))
load(file = paste0(current_path, "loc4_comp_BUCB.Rda"))
load(file = paste0(current_path, "loc4_comp_BUCB_horizon.Rda"))
load(file = paste0(current_path, "loc4_comp_BUCB_horizon5.Rda"))
load(file = paste0(current_path, "loc4_comp_BUCB_squared.Rda"))

load(paste0(current_path, "loc4_comp_BUCB_horizon_8000.Rda"))
load(file = paste0(current_path, "loc4_comp_BUCB_horizon_long.Rda"))
load(file = paste0(current_path, "loc4_comp_BUCB_squared_long.Rda"))
load(file = paste0(current_path, "loc4_comp_BUCB_linearhorizon_long.Rda"))
load(paste0(current_path, "loc4_comp_BUCB_linearhorizonC_long.Rda"))
load(paste0(current_path, "loc4_comp_BUCB_horizon_uninfprior_long.Rda"))
load(paste0(current_path, "loc4_comp_BUCB_linearhorizonC_uninf_long.Rda"))

load(paste0(current_path, "loc4_comp_KLUCB12.Rda"))
load(paste0(current_path, "loc4_comp_KLUCB70.Rda"))
load(paste0(current_path, "loc4_comp_KLUCB90.Rda"))
load(paste0(current_path, "loc4_comp_KLUCB24100.Rda"))
load(paste0(current_path, "loc4_comp_AugUCB.Rda"))
load(paste0(current_path, "loc4_comp_AugUCB_8000.Rda"))
load(paste0(current_path, "loc4_KLUCB90_long.Rda"))

plot(c(0,8000), c(0, -5), type = "n")
lines(log(loc4_comp_BUCB_horizon_8000), col = "red")
lines(log(loc4_comp_BUCB_horizon_long), col = "red")
lines(log(loc4_comp_BUCB_horizon_uninfprior_long), col = "darkred")
lines(log(loc4_comp_BUCB_squared_long), col = "pink")
lines(log(loc4_comp_KLUCB180), col = "blue")
lines(log(loc4_comp_KLUCB90), col = "darkblue")
lines(log(loc4_comp_KLUCB90), col = "darkblue")
lines(log(loc4_comp_APT), col = "green")
lines(log(loc4_comp_AugUCB), col = "grey")
lines(log(loc4_comp_AugUCB_8000), col = "black")
lines(log(loc4_comp_BUCB_linearhorizon_long), col = "blue")
lines(log(loc4_comp_BUCB_linearhorizonC_long), col = "green")
lines(log(loc4_comp_BUCB_linearhorizonC_uninf_long), col = "darkgreen")
lines(log(loc4_comp_KLUCB90), col = "black")
lines(log(loc4_comp_KLUCB90_long), col = "black")
lines(log(loc4_comp_KLUCB24100), col = "blue")
lines(log(loc4_comp_KLUCB12), col = "blue")
lines(log(loc4_comp_KLUCB70), col = "blue")
lines(log(loc4_comp_KLUCB90_2), col = "orange")
lines(log(loc4_comp_KLUCB110), col = "blue")

########################################################################

library(plyr)
loc4_UNIFORM_as <- colMeans(ldply(loc4_UNIFORM, function(x) table(x$arm_sequence)))
loc4_APT_as <- colMeans(ldply(loc4_APT, function(x) table(x$arm_sequence)))
loc4_PI_as <- colMeans(ldply(loc4_PI, function(x) table(x$arm_sequence)))
loc4_PI_tadj_as <- colMeans(ldply(loc4_PI_tadj, function(x) table(x$arm_sequence)))
loc4_TTS_as <- colMeans(ldply(loc4_TTS, function(x) table(x$arm_sequence)))
loc4_BUCB_as <- colMeans(ldply(loc4_BUCB, function(x) table(x$arm_sequence)))
loc4_BUCB_squared_as <- colMeans(ldply(loc4_BUCB_squared, function(x) table(x$arm_sequence)))
loc4_BUCB_cubic_as <- colMeans(ldply(loc4_BUCB_cubic, function(x) table(x$arm_sequence)))
loc4_KL_as <- colMeans(ldply(loc4_KL, function(x) table(x$arm_sequence)))
loc4_KL_not_tau_as <- colMeans(ldply(loc4_KL_not_tau, function(x) table(x$arm_sequence)))
loc4_KL_horizon1_as <- colMeans(ldply(loc4_KL_horizon1, function(x) table(x$arm_sequence)))

round(data.frame(loc4_UNIFORM_as, loc4_APT_as, 
                 loc4_BUCB_as, loc4_BUCB_squared_as,
                 loc4_TTS_as#, loc4_KL_as, loc4_KL_horizon1_as
))


loc4_UNIFORM_as_1 <- colMeans(ldply(loc4_UNIFORM, function(x) table(x$arm_sequence[1:500])))
loc4_APT_as_1 <- colMeans(ldply(loc4_APT, function(x) table(x$arm_sequence[1:500])))
loc4_PI_as_1 <- colMeans(ldply(loc4_PI, function(x) table(x$arm_sequence[1:500])))
loc4_PI_tadj_as_1 <- colMeans(ldply(loc4_PI_tadj, function(x) table(x$arm_sequence[1:500])))
loc4_TTS_as_1 <- colMeans(ldply(loc4_TTS, function(x) table(x$arm_sequence[1:500])))
loc4_BUCB_as_1 <- colMeans(ldply(loc4_BUCB, function(x) table(x$arm_sequence[1:500])))
loc4_KL_as_1 <- colMeans(ldply(loc4_KL, function(x) table(x$arm_sequence[1:500])))
loc4_KL_not_tau_as_1 <- colMeans(ldply(loc4_KL_not_tau, function(x) table(x$arm_sequence[1:500])))
loc4_KL_horizon1_as_1 <- colMeans(ldply(loc4_KL_horizon1, function(x) table(x$arm_sequence[1:500])))

round(data.frame(loc4_UNIFORM_as_1, loc4_APT_as_1, loc4_BUCB_as_1,
                 loc4_PI_as_1, loc4_PI_tadj_as_1, loc4_TTS_as_1,
                 loc4_KL_as_1, loc4_KL_horizon1_as_1))/500


loc4_UNIFORM_as_2 <- colMeans(ldply(loc4_UNIFORM, function(x) table(x$arm_sequence[501:1000])))
loc4_APT_as_2 <- colMeans(ldply(loc4_APT, function(x) table(x$arm_sequence[501:1000])))
loc4_PI_as_2 <- colMeans(ldply(loc4_PI, function(x) table(x$arm_sequence[501:1000])))
loc4_PI_tadj_as_2 <- colMeans(ldply(loc4_PI_tadj, function(x) table(x$arm_sequence[501:1000])))
loc4_TTS_as_2 <- colMeans(ldply(loc4_TTS, function(x) table(x$arm_sequence[501:1000])))
loc4_BUCB_as_2 <- colMeans(ldply(loc4_BUCB, function(x) table(x$arm_sequence[501:1000])))
loc4_KL_as_2 <- colMeans(ldply(loc4_KL, function(x) table(x$arm_sequence[501:1000])))
loc4_KL_not_tau_as_2 <- colMeans(ldply(loc4_KL_not_tau, function(x) table(x$arm_sequence[501:1000])))
loc4_KL_horizon1_as_2 <- colMeans(ldply(loc4_KL_horizon1, function(x) table(x$arm_sequence[501:1000])))

round(data.frame(loc4_UNIFORM_as, loc4_APT_as, loc4_BUCB_as,
                 loc4_PI_as, loc4_PI_tadj_as, loc4_TTS_as,
                 loc4_KL_as, loc4_KL_horizon1_as))

########################################################################

plot(loc4_BUCB_squared[[1]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc4_BUCB_squared[[2]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc4_BUCB_squared[[3]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc4_BUCB_squared[[4]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc4_BUCB_squared[[5]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc4_BUCB_squared[[6]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc4_BUCB_squared[[7]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc4_BUCB_squared[[8]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc4_BUCB_squared[[9]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc4_BUCB_squared[[10]]$arm_sequence, pch = 19, cex = 0.3)

plot(loc4_APT[[1]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc4_APT[[2]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc4_APT[[3]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc4_APT[[4]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc4_APT[[5]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc4_APT[[6]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc4_APT[[7]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc4_APT[[8]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc4_APT[[9]]$arm_sequence, pch = 19, cex = 0.3)
plot(loc4_APT[[10]]$arm_sequence, pch = 19, cex = 0.3)

library(plyr)
library(tidyr)
library(dplyr)
#arm_seq_table <- arm_seq_res %>% tbl_df() %>% mutate(index = 1:1000) %>%
#  gather(key = iter, value = arm, -index) %>%
#  group_by(index, arm) %>% summarize(count = n())

library(ggplot2)
arm_seq_res_APT <- data.frame(t(laply(loc4_APT, function(x) x$arm_sequence)))
arm_seq_res_APT %>% tbl_df() %>% mutate(index = 1:2000) %>%
  gather(key = iter, value = arm, -index) %>%
  filter(index > 11) %>%
  ggplot(aes(index, arm)) + geom_count() +
  scale_size_area(max_size = 10)
#ggplot(arm_seq_table[11:9910,], aes(index, arm)) + geom_tile(aes(fill = count))

arm_seq_res_BUCB <- data.frame(t(laply(loc4_BUCB, function(x) x$arm_sequence)))
arm_seq_res_BUCB %>% tbl_df() %>% mutate(index = 1:2000) %>%
  gather(key = iter, value = arm, -index) %>%
  filter(index > 11) %>%
  ggplot(aes(index, arm)) + geom_count() +
  scale_size_area(max_size = 10)

arm_seq_res_BUCB_squared <- data.frame(t(laply(loc4_BUCB_squared, function(x) x$arm_sequence)))
arm_seq_res_BUCB_squared %>% tbl_df() %>% mutate(index = 1:2000) %>%
  gather(key = iter, value = arm, -index) %>%
  filter(index > 11) %>%
  ggplot(aes(index, arm)) + geom_count() +
  scale_size_area(max_size = 10)

arm_seq_res_TTS <- data.frame(t(laply(loc4_TTS, function(x) x$arm_sequence)))
arm_seq_res_TTS %>% tbl_df() %>% mutate(index = 1:2000) %>%
  gather(key = iter, value = arm, -index) %>%
  filter(index > 11) %>%
  ggplot(aes(index, arm)) + geom_count() +
  scale_size_area(max_size = 10)



load(paste0(current_path, "loc4_BUCB_horizon_long.Rda"))
arm_seq_res_BUCB_horizon <- data.frame(t(laply(loc4_BUCB_horizon_long, function(x) x$arm_sequence)))
arm_seq_res_BUCB_horizon %>% tbl_df() %>% mutate(index = 1:6000) %>%
  gather(key = iter, value = arm, -index) %>%
  filter(index > 11) %>%
  ggplot(aes(index, arm)) + geom_count() +
  scale_size_area(max_size = 10) + labs(title = "BUCB Log Horizon")

load(paste0(current_path, "loc4_BUCB_linearhorizon_long.Rda"))
arm_seq_res_BUCB_linearhorizon_long <- data.frame(t(laply(loc4_BUCB_linearhorizon_long, function(x) x$arm_sequence)))
arm_seq_res_BUCB_linearhorizon_long %>% tbl_df() %>% mutate(index = 1:6000) %>%
  gather(key = iter, value = arm, -index) %>%
  filter(index > 11) %>%
  ggplot(aes(index, arm)) + geom_count() +
  scale_size_area(max_size = 10) + labs(title = "BUCB Linear Horizon")

load(paste0(current_path, "loc4_BUCB_linearhorizonC_long.Rda"))
arm_seq_res_BUCB_linearhorizonC_long <- data.frame(t(laply(loc4_BUCB_linearhorizonC_long, function(x) x$arm_sequence)))
arm_seq_res_BUCB_linearhorizonC_long %>% tbl_df() %>% mutate(index = 1:6000) %>%
  gather(key = iter, value = arm, -index) %>%
  filter(index > 11) %>%
  ggplot(aes(index, arm)) + geom_count() +
  scale_size_area(max_size = 10) + labs(title = "BUCB Linear Horizon with C")

arm_count_horizon <- arm_seq_res_BUCB_horizon %>% tbl_df() %>% mutate(index = 1:6000) %>%
  gather(key = iter, value = arm, -index) %>%
  filter(index > 11)
arm_count_linearhorizon <- arm_seq_res_BUCB_linearhorizon_long %>% tbl_df() %>% mutate(index = 1:6000) %>%
  gather(key = iter, value = arm, -index) %>%
  filter(index > 11)
arm_count_linearhorizonC <- arm_seq_res_BUCB_linearhorizonC_long %>% tbl_df() %>% mutate(index = 1:6000) %>%
  gather(key = iter, value = arm, -index) %>%
  filter(index > 11)

rbind(mean_loc4,
      round(table(arm_count_horizon$arm)/5000),
      round(table(arm_count_linearhorizon$arm)/5000),
      round(table(arm_count_linearhorizonC$arm)/5000))
