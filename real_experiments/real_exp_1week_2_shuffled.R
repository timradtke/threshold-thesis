####################################################################

# This script loads real data for 10 products out of 197 products.
# We look at 10080 rounds per iteration so as to cover a whole week.
# This is useful to account for seasonality effects in the data (if
# looking at uniformly sampled data).

####################################################################

sessions_threshold <- 0
set.seed(214237)
# shuffle the first day of data
shuffle_idx <- sample(1:(24*60), 24*60)

set.seed(357458)
selected_products <- sample(1:197, 10)
niter <- 5000
nrounds <- 10080 # These two values give the maximum size data 
nahead <- 10080  # possible given the data set size we have
source("/Users/timradtke/Desktop/thesis_data/get_dataset_shuffled.R")

dim(pv_list[[1]])
pv_list[[1]]

pv_list[[1]]
pv_list[[5000]]
tail(pv_list[[5000]])
round(pv_next_means[[5000]],4)
round(pv_next_means[[1]],4)

####################################################################

current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/AugUCB/"
source(paste0(current_path, "run_from_data.R"))
source(paste0(current_path, "run_parallelized_from_data.R"))
current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/real_experiments/"

####################################################################

data_amo4_shuffled <- pv_list
data_amo4_shuffled_wide <- pv_products_wide
data_amo4_shuffled_own_means <- pv_own_means
data_amo4_shuffled_next_means <- pv_next_means
data_amo4_shuffled_mean_secondhalf <- colMeans(data_amo4_shuffled_wide[15081:25160,])
data_amo4_shuffled_mean_middlehalf <- colMeans(data_amo4_shuffled_wide[10081:20160,])
data_amo4_shuffled_mean_firsthalf <- colMeans(data_amo4_shuffled_wide[1:15080,])
save(data_amo4_shuffled, file = paste0(current_path, "data_amo4_shuffled.Rda"))
save(data_amo4_shuffled_wide, file = paste0(current_path, "data_amo4_shuffled_wide.Rda"))
save(data_amo4_shuffled_own_means, file = paste0(current_path, "data_amo4_shuffled_own_means.Rda"))
save(data_amo4_shuffled_next_means, file = paste0(current_path, "data_amo4_shuffled_next_means.Rda"))
save(data_amo4_shuffled_mean_firsthalf, file = paste0(current_path, "data_amo4_shuffled_mean_firsthalf.Rda"))
save(data_amo4_shuffled_mean_secondhalf, file = paste0(current_path, "data_amo4_shuffled_mean_secondhalf.Rda"))
rm(pv_list, pv_own_means, pv_next_means, pv_products_wide)
gc()

load(paste0(current_path, "data_amo4_shuffled_mean_firsthalf.Rda"))
load(paste0(current_path, "data_amo4_shuffled_mean_secondhalf.Rda"))
load(paste0(current_path, "data_amo4_shuffled_own_means.Rda"))
load(paste0(current_path, "data_amo4_shuffled.Rda"))
tau_amo4_shuffled <- 5/60
epsilon_amo4_shuffled <- 0

########################################################################
# Standard Likelihood Ratio
system.time(amo4_shuffled_LR <- para_bandit_sim_LR(data = data_amo4_shuffled, 
                                          rounds = 10080, 
                                          tau = tau_amo4_shuffled, 
                                          epsilon = epsilon_amo4_shuffled))

save(amo4_shuffled_LR, file = paste0(current_path, "amo4_shuffled_LR.Rda"))
#load(file = paste0(current_path, "amo4_shuffled_LR.Rda"))
amo4_shuffled_comp_LR <- compare_to_ground_truth(data_amo4_shuffled_mean_firsthalf, 
                                        amo4_shuffled_LR, 
                                        tau_amo4_shuffled, epsilon_amo4_shuffled)$mean
amo4_shuffled_compholdout_LR <- compare_to_ground_truth(data_amo4_shuffled_mean_secondhalf, 
                                               amo4_shuffled_LR, 
                                               tau_amo4_shuffled, epsilon_amo4_shuffled)$mean
amo4_shuffled_compown_LR <- compare_to_cv_data(data_amo4_shuffled_own_means, amo4_shuffled_LR, 
                                      tau_amo4_shuffled, epsilon_amo4_shuffled)$mean
amo4_shuffled_simpleregret_LR <- get_simple_regret(data_amo4_shuffled_mean_firsthalf, 
                                          amo4_shuffled_LR, 
                                          tau_amo4_shuffled, epsilon_amo4_shuffled)$mean
save(amo4_shuffled_comp_LR, file = paste0(current_path, "amo4_shuffled_comp_LR.Rda"))
save(amo4_shuffled_compholdout_LR, file = paste0(current_path, 
                                        "amo4_shuffled_compholdout_LR.Rda"))
save(amo4_shuffled_compown_LR, file = paste0(current_path, 
                                    "amo4_shuffled_compown_LR.Rda"))
save(amo4_shuffled_simpleregret_LR, file = paste0(current_path, 
                                         "amo4_shuffled_simpleregret_LR.Rda"))
rm(amo4_shuffled_LR, amo4_shuffled_comp_LR)
gc()

########################################################################
# Empirical Variance Guided Algorithm (Zhong et al., 2017)

system.time(amo4_shuffled_EVT <- para_bandit_sim_EVT(data = data_amo4_shuffled, 
                                            rounds = 10080, 
                                            tau = tau_amo4_shuffled, 
                                            epsilon = epsilon_amo4_shuffled))
save(amo4_shuffled_EVT, file = paste0(current_path, "amo4_shuffled_EVT.Rda"))
#load(file = paste0(current_path, "amo4_shuffled_EVT.Rda"))
amo4_shuffled_comp_EVT <- compare_to_ground_truth(data_amo4_shuffled_mean_firsthalf, 
                                         amo4_shuffled_EVT, 
                                         tau_amo4_shuffled, epsilon_amo4_shuffled)$mean
amo4_shuffled_compholdout_EVT <- compare_to_ground_truth(data_amo4_shuffled_mean_secondhalf, 
                                                amo4_shuffled_EVT, 
                                                tau_amo4_shuffled, epsilon_amo4_shuffled)$mean
amo4_shuffled_simpleregret_EVT <- get_simple_regret(data_amo4_shuffled_mean_firsthalf, 
                                           amo4_shuffled_EVT, 
                                           tau_amo4_shuffled, epsilon_amo4_shuffled)$mean
save(amo4_shuffled_comp_EVT, file = paste0(current_path, "amo4_shuffled_comp_EVT.Rda"))
save(amo4_shuffled_compholdout_EVT, file = paste0(current_path, 
                                         "amo4_shuffled_compholdout_EVT.Rda"))
save(amo4_shuffled_simpleregret_EVT, file = paste0(current_path, 
                                          "amo4_shuffled_simpleregret_EVT.Rda"))
rm(amo4_shuffled_EVT, amo4_shuffled_comp_EVT, amo4_shuffled_compholdout_EVT)
gc()

########################################################################
# Standard AugUCB
system.time(amo4_shuffled_AugUCB <- para_bandit_sim_AugUCB(data = data_amo4_shuffled, 
                                                  rounds = 10080, 
                                                  tau = tau_amo4_shuffled))
save(amo4_shuffled_AugUCB, file = paste0(current_path, "amo4_shuffled_AugUCB.Rda"))
load(file = paste0(current_path, "amo4_shuffled_AugUCB.Rda"))

amo4_shuffled_comp_AugUCB <- compare_to_ground_truth(data_amo4_shuffled_mean_firsthalf, 
                                            amo4_shuffled_AugUCB, 
                                            tau_amo4_shuffled, epsilon_amo4_shuffled)$mean
amo4_shuffled_compholdout_AugUCB <- compare_to_ground_truth(data_amo4_shuffled_mean_secondhalf, 
                                                   amo4_shuffled_AugUCB, 
                                                   tau_amo4_shuffled, epsilon_amo4_shuffled)$mean
amo4_shuffled_simpleregret_AugUCB <- get_simple_regret(data_amo4_shuffled_mean_firsthalf, 
                                              amo4_shuffled_AugUCB, 
                                              tau_amo4_shuffled, epsilon_amo4_shuffled)$mean
save(amo4_shuffled_comp_AugUCB, file = paste0(current_path, "amo4_shuffled_comp_AugUCB.Rda"))
save(amo4_shuffled_compholdout_AugUCB, file = paste0(current_path, 
                                            "amo4_shuffled_compholdout_AugUCB.Rda"))
save(amo4_shuffled_simpleregret_AugUCB, file = paste0(current_path, 
                                             "amo4_shuffled_simpleregret_AugUCB.Rda"))
rm(amo4_shuffled_AugUCB, amo4_shuffled_comp_AugUCB, amo4_shuffled_compholdout_AugUCB)
gc()

########################################################################
# Bayes-UCB

amo4_shuffled_BUCB <- para_bandit_sim_bucb(data = data_amo4_shuffled, rounds = 10080, 
                                  rate = "inverse_horizon",
                                  tau = tau_amo4_shuffled, epsilon = epsilon_amo4_shuffled, 
                                  alpha = tau_amo4_shuffled, beta = 1-tau_amo4_shuffled)
save(amo4_shuffled_BUCB, 
     file = paste0(current_path, "amo4_shuffled_BUCB.Rda"))
load(file = paste0(current_path, "amo4_shuffled_BUCB.Rda"))
amo4_shuffled_comp_BUCB <- compare_to_ground_truth(data_amo4_shuffled_mean_firsthalf, amo4_shuffled_BUCB, 
                                          tau_amo4_shuffled, epsilon_amo4_shuffled)$mean
amo4_shuffled_compholdout_BUCB <- compare_to_ground_truth(data_amo4_shuffled_mean_secondhalf, 
                                                 amo4_shuffled_BUCB, 
                                                 tau_amo4_shuffled, epsilon_amo4_shuffled)$mean
amo4_shuffled_simpleregret_BUCB <- get_simple_regret(data_amo4_shuffled_mean_firsthalf, 
                                            amo4_shuffled_BUCB, 
                                            tau_amo4_shuffled, epsilon_amo4_shuffled)$mean
save(amo4_shuffled_comp_BUCB, file = paste0(current_path, "amo4_shuffled_comp_BUCB.Rda"))
save(amo4_shuffled_compholdout_BUCB, file = paste0(current_path, "amo4_shuffled_compholdout_BUCB.Rda"))
save(amo4_shuffled_simpleregret_BUCB, file = paste0(current_path, "amo4_shuffled_simpleregret_BUCB.Rda"))
rm(amo4_shuffled_BUCB)
gc()

########################################################################
# Standard APT Algorithm
system.time(amo4_shuffled_APT <- para_bandit_sim_APT(data = data_amo4_shuffled, rounds = 10080, 
                                            tau = tau_amo4_shuffled, epsilon = epsilon_amo4_shuffled))
save(amo4_shuffled_APT, file = paste0(current_path, "amo4_shuffled_APT.Rda"))
#load(paste0(current_path, "amo4_shuffled_APT.Rda"))
amo4_shuffled_comp_APT <- compare_to_ground_truth(data_amo4_shuffled_mean_firsthalf, amo4_shuffled_APT, 
                                         tau_amo4_shuffled, epsilon_amo4_shuffled)$mean
amo4_shuffled_compholdout_APT <- compare_to_ground_truth(data_amo4_shuffled_mean_secondhalf, amo4_shuffled_APT, 
                                                tau_amo4_shuffled, epsilon_amo4_shuffled)$mean
amo4_shuffled_simpleregret_APT <- get_simple_regret(data_amo4_shuffled_mean_firsthalf, 
                                           amo4_shuffled_APT, 
                                           tau_amo4_shuffled, epsilon_amo4_shuffled)$mean
save(amo4_shuffled_comp_APT, file = paste0(current_path, "amo4_shuffled_comp_APT.Rda"))
save(amo4_shuffled_compholdout_APT, file = paste0(current_path, "amo4_shuffled_compholdout_APT.Rda"))
save(amo4_shuffled_simpleregret_APT, file = paste0(current_path, "amo4_shuffled_simpleregret_APT.Rda"))
rm(amo4_shuffled_APT)
gc()

########################################################################
# Standard Uniform
system.time(amo4_shuffled_UNIFORM <- para_bandit_sim_uniform(data = data_amo4_shuffled, 
                                                    rounds = 10080))
save(amo4_shuffled_UNIFORM, file = paste0(current_path, "amo4_shuffled_UNIFORM.Rda"))
#load(file = paste0(current_path, "amo4_shuffled_UNIFORM.Rda"))

amo4_shuffled_comp_UNIFORM <- compare_to_ground_truth(data_amo4_shuffled_mean_firsthalf, 
                                             amo4_shuffled_UNIFORM, 
                                             tau_amo4_shuffled, epsilon_amo4_shuffled)$mean
amo4_shuffled_compholdout_UNIFORM <- compare_to_ground_truth(data_amo4_shuffled_mean_secondhalf, 
                                                    amo4_shuffled_UNIFORM, 
                                                    tau_amo4_shuffled, epsilon_amo4_shuffled)$mean
amo4_shuffled_simpleregret_UNIFORM <- get_simple_regret(data_amo4_shuffled_mean_firsthalf, 
                                               amo4_shuffled_UNIFORM, 
                                               tau_amo4_shuffled, epsilon_amo4_shuffled)$mean

save(amo4_shuffled_comp_UNIFORM, file = paste0(current_path, "amo4_shuffled_comp_UNIFORM.Rda"))
save(amo4_shuffled_compholdout_UNIFORM, file = paste0(current_path, 
                                             "amo4_shuffled_compholdout_UNIFORM.Rda"))
save(amo4_shuffled_simpleregret_UNIFORM, file = paste0(current_path, 
                                              "amo4_shuffled_simpleregret_UNIFORM.Rda"))
rm(amo4_shuffled_UNIFORM)
gc()

########################################################################
# Standard KL-UCB
Hholdout_amo4_shuffled <- get_complexity(data_amo4_shuffled_mean_secondhalf, tau_amo4_shuffled, epsilon_amo4_shuffled)
H_amo4_shuffled <- get_complexity(data_amo4_shuffled_mean_firsthalf, tau_amo4_shuffled, epsilon_amo4_shuffled)

system.time(amo4_shuffled_KLUCB <- para_bandit_sim_KLUCB(data = data_amo4_shuffled, 
                                                rounds = 7580, 
                                                tau = tau_amo4_shuffled, 
                                                epsilon = epsilon_amo4_shuffled,
                                                horizon = 5.5*H_amo4_shuffled,
                                                H = H_amo4_shuffled))
save(amo4_shuffled_KLUCB, file = paste0(current_path, "amo4_shuffled_KLUCB.Rda"))
#load(file = paste0(current_path, "amo4_shuffled_KLUCB.Rda"))

amo4_shuffled_comp_KLUCB <- compare_to_ground_truth(data_amo4_shuffled_mean_firsthalf, 
                                           amo4_shuffled_KLUCB, 
                                           tau_amo4_shuffled, epsilon_amo4_shuffled)$mean
amo4_shuffled_compholdout_KLUCB <- compare_to_ground_truth(data_amo4_shuffled_mean_secondhalf, 
                                                  amo4_shuffled_KLUCB, 
                                                  tau_amo4_shuffled, epsilon_amo4_shuffled)$mean
save(amo4_shuffled_comp_KLUCB, file = paste0(current_path, "amo4_shuffled_comp_KLUCB.Rda"))
save(amo4_shuffled_compholdout_KLUCB, file = paste0(current_path, 
                                           "amo4_shuffled_compholdout_KLUCB.Rda"))
rm(amo4_shuffled_KLUCB)
gc()

########################################################################

amo4_shuffled_TTS <- para_bandit_sim_TTS(data = data_amo4_shuffled, rounds = 7580,
                                tau = tau_amo4_shuffled, epsilon = epsilon_amo4_shuffled,
                                alpha = tau_amo4_shuffled, beta = 1 - tau_amo4_shuffled)
save(amo4_shuffled_TTS, file = paste0(current_path, "amo4_shuffled_TTS.Rda"))

amo4_shuffled_comp_TTS <- compare_to_ground_truth(data_amo4_shuffled_mean_firsthalf, 
                                         amo4_shuffled_TTS, 
                                         tau_amo4_shuffled, epsilon_amo4_shuffled)$mean
amo4_shuffled_compholdout_TTS <- compare_to_ground_truth(data_amo4_shuffled_mean_secondhalf, 
                                                amo4_shuffled_TTS, 
                                                tau_amo4_shuffled, epsilon_amo4_shuffled)$mean
save(amo4_shuffled_comp_TTS, file = paste0(current_path, "amo4_shuffled_comp_TTS.Rda"))
save(amo4_shuffled_compholdout_TTS, file = paste0(current_path, 
                                         "amo4_shuffled_compholdout_TTS.Rda"))
rm(amo4_shuffled_TTS)
gc()

########################################################################

round(data_amo4_shuffled_mean_firsthalf - tau_amo4_shuffled,2)

load(paste0(current_path, "amo4_shuffled_comp_UNIFORM.Rda"))
load(paste0(current_path, "amo4_shuffled_comp_APT.Rda"))
load(paste0(current_path, "amo4_shuffled_comp_LR.Rda"))
load(paste0(current_path, "amo4_shuffled_comp_AugUCB.Rda"))
load(paste0(current_path, "amo4_shuffled_comp_EVT.Rda"))
load(paste0(current_path, "amo4_shuffled_comp_BUCB.Rda"))

plot(c(0,10080), c(0, -1.5), type = "n")
lines(log(amo4_shuffled_comp_UNIFORM), col = "black")
lines(log(amo4_shuffled_comp_APT), col = "blue")
lines(log(amo4_shuffled_comp_AugUCB), col = "green")
lines(log(amo4_shuffled_comp_EVT), col = "darkgreen")
lines(log(amo4_shuffled_comp_LR), col = "red")
lines(log(amo4_shuffled_compown_LR), col = "blue")
lines(log(amo4_shuffled_comp_BUCB), col = "violet")
abline(h = log(0.1), lty = 2)

plot(c(0,10080), c(0, -8), type = "n")
lines(log(amo4_shuffled_simpleregret_UNIFORM), col = "black")
lines(log(amo4_shuffled_simpleregret_APT), col = "blue")
lines(log(amo4_shuffled_simpleregret_AugUCB), col = "green")
lines(log(amo4_shuffled_simpleregret_EVT), col = "darkgreen")
lines(log(amo4_shuffled_simpleregret_LR), col = "red")
lines(log(amo4_shuffled_simpleregret_BUCB), col = "violet")

plot(c(0,10080), c(0, 0.7), type = "n")
lines((amo4_shuffled_simpleregret_UNIFORM), col = "black")
lines((amo4_shuffled_simpleregret_APT), col = "blue")
lines((amo4_shuffled_simpleregret_AugUCB), col = "green")
lines((amo4_shuffled_simpleregret_EVT), col = "darkgreen")
lines((amo4_shuffled_simpleregret_LR), col = "red")
lines((amo4_shuffled_simpleregret_BUCB), col = "violet")
