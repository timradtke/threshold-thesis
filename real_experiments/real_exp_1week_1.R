####################################################################

# This script loads real data for 10 products out of 197 products.
# We look at 10080 rounds per iteration so as to cover a whole week.
# This is useful to account for seasonality effects in the data (if
# looking at uniformly sampled data).

####################################################################

sessions_threshold <- 0
set.seed(1068)
selected_products <- sample(1:197, 10)
niter <- 5000
nrounds <- 10080 # These two values give the maximum size data 
nahead <- 10080  # possible given the data set size we have
source("/Users/timradtke/Desktop/thesis_data/get_dataset.R")

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

data_amo3 <- pv_list
data_amo3_wide <- pv_products_wide
data_amo3_own_means <- pv_own_means
data_amo3_next_means <- pv_next_means
data_amo3_mean_secondhalf <- colMeans(data_amo3_wide[15081:25160,])
data_amo3_mean_middlehalf <- colMeans(data_amo3_wide[10081:20160,])
data_amo3_mean_firsthalf <- colMeans(data_amo3_wide[1:15080,])
save(data_amo3, file = paste0(current_path, "data_amo3.Rda"))
save(data_amo3_wide, file = paste0(current_path, "data_amo3_wide.Rda"))
save(data_amo3_own_means, file = paste0(current_path, "data_amo3_own_means.Rda"))
save(data_amo3_next_means, file = paste0(current_path, "data_amo3_next_means.Rda"))
save(data_amo3_mean_firsthalf, file = paste0(current_path, "data_amo3_mean_firsthalf.Rda"))
save(data_amo3_mean_secondhalf, file = paste0(current_path, "data_amo3_mean_secondhalf.Rda"))
rm(pv_list, pv_own_means, pv_next_means, pv_products_wide)
gc()

load(paste0(current_path, "data_amo3_mean_firsthalf.Rda"))
load(paste0(current_path, "data_amo3_mean_secondhalf.Rda"))
load(paste0(current_path, "data_amo3.Rda"))
tau_amo3 <- 2/60
epsilon_amo3 <- 0

########################################################################
# Bayes-UCB

amo3_BUCB <- para_bandit_sim_bucb(data = data_amo3, rounds = 10080, 
                                  rate = "inverse_horizon",
                                  tau = tau_amo3, epsilon = epsilon_amo3, 
                                  alpha = tau_amo3, beta = 1-tau_amo3)
save(amo3_BUCB, 
     file = paste0(current_path, "amo3_BUCB.Rda"))
amo3_comp_BUCB <- compare_to_ground_truth(data_amo3_mean_firsthalf, amo3_BUCB, 
                                          tau_amo3, epsilon_amo3)$mean
amo3_compholdout_BUCB <- compare_to_ground_truth(data_amo3_mean_secondhalf, 
                                                 amo3_BUCB, 
                                                 tau_amo3, epsilon_amo3)$mean
save(amo3_comp_BUCB, file = paste0(current_path, "amo3_comp_BUCB.Rda"))
save(amo3_compholdout_BUCB, file = paste0(current_path, "amo3_compholdout_BUCB.Rda"))
rm(amo3_BUCB)
gc()

########################################################################
# Standard Likelihood Ratio
system.time(amo3_LR <- para_bandit_sim_LR(data = data_amo3, 
                                          rounds = 10080, 
                                          tau = tau_amo3, 
                                          epsilon = epsilon_amo3))

save(amo3_LR, file = paste0(current_path, "amo3_LR.Rda"))
#load(file = paste0(current_path, "amo3_LR.Rda"))
amo3_comp_LR <- compare_to_ground_truth(data_amo3_mean_firsthalf, 
                                             amo3_LR, 
                                             tau_amo3, epsilon_amo3)$mean
amo3_compholdout_LR <- compare_to_ground_truth(data_amo3_mean_secondhalf, 
                                                    amo3_LR, 
                                                    tau_amo3, epsilon_amo3)$mean
save(amo3_comp_LR, file = paste0(current_path, "amo3_comp_LR.Rda"))
save(amo3_compholdout_LR, file = paste0(current_path, 
                                        "amo3_compholdout_LR.Rda"))
rm(amo3_LR, amo3_comp_LR)
gc()

########################################################################
# Empirical Variance Guided Algorithm (Zhong et al., 2017)

system.time(amo3_EVT <- para_bandit_sim_EVT(data = data_amo3, 
                                            rounds = 10080, 
                                            tau = tau_amo3, 
                                            epsilon = epsilon_amo3))
save(amo3_EVT, file = paste0(current_path, "amo3_EVT.Rda"))
amo3_comp_EVT <- compare_to_ground_truth(data_amo3_mean_firsthalf, 
                                        amo3_EVT, 
                                        tau_amo3, epsilon_amo3)$mean
amo3_compholdout_EVT <- compare_to_ground_truth(data_amo3_mean_secondhalf, 
                                               amo3_EVT, 
                                               tau_amo3, epsilon_amo3)$mean
save(amo3_comp_EVT, file = paste0(current_path, "amo3_comp_EVT.Rda"))
save(amo3_compholdout_EVT, file = paste0(current_path, 
                                        "amo3_compholdout_EVT.Rda"))
rm(amo3_EVT, amo3_comp_EVT, amo3_compholdout_EVT)
gc()

########################################################################
# Standard AugUCB
system.time(amo3_AugUCB <- para_bandit_sim_AugUCB(data = data_amo3, 
                                                  rounds = 10080, 
                                                  tau = tau_amo3))
save(amo3_AugUCB, file = paste0(current_path, "amo3_AugUCB.Rda"))
#load(file = paste0(current_path, "amo3_AugUCB.Rda"))

amo3_comp_AugUCB <- compare_to_ground_truth(data_amo3_mean_firsthalf, 
                                            amo3_AugUCB, 
                                            tau_amo3, epsilon_amo3)$mean
amo3_compholdout_AugUCB <- compare_to_ground_truth(data_amo3_mean_secondhalf, 
                                                   amo3_AugUCB, 
                                                   tau_amo3, epsilon_amo3)$mean
save(amo3_comp_AugUCB, file = paste0(current_path, "amo3_comp_AugUCB.Rda"))
save(amo3_compholdout_AugUCB, file = paste0(current_path, 
                                            "amo3_compholdout_AugUCB.Rda"))
rm(amo3_AugUCB, amo3_comp_AugUCB, amo3_compholdout_AugUCB)
gc()

########################################################################
# Standard APT Algorithm
system.time(amo3_APT <- para_bandit_sim_APT(data = data_amo3, rounds = 10080, 
                                            tau = tau_amo3, epsilon = epsilon_amo3))
save(amo3_APT, file = paste0(current_path, "amo3_APT.Rda"))
#load(paste0(current_path, "amo3_APT.Rda"))
amo3_comp_APT <- compare_to_ground_truth(data_amo3_mean_firsthalf, amo3_APT, 
                                         tau_amo3, epsilon_amo3)$mean
amo3_compholdout_APT <- compare_to_ground_truth(data_amo3_mean_secondhalf, amo3_APT, 
                                                tau_amo3, epsilon_amo3)$mean

save(amo3_comp_APT, file = paste0(current_path, "amo3_comp_APT.Rda"))
save(amo3_compholdout_APT, file = paste0(current_path, "amo3_compholdout_APT.Rda"))
rm(amo3_APT)
gc()

########################################################################
# Standard Uniform
system.time(amo3_UNIFORM <- para_bandit_sim_uniform(data = data_amo3, 
                                                    rounds = 10080))
save(amo3_UNIFORM, file = paste0(current_path, "amo3_UNIFORM.Rda"))
#load(file = paste0(current_path, "amo3_UNIFORM.Rda"))

amo3_comp_UNIFORM <- compare_to_ground_truth(data_amo3_mean_firsthalf, 
                                             amo3_UNIFORM, 
                                             tau_amo3, epsilon_amo3)$mean
amo3_compholdout_UNIFORM <- compare_to_ground_truth(data_amo3_mean_secondhalf, 
                                                    amo3_UNIFORM, 
                                                    tau_amo3, epsilon_amo3)$mean
save(amo3_comp_UNIFORM, file = paste0(current_path, "amo3_comp_UNIFORM.Rda"))
save(amo3_compholdout_UNIFORM, file = paste0(current_path, 
                                             "amo3_compholdout_UNIFORM.Rda"))
rm(amo3_UNIFORM)
gc()

########################################################################
# Standard KL-UCB
Hholdout_amo3 <- get_complexity(data_amo3_mean_secondhalf, tau_amo3, epsilon_amo3)
H_amo3 <- get_complexity(data_amo3_mean_firsthalf, tau_amo3, epsilon_amo3)

system.time(amo3_KLUCB <- para_bandit_sim_KLUCB(data = data_amo3, 
                                                rounds = 7580, 
                                                tau = tau_amo3, 
                                                epsilon = epsilon_amo3,
                                                horizon = 5.5*H_amo3,
                                                H = H_amo3))
save(amo3_KLUCB, file = paste0(current_path, "amo3_KLUCB.Rda"))
#load(file = paste0(current_path, "amo3_KLUCB.Rda"))

amo3_comp_KLUCB <- compare_to_ground_truth(data_amo3_mean_firsthalf, 
                                           amo3_KLUCB, 
                                           tau_amo3, epsilon_amo3)$mean
amo3_compholdout_KLUCB <- compare_to_ground_truth(data_amo3_mean_secondhalf, 
                                                  amo3_KLUCB, 
                                                  tau_amo3, epsilon_amo3)$mean
save(amo3_comp_KLUCB, file = paste0(current_path, "amo3_comp_KLUCB.Rda"))
save(amo3_compholdout_KLUCB, file = paste0(current_path, 
                                           "amo3_compholdout_KLUCB.Rda"))
rm(amo3_KLUCB)
gc()

########################################################################

amo3_TTS <- para_bandit_sim_TTS(data = data_amo3, rounds = 7580,
                                tau = tau_amo3, epsilon = epsilon_amo3,
                                alpha = tau_amo3, beta = 1 - tau_amo3)
save(amo3_TTS, file = paste0(current_path, "amo3_TTS.Rda"))

amo3_comp_TTS <- compare_to_ground_truth(data_amo3_mean_firsthalf, 
                                         amo3_TTS, 
                                         tau_amo3, epsilon_amo3)$mean
amo3_compholdout_TTS <- compare_to_ground_truth(data_amo3_mean_secondhalf, 
                                                amo3_TTS, 
                                                tau_amo3, epsilon_amo3)$mean
save(amo3_comp_TTS, file = paste0(current_path, "amo3_comp_TTS.Rda"))
save(amo3_compholdout_TTS, file = paste0(current_path, 
                                         "amo3_compholdout_TTS.Rda"))
rm(amo3_TTS)
gc()

########################################################################

round(data_amo3_mean_firsthalf - tau_amo3,2)

load(paste0(current_path, "amo3_comp_UNIFORM.Rda"))
load(paste0(current_path, "amo3_comp_APT.Rda"))
load(paste0(current_path, "amo3_comp_LR.Rda"))
load(paste0(current_path, "amo3_comp_AugUCB.Rda"))
load(paste0(current_path, "amo3_comp_EVT.Rda"))
load(paste0(current_path, "amo3_comp_BUCB.Rda"))

plot(c(0,10080), c(0, -3), type = "n")
lines(log(amo3_comp_UNIFORM), col = "black")
lines(log(amo3_comp_APT), col = "blue")
lines(log(amo3_comp_AugUCB), col = "green")
lines(log(amo3_comp_EVT), col = "darkgreen")
lines(log(amo3_comp_LR), col = "red")
abline(h = log(0.1), lty = 2)

round(data_amo3_mean_secondhalf - tau_amo3,2)
load(paste0(current_path, "amo3_compholdout_BUCB.Rda"))
load(paste0(current_path, "amo3_compholdout_APT.Rda"))
load(paste0(current_path, "amo3_compholdout_UNIFORM.Rda"))
plot(c(0,10080), c(0, -4), type = "n")
lines(log(amo3_compholdout_UNIFORM), col = "black")
lines(log(amo3_compholdout_APT), col = "blue")
lines(log(amo3_compholdout_AugUCB), col = "red")
lines(log(amo3_compholdout_BUCB), col = "green")
lines(log(amo3_compholdout_KLUCB), col = "darkgreen")
lines(log(amo3_compholdout_TTS), col = "violet")
abline(h = log(0.1), lty = 2)
