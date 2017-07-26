sessions_threshold <- 0
set.seed(1068)
selected_products <- sample(1:196, 10)
niter <- 5000
nrounds <- 7580 # These two values give the maximum size data 
nahead <- 7580  # possible given the data set size we have
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

data_amo2 <- pv_list
data_amo2_wide <- pv_products_wide
data_amo2_own_means <- pv_own_means
data_amo2_next_means <- pv_next_means
data_amo2_mean_secondhalf <- colMeans(data_amo2_wide[12581:20160,])
data_amo2_mean_middlehalf <- colMeans(data_amo2_wide[7581:15161,])
data_amo2_mean_firsthalf <- colMeans(data_amo2_wide[1:12580,])
save(data_amo2, file = paste0(current_path, "data_amo2.Rda"))
save(data_amo2_wide, file = paste0(current_path, "data_amo2_wide.Rda"))
save(data_amo2_own_means, file = paste0(current_path, "data_amo2_own_means.Rda"))
save(data_amo2_next_means, file = paste0(current_path, "data_amo2_next_means.Rda"))
save(data_amo2_mean_firsthalf, file = paste0(current_path, "data_amo2_mean_firsthalf.Rda"))
save(data_amo2_mean_secondhalf, file = paste0(current_path, "data_amo2_mean_secondhalf.Rda"))
rm(pv_list, pv_own_means, pv_next_means, pv_products_wide)
gc()

load(paste0(current_path, "data_amo2_mean_firsthalf.Rda"))
load(paste0(current_path, "data_amo2_mean_secondhalf.Rda"))
load(paste0(current_path, "data_amo2.Rda"))
tau_amo2 <- 4/60
epsilon_amo2 <- 1/60

########################################################################
# Bayes-UCB

amo2_BUCB <- para_bandit_sim_bucb(data = data_amo2, rounds = 7580, 
                                  rate = "inverse_horizon",
                                  tau = tau_amo2, epsilon = epsilon_amo2, 
                                  alpha = tau_amo2, beta = 1-tau_amo2)
save(amo2_BUCB, 
     file = paste0(current_path, "amo2_BUCB.Rda"))
amo2_compnext_BUCB <- compare_to_cv_data(data_amo2_next_means, amo2_BUCB, 
                                         tau_amo2, epsilon_amo2)$mean
amo2_compown_BUCB <- compare_to_cv_data(data_amo2_own_means, amo2_BUCB, 
                                        tau_amo2, epsilon_amo2)$mean
amo2_comp_BUCB <- compare_to_ground_truth(data_amo2_mean_firsthalf, amo2_BUCB, 
                                          tau_amo2, epsilon_amo2)$mean
amo2_compholdout_BUCB <- compare_to_ground_truth(data_amo2_mean_secondhalf, 
                                                 amo2_BUCB, 
                                                 tau_amo2, epsilon_amo2)$mean
save(amo2_compnext_BUCB, file = paste0(current_path, "amo2_compnext_BUCB.Rda"))
save(amo2_compown_BUCB, file = paste0(current_path, "amo2_compown_BUCB.Rda"))
save(amo2_comp_BUCB, file = paste0(current_path, "amo2_comp_BUCB.Rda"))
save(amo2_compholdout_BUCB, file = paste0(current_path, "amo2_compholdout_BUCB.Rda"))
rm(amo2_BUCB)
gc()


########################################################################
# Standard Uniform
system.time(amo2_UNIFORM <- para_bandit_sim_uniform(data = data_amo2, 
                                                    rounds = 7580))
save(amo2_UNIFORM, file = paste0(current_path, "amo2_UNIFORM.Rda"))
load(file = paste0(current_path, "amo2_UNIFORM.Rda"))

amo2_comp_UNIFORM <- compare_to_ground_truth(data_amo2_mean_firsthalf, 
                                             amo2_UNIFORM, 
                                             tau_amo2, epsilon_amo2)$mean
amo2_compholdout_UNIFORM <- compare_to_ground_truth(data_amo2_mean_secondhalf, 
                                                    amo2_UNIFORM, 
                                                    tau_amo2, epsilon_amo2)$mean
save(amo2_comp_UNIFORM, file = paste0(current_path, "amo2_comp_UNIFORM.Rda"))
save(amo2_compholdout_UNIFORM, file = paste0(current_path, 
                                             "amo2_compholdout_UNIFORM.Rda"))
rm(amo2_UNIFORM)
gc()

########################################################################
# Standard APT Algorithm
system.time(amo2_APT <- para_bandit_sim_APT(data = data_amo2, rounds = 7580, 
                                            tau = tau_amo2, epsilon = epsilon_amo2))
save(amo2_APT, file = paste0(current_path, "amo2_APT.Rda"))
#load(paste0(current_path, "amo2_APT.Rda"))
amo2_comp_APT <- compare_to_ground_truth(data_amo2_mean_firsthalf, amo2_APT, 
                                         tau_amo2, epsilon_amo2)$mean
amo2_compholdout_APT <- compare_to_ground_truth(data_amo2_mean_secondhalf, amo2_APT, 
                                                tau_amo2, epsilon_amo2)$mean

save(amo2_comp_APT, file = paste0(current_path, "amo2_comp_APT.Rda"))
save(amo2_compholdout_APT, file = paste0(current_path, "amo2_compholdout_APT.Rda"))
rm(amo2_APT)
gc()

########################################################################
# Standard KL-UCB
Hholdout_amo2 <- get_complexity(data_amo2_mean_secondhalf, tau_amo2, epsilon_amo2)
H_amo2 <- get_complexity(data_amo2_mean_firsthalf, tau_amo2, epsilon_amo2)

system.time(amo2_KLUCB <- para_bandit_sim_KLUCB(data = data_amo2, 
                                                rounds = 7580, 
                                                tau = tau_amo2, 
                                                epsilon = epsilon_amo2,
                                                horizon = 5.5*H_amo2,
                                                H = H_amo2))
save(amo2_KLUCB, file = paste0(current_path, "amo2_KLUCB.Rda"))
#load(file = paste0(current_path, "amo2_KLUCB.Rda"))

amo2_comp_KLUCB <- compare_to_ground_truth(data_amo2_mean_firsthalf, 
                                           amo2_KLUCB, 
                                           tau_amo2, epsilon_amo2)$mean
amo2_compholdout_KLUCB <- compare_to_ground_truth(data_amo2_mean_secondhalf, 
                                                  amo2_KLUCB, 
                                                  tau_amo2, epsilon_amo2)$mean
save(amo2_comp_KLUCB, file = paste0(current_path, "amo2_comp_KLUCB.Rda"))
save(amo2_compholdout_KLUCB, file = paste0(current_path, 
                                           "amo2_compholdout_KLUCB.Rda"))
rm(amo2_KLUCB)
gc()

########################################################################

amo2_TTS <- para_bandit_sim_TTS(data = data_amo2, rounds = 7580,
                                tau = tau_amo2, epsilon = epsilon_amo2,
                                alpha = tau_amo2, beta = 1 - tau_amo2)
save(amo2_TTS, file = paste0(current_path, "amo2_TTS.Rda"))

amo2_comp_TTS <- compare_to_ground_truth(data_amo2_mean_firsthalf, 
                                         amo2_TTS, 
                                         tau_amo2, epsilon_amo2)$mean
amo2_compholdout_TTS <- compare_to_ground_truth(data_amo2_mean_secondhalf, 
                                                amo2_TTS, 
                                                tau_amo2, epsilon_amo2)$mean
save(amo2_comp_TTS, file = paste0(current_path, "amo2_comp_TTS.Rda"))
save(amo2_compholdout_TTS, file = paste0(current_path, 
                                         "amo2_compholdout_TTS.Rda"))
rm(amo2_TTS)
gc()
########################################################################
# Standard AugUCB
system.time(amo2_AugUCB <- para_bandit_sim_AugUCB(data = data_amo2, 
                                                  rounds = 7580, 
                                                  tau = tau_amo2))
save(amo2_AugUCB, file = paste0(current_path, "amo2_AugUCB.Rda"))
load(file = paste0(current_path, "amo2_AugUCB.Rda"))

amo2_comp_AugUCB <- compare_to_ground_truth(data_amo2_mean_firsthalf, 
                                            amo2_AugUCB, 
                                            tau_amo2, epsilon_amo2)$mean
amo2_compholdout_AugUCB <- compare_to_ground_truth(data_amo2_mean_secondhalf, 
                                                   amo2_AugUCB, 
                                                   tau_amo2, epsilon_amo2)$mean
save(amo2_comp_AugUCB, file = paste0(current_path, "amo2_comp_AugUCB.Rda"))
save(amo2_compholdout_AugUCB, file = paste0(current_path, 
                                            "amo2_compholdout_AugUCB.Rda"))
rm(amo2_AugUCB)
gc()

########################################################################

data_amo2_mean_firsthalf > 3/60
data_amo2_mean_firsthalf < 1/60
load(paste0(current_path, "amo2_comp_APT.Rda"))
load(paste0(current_path, "amo2_comp_BUCB.Rda"))
plot(c(0,8000), c(0, -8), type = "n")
lines(log(amo2_comp_UNIFORM), col = "black")
lines(log(amo2_comp_APT), col = "blue")
lines(log(amo2_comp_AugUCB), col = "red")
lines(log(amo2_comp_BUCB), col = "green")
lines(log(amo2_comp_KLUCB), col = "darkgreen")
lines(log(amo2_comp_TTS), col = "violet")
abline(h = log(0.1), lty = 2)

data_amo2_mean_secondhalf > 3/60
data_amo2_mean_secondhalf < 1/60
load(paste0(current_path, "amo2_compholdout_BUCB.Rda"))
load(paste0(current_path, "amo2_compholdout_APT.Rda"))
plot(c(0,8000), c(0, -6), type = "n")
lines(log(amo2_compholdout_UNIFORM), col = "black")
lines(log(amo2_compholdout_APT), col = "blue")
lines(log(amo2_compholdout_AugUCB), col = "red")
lines(log(amo2_compholdout_BUCB), col = "green")
lines(log(amo2_compholdout_KLUCB), col = "darkgreen")
lines(log(amo2_compholdout_TTS), col = "violet")
abline(h = log(0.1), lty = 2)
