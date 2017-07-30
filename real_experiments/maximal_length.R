sessions_threshold <- 0
set.seed(512)
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

data_amo1 <- pv_list
data_amo1_wide <- pv_products_wide
data_amo1_own_means <- pv_own_means
data_amo1_next_means <- pv_next_means
data_amo1_mean_secondhalf <- colMeans(data_amo1_wide[12581:20160,])
data_amo1_mean_middlehalf <- colMeans(data_amo1_wide[7581:15161,])
data_amo1_mean_firsthalf <- colMeans(data_amo1_wide[1:12580,])
save(data_amo1, file = paste0(current_path, "data_amo1.Rda"))
save(data_amo1_wide, file = paste0(current_path, "data_amo1_wide.Rda"))
save(data_amo1_own_means, file = paste0(current_path, "data_amo1_own_means.Rda"))
save(data_amo1_next_means, file = paste0(current_path, "data_amo1_next_means.Rda"))
save(data_amo1_mean_firsthalf, file = paste0(current_path, "data_amo1_mean_firsthalf.Rda"))
save(data_amo1_mean_secondhalf, file = paste0(current_path, "data_amo1_mean_secondhalf.Rda"))
save(data_amo1_mean_middlehalf, file = paste0(current_path, "data_amo1_mean_middlehalf.Rda"))
rm(pv_list, pv_own_means, pv_next_means, pv_products_wide)
gc()

load(paste0(current_path, "data_amo1_mean_firsthalf"))
load(paste0(current_path, "data_amo1_mean_secondhalf"))
load(paste0(current_path, "data_amo1.Rda"))
load(paste0(current_path, "data_amo1_wide.Rda"))
tau_amo1 <- 3/60
epsilon_amo1 <- 1/60

########################################################################
# Standard APT Algorithm
system.time(amo1_APT <- para_bandit_sim_APT(data = data_amo1, rounds = 7580, 
                                            tau = tau_amo1, epsilon = epsilon_amo1))
save(amo1_APT, file = paste0(current_path, "amo1_APT.Rda"))
load(paste0(current_path, "amo1_APT.Rda"))
amo1_compnext_APT <- compare_to_cv_data(data_amo1_next_means, amo1_APT, 
                                        tau_amo1, epsilon_amo1)$mean
amo1_compown_APT <- compare_to_cv_data(data_amo1_own_means, amo1_APT, 
                                        tau_amo1, epsilon_amo1)$mean
amo1_comp_APT <- compare_to_ground_truth(data_amo1_mean_firsthalf, amo1_APT, 
                                         tau_amo1, epsilon_amo1)$mean
amo1_compholdout_APT <- compare_to_ground_truth(data_amo1_mean_secondhalf, amo1_APT, 
                                            tau_amo1, epsilon_amo1)$mean

save(amo1_compnext_APT, file = paste0(current_path, "amo1_compnext_APT.Rda"))
save(amo1_compown_APT, file = paste0(current_path, "amo1_compown_APT.Rda"))
save(amo1_comp_APT, file = paste0(current_path, "amo1_comp_APT.Rda"))
save(amo1_compholdout_APT, file = paste0(current_path, "amo1_compholdout_APT.Rda"))
rm(amo1_APT)
gc()


load(paste0(current_path, "amo1_compown_APT.Rda"))
load(paste0(current_path, "amo1_compnext_APT.Rda"))
plot(c(0,8000), c(0, -4), type = "n")
#lines(log(amo1_compnext_APT), col = "red")
#lines(log(amo1_compown_APT), col = "blue")
lines(log(amo1_comp_APT), col = "blue")
lines(log(amo1_compholdout_APT), col = "darkblue")
lines(log(amo1_comp_BUCB), col = "green")
lines(log(amo1_compholdout_BUCB), col = "darkgreen")

plot(c(0,8000), c(0, 1), type = "n")
lines((amo1_comp_APT), col = "blue")
lines((amo1_compholdout_APT), col = "darkblue")
lines((amo1_comp_BUCB), col = "green")
lines((amo1_compholdout_BUCB), col = "darkgreen")
abline(h = 0.1)

########################################################################
# Bayes-UCB

amo1_BUCB <- para_bandit_sim_bucb(data = data_amo1, rounds = 7580, 
                                   rate = "inverse_horizon",
                                   tau = tau_amo1, epsilon = epsilon_amo1, 
                                   alpha = tau_amo1, beta = 1-tau_amo1)
save(amo1_BUCB, 
     file = paste0(current_path, "amo1_BUCB.Rda"))
amo1_compnext_BUCB <- compare_to_cv_data(data_amo1_next_means, amo1_BUCB, 
                                        tau_amo1, epsilon_amo1)$mean
amo1_compown_BUCB <- compare_to_cv_data(data_amo1_own_means, amo1_BUCB, 
                                       tau_amo1, epsilon_amo1)$mean
amo1_comp_BUCB <- compare_to_ground_truth(data_amo1_mean_firsthalf, amo1_BUCB, 
                                         tau_amo1, epsilon_amo1)$mean
amo1_compholdout_BUCB <- compare_to_ground_truth(data_amo1_mean_secondhalf, 
                                                 amo1_BUCB, 
                                                 tau_amo1, epsilon_amo1)$mean
save(amo1_compnext_BUCB, file = paste0(current_path, "amo1_compnext_BUCB.Rda"))
save(amo1_compown_BUCB, file = paste0(current_path, "amo1_compown_BUCB.Rda"))
save(amo1_comp_BUCB, file = paste0(current_path, "amo1_comp_BUCB.Rda"))
save(amo1_compholdout_BUCB, file = paste0(current_path, "amo1_compholdout_BUCB.Rda"))
rm(amo1_BUCB)
gc()


########################################################################
# LR Algorithm
system.time(amo1_LR <- para_bandit_sim_LR(data = data_amo1, rounds = 7580, 
                                            tau = tau_amo1, epsilon = epsilon_amo1))
save(amo1_LR, file = paste0(current_path, "amo1_LR.Rda"))
#load(paste0(current_path, "amo1_LR.Rda"))
amo1_comp_LR <- compare_to_ground_truth(data_amo1_mean_firsthalf, amo1_LR, 
                                         tau_amo1, epsilon_amo1)$mean
amo1_compholdout_LR <- compare_to_ground_truth(data_amo1_mean_secondhalf, amo1_LR, 
                                                tau_amo1, epsilon_amo1)$mean

save(amo1_comp_LR, file = paste0(current_path, "amo1_comp_LR.Rda"))
save(amo1_compholdout_LR, file = paste0(current_path, "amo1_compholdout_LR.Rda"))
rm(amo1_LR)
gc()

########################################################################
# Standard Uniform
system.time(amo1_UNIFORM <- para_bandit_sim_uniform(data = data_amo1, 
                                                    rounds = 7580))
save(amo1_UNIFORM, file = paste0(current_path, "amo1_UNIFORM.Rda"))
load(file = paste0(current_path, "amo1_UNIFORM.Rda"))

amo1_comp_UNIFORM <- compare_to_ground_truth(data_amo1_mean_firsthalf, 
                                             amo1_UNIFORM, 
                                             tau_amo1, epsilon_amo1)$mean
amo1_compholdout_UNIFORM <- compare_to_ground_truth(data_amo1_mean_secondhalf, 
                                                    amo1_UNIFORM, 
                                                    tau_amo1, epsilon_amo1)$mean
save(amo1_comp_UNIFORM, file = paste0(current_path, "amo1_comp_UNIFORM.Rda"))
save(amo1_compholdout_UNIFORM, file = paste0(current_path, 
                                             "amo1_compholdout_UNIFORM.Rda"))
rm(amo1_UNIFORM)
gc()

########################################################################
# Standard AugUCB
system.time(amo1_AugUCB <- para_bandit_sim_AugUCB(data = data_amo1, 
                                                  rounds = 7580, 
                                                  tau = tau_amo1))
save(amo1_AugUCB, file = paste0(current_path, "amo1_AugUCB.Rda"))
load(file = paste0(current_path, "amo1_AugUCB.Rda"))

amo1_comp_AugUCB <- compare_to_ground_truth(data_amo1_mean_firsthalf, 
                                             amo1_AugUCB, 
                                             tau_amo1, epsilon_amo1)$mean
amo1_compholdout_AugUCB <- compare_to_ground_truth(data_amo1_mean_secondhalf, 
                                                   amo1_AugUCB, 
                                                   tau_amo1, epsilon_amo1)$mean
save(amo1_comp_AugUCB, file = paste0(current_path, "amo1_comp_AugUCB.Rda"))
save(amo1_compholdout_AugUCB, file = paste0(current_path, 
                                             "amo1_compholdout_AugUCB.Rda"))
rm(amo1_AugUCB)
gc()

########################################################################
# Standard KL-UCB
Hholdout_amo1 <- get_complexity(data_amo1_mean_secondhalf, tau_amo1, epsilon_amo1)
H_amo1 <- get_complexity(data_amo1_mean_firsthalf, tau_amo1, epsilon_amo1)

system.time(amo1_KLUCB <- para_bandit_sim_KLUCB(data = data_amo1, 
                                                       rounds = 7580, 
                                                       tau = tau_amo1, 
                                                       epsilon = epsilon_amo1,
                                                       horizon = 5.5*H_amo1,
                                                       H = H_amo1))
save(amo1_KLUCB, file = paste0(current_path, "amo1_KLUCB.Rda"))
load(file = paste0(current_path, "amo1_KLUCB.Rda"))

amo1_comp_KLUCB <- compare_to_ground_truth(data_amo1_mean_firsthalf, 
                                            amo1_KLUCB, 
                                            tau_amo1, epsilon_amo1)$mean
amo1_compholdout_KLUCB <- compare_to_ground_truth(data_amo1_mean_secondhalf, 
                                                   amo1_KLUCB, 
                                                   tau_amo1, epsilon_amo1)$mean
save(amo1_comp_KLUCB, file = paste0(current_path, "amo1_comp_KLUCB.Rda"))
save(amo1_compholdout_KLUCB, file = paste0(current_path, 
                                            "amo1_compholdout_KLUCB.Rda"))
rm(amo1_KLUCB)
gc()

########################################################################

amo1_TTS <- para_bandit_sim_TTS(data = data_amo1, rounds = 7580,
                                tau = tau_amo1, epsilon = epsilon_amo1,
                                alpha = tau_amo1, beta = 1 - tau_amo1)
save(amo1_TTS, file = paste0(current_path, "amo1_TTS.Rda"))

amo1_comp_TTS <- compare_to_ground_truth(data_amo1_mean_firsthalf, 
                                           amo1_TTS, 
                                           tau_amo1, epsilon_amo1)$mean
amo1_compholdout_TTS <- compare_to_ground_truth(data_amo1_mean_secondhalf, 
                                                  amo1_TTS, 
                                                  tau_amo1, epsilon_amo1)$mean
save(amo1_comp_TTS, file = paste0(current_path, "amo1_comp_TTS.Rda"))
save(amo1_compholdout_TTS, file = paste0(current_path, 
                                           "amo1_compholdout_TTS.Rda"))


########################################################################

data_amo1_mean_firsthalf > 3/60
data_amo1_mean_firsthalf < 1/60
load(paste0(current_path, "amo1_comp_APT.Rda"))
load(paste0(current_path, "amo1_comp_BUCB.Rda"))
load(paste0(current_path, "amo1_comp_KLUCB.Rda"))
load(paste0(current_path, "amo1_comp_AugUCB.Rda"))
load(paste0(current_path, "amo1_comp_UNIFORM.Rda"))
load(paste0(current_path, "amo1_comp_TTS.Rda"))
load(paste0(current_path, "amo1_comp_LR.Rda"))
plot(c(0,8000), c(0, -8), type = "n")
lines(log(amo1_comp_UNIFORM), col = "black")
lines(log(amo1_comp_APT), col = "blue")
lines(log(amo1_comp_AugUCB), col = "red")
lines(log(amo1_comp_BUCB), col = "green")
lines(log(amo1_comp_KLUCB), col = "darkgreen")
lines(log(amo1_comp_TTS), col = "violet")
lines(log(amo1_comp_LR), col = "darkblue")
abline(h = log(0.1), lty = 2)

data_amo1_mean_secondhalf > 3/60
data_amo1_mean_secondhalf < 1/60
load(paste0(current_path, "amo1_compholdout_BUCB.Rda"))
load(paste0(current_path, "amo1_compholdout_APT.Rda"))
load(paste0(current_path, "amo1_compholdout_UNIFORM.Rda"))
load(paste0(current_path, "amo1_compholdout_KLUCB.Rda"))
load(paste0(current_path, "amo1_compholdout_AugUCB.Rda"))
load(paste0(current_path, "amo1_compholdout_TTS.Rda"))
load(paste0(current_path, "amo1_compholdout_LR.Rda"))
plot(c(0,8000), c(0, -4), type = "n")
lines(log(amo1_compholdout_UNIFORM), col = "black")
lines(log(amo1_compholdout_APT), col = "blue")
lines(log(amo1_compholdout_AugUCB), col = "red")
lines(log(amo1_compholdout_BUCB), col = "green")
lines(log(amo1_compholdout_KLUCB), col = "darkgreen")
lines(log(amo1_compholdout_TTS), col = "violet")
lines(log(amo1_compholdout_LR), col = "darkblue")
abline(h = log(0.1), lty = 2)
