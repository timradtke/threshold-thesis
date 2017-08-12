# What happens if we use synthetic data with the same mean as in the real data?

current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/AugUCB/"
source(paste0(current_path, "run_from_data.R"))
source(paste0(current_path, "run_parallelized_from_data.R"))
current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/real_experiments/"

load(paste0(current_path, "data_amo4_mean_firsthalf.Rda"))
tau_amo4sim <- 5/60
epsilon_amo4sim <- 0

data_amo4sim <- list()
set.seed(458693)
for(j in 1:5000) {
  curr_data <- data.frame(rep(NA, times = 10080))
  for(i in 1:length(data_amo4_mean_firsthalf)) {
    curr_data[[i]] <- as.numeric(purrr::rbernoulli(10080, p  = data_amo4_mean_firsthalf[i]))
  }
  names(curr_data) <- paste0("V", rep(1:length(data_amo4_mean_firsthalf)))
  data_amo4sim[[j]] <- curr_data
}
gc()




########################################################################
# Standard Likelihood Ratio
system.time(amo4sim_LR <- para_bandit_sim_LR(data = data_amo4sim, 
                                          rounds = 10080, 
                                          tau = tau_amo4sim, 
                                          epsilon = epsilon_amo4sim))

save(amo4sim_LR, file = paste0(current_path, "amo4sim_LR.Rda"))
#load(file = paste0(current_path, "amo4sim_LR.Rda"))
amo4sim_comp_LR <- compare_to_ground_truth(data_amo4_mean_firsthalf, 
                                           amo4sim_LR, 
                                           tau_amo4sim, epsilon_amo4sim)$mean
save(amo4sim_comp_LR, file = paste0(current_path, "amo4sim_comp_LR.Rda"))
rm(amo4sim_LR, amo4sim_comp_LR)
gc()

########################################################################
# Standard Uniform
system.time(amo4sim_UNIFORM <- para_bandit_sim_uniform(data = data_amo4sim, 
                                                       rounds = 10080))
save(amo4sim_UNIFORM, file = paste0(current_path, "amo4sim_UNIFORM.Rda"))
#load(file = paste0(current_path, "amo4sim_UNIFORM.Rda"))

amo4sim_comp_UNIFORM <- compare_to_ground_truth(data_amo4_mean_firsthalf, 
                                                amo4sim_UNIFORM, 
                                                tau_amo4sim, epsilon_amo4sim)$mean
save(amo4sim_comp_UNIFORM, file = paste0(current_path, "amo4sim_comp_UNIFORM.Rda"))
rm(amo4sim_UNIFORM)
gc()

########################################################################
# Empirical Variance Guided Algorithm (Zhong et al., 2017)

system.time(amo4sim_EVT <- para_bandit_sim_EVT(data = data_amo4sim, 
                                            rounds = 10080, 
                                            tau = tau_amo4sim, 
                                            epsilon = epsilon_amo4sim))
save(amo4sim_EVT, file = paste0(current_path, "amo4sim_EVT.Rda"))
amo4sim_comp_EVT <- compare_to_ground_truth(data_amo4_mean_firsthalf, 
                                            amo4sim_EVT, 
                                            tau_amo4sim, epsilon_amo4sim)$mean
save(amo4sim_comp_EVT, file = paste0(current_path, "amo4sim_comp_EVT.Rda"))
rm(amo4sim_EVT, amo4sim_comp_EVT)
gc()

########################################################################
# Bayes-UCB

amo4sim_BUCB <- para_bandit_sim_bucb(data = data_amo4sim, rounds = 10080, 
                                     rate = "inverse_horizon",
                                     tau = tau_amo4sim, epsilon = epsilon_amo4sim, 
                                     alpha = tau_amo4sim, beta = 1-tau_amo4sim)
save(amo4sim_BUCB, 
     file = paste0(current_path, "amo4sim_BUCB.Rda"))
amo4sim_comp_BUCB <- compare_to_ground_truth(data_amo4_mean_firsthalf, amo4sim_BUCB, 
                                             tau_amo4sim, epsilon_amo4sim)$mean
save(amo4sim_comp_BUCB, file = paste0(current_path, "amo4sim_comp_BUCB.Rda"))
rm(amo4sim_BUCB)
gc()

########################################################################
# Standard APT Algorithm
system.time(amo4sim_APT <- para_bandit_sim_APT(data = data_amo4sim, rounds = 10080, 
                                               tau = tau_amo4sim, epsilon = epsilon_amo4sim))
save(amo4sim_APT, file = paste0(current_path, "amo4sim_APT.Rda"))
#load(paste0(current_path, "amo4sim_APT.Rda"))
amo4sim_comp_APT <- compare_to_ground_truth(data_amo4_mean_firsthalf, amo4sim_APT, 
                                            tau_amo4sim, epsilon_amo4sim)$mean

save(amo4sim_comp_APT, file = paste0(current_path, "amo4sim_comp_APT.Rda"))
rm(amo4sim_APT)
gc()

########################################################################
# Standard AugUCB
system.time(amo4sim_AugUCB <- para_bandit_sim_AugUCB(data = data_amo4sim, 
                                                  rounds = 10080, 
                                                  tau = tau_amo4sim))
save(amo4sim_AugUCB, file = paste0(current_path, "amo4sim_AugUCB.Rda"))
#load(file = paste0(current_path, "amo4sim_AugUCB.Rda"))

amo4sim_comp_AugUCB <- compare_to_ground_truth(data_amo4_mean_firsthalf, 
                                            amo4sim_AugUCB, 
                                            tau_amo4sim, epsilon_amo4sim)$mean
save(amo4sim_comp_AugUCB, file = paste0(current_path, "amo4sim_comp_AugUCB.Rda"))
rm(amo4sim_AugUCB, amo4sim_comp_AugUCB)
gc()




round(data_amo4_mean_firsthalf - tau_amo4,2)

load(paste0(current_path, "amo4_comp_UNIFORM.Rda"))
load(paste0(current_path, "amo4_comp_APT.Rda"))
load(paste0(current_path, "amo4sim_comp_LR.Rda"))
load(paste0(current_path, "amo4_comp_AugUCB.Rda"))
load(paste0(current_path, "amo4_comp_EVT.Rda"))
load(paste0(current_path, "amo4_comp_BUCB.Rda"))

plot(c(0,10080), c(0, -5), type = "n")
lines(log(amo4sim_comp_UNIFORM), col = "black")
lines(log(amo4_comp_APT), col = "blue")
lines(log(amo4_comp_AugUCB), col = "green")
lines(log(amo4_comp_EVT), col = "darkgreen")
lines(log(amo4sim_comp_LR), col = "red")
lines(log(amo4_comp_BUCB), col = "violet")
abline(h = log(0.1), lty = 2)
