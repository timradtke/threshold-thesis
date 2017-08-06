########################################################################
# load functions

current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/AugUCB/"
source(paste0(current_path, "run_from_data.R"))
source(paste0(current_path, "run_parallelized_from_data.R"))
current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/synthetic_experiments/exponential/"

########################################################################
# Create the data from exponential distributions

set.seed(834795)
tau_pois1 <- 1
epsilon_pois1 <- 0
mu_pois1 <- rexp(20, 1/5)

plot(c(0,8), c(0,2), type = "n")
for(i in order(mu_pois1)) {
  lines(seq(0,8,0.001), dpois(seq(0,8,0.001), mu_pois1[i]), col = rainbow(length(mu_pois1))[i])
}

data_list_pois1 <- list()
for(j in 1:5000) {
  curr_data <- data.frame(rep(NA, times = 10000))
  for(i in 1:length(mu_pois1)) {
    curr_data[[i]] <- as.numeric(rpois(10000, 1/mu_pois1[i]))
  }
  names(curr_data) <- paste0("V", rep(1:length(mu_pois1)))
  data_list_pois1[[j]] <- curr_data
}
gc()

########################################################################

plot(c(0,10000), c(0, -8), type = "n")
abline(h=log(0.01), lty = 2)
lines(log(pois1_comp_APT), col = "black")
lines(log(pois1_comp_LR), col = "blue")
lines(log(pois1_comp_EVT), col = "red")

########################################################################
# Gaussian Likelihood Ratio
system.time(pois1_LR <- para_bandit_sim_LR_poisson(data = data_list_pois1, 
                                                   rounds = 10000, 
                                                   tau = tau_pois1, 
                                                   epsilon = epsilon_pois1))

save(pois1_LR, file = paste0(current_path, "pois1_LR.Rda"))
pois1_comp_LR <- compare_to_ground_truth(mu_pois1, pois1_LR, 
                                        tau_pois1, epsilon_pois1)$mean
save(pois1_comp_LR, file = paste0(current_path, "pois1_comp_LR.Rda"))
gc()

########################################################################

system.time(pois1_APT <- para_bandit_sim_APT(data = data_list_pois1, 
                                            rounds = 10000, 
                                            tau = tau_pois1, 
                                            epsilon = epsilon_pois1))

save(pois1_APT, file = paste0(current_path, "pois1_APT.Rda"))
pois1_comp_APT <- compare_to_ground_truth(mu_pois1, pois1_APT, 
                                         tau_pois1, epsilon_pois1)$mean
save(pois1_comp_APT, file = paste0(current_path, "pois1_comp_APT.Rda"))
gc()

########################################################################

system.time(pois1_EVT <- para_bandit_sim_EVT(data = data_list_pois1, 
                                            rounds = 10000, 
                                            tau = tau_pois1, 
                                            epsilon = epsilon_pois1))

save(pois1_EVT, file = paste0(current_path, "pois1_EVT.Rda"))
pois1_comp_EVT <- compare_to_ground_truth(mu_pois1, pois1_EVT, 
                                          tau_pois1, epsilon_pois1)$mean
save(pois1_comp_EVT, file = paste0(current_path, "pois1_comp_EVT.Rda"))
gc()

########################################################################
# Uniform Sampling

system.time(pois1_UNIFORM <- para_bandit_sim_uniform(data = data_list_pois1, 
                                                    rounds = 10000))
save(pois1_UNIFORM, file = paste0(current_path, "pois1_UNIFORM.Rda"))
pois1_comp_UNIFORM <- compare_to_ground_truth(mean_pois1, 
                                             pois1_UNIFORM, 
                                             tau_pois1, 
                                             epsilon_pois1)$mean
save(pois1_comp_UNIFORM, file = paste0(current_path, "pois1_comp_UNIFORM.Rda"))

########################################################################
# Augmented-UCB (Mukherjee et al., 2017)

pois1_AugUCB <- para_bandit_sim_AugUCB(data = data_list_pois1, 
                                      rounds = 10000, 
                                      tau = tau_pois1)
save(pois1_AugUCB, file = paste0(current_path, "pois1_AugUCB.Rda"))
pois1_comp_AugUCB <- compare_to_ground_truth(mean_pois1, pois1_AugUCB, 
                                            tau_pois1, epsilon_pois1)$mean
save(pois1_comp_AugUCB, file = paste0(current_path, "pois1_comp_AugUCB.Rda"))

