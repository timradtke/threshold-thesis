########################################################################
# load functions

current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/AugUCB/"
source(paste0(current_path, "run_from_data.R"))
source(paste0(current_path, "run_parallelized_from_data.R"))
current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/synthetic_experiments/poisson/"

########################################################################
# Create the data from exponential distributions

set.seed(3795969)
tau_pois2 <- 0.4
epsilon_pois2 <- 0
mu_pois2 <- rexp(20, 1/0.7)

#plot(c(0,8), c(0,2), type = "n")
#for(i in order(mu_pois2)) {
#  lines(seq(0,8,0.001), dpois(seq(0,8,0.001), mu_pois2[i]), col = rainbow(length(mu_pois2))[i])
#}

data_list_pois2 <- list()
for(j in 1:5000) {
  curr_data <- data.frame(rep(NA, times = 10000))
  for(i in 1:length(mu_pois2)) {
    curr_data[[i]] <- as.numeric(rpois(10000, mu_pois2[i]))
  }
  names(curr_data) <- paste0("V", rep(1:length(mu_pois2)))
  data_list_pois2[[j]] <- curr_data
}
gc()

########################################################################

load(paste0(current_path, "pois2_comp_LR.Rda"))
load(paste0(current_path, "pois2_comp_APT.Rda"))
load(paste0(current_path, "pois2_comp_EVT.Rda"))

plot(c(0,10000), c(0, -7), type = "n")
abline(h=log(0.1), lty = 2)
abline(h=log(0.01), lty = 2)
lines(log(pois2_comp_UNIFORM), col = "black")
lines(log(pois2_comp_APT), col = "darkgreen")
lines(log(pois2_comp_LR), col = "blue")
lines(log(pois2_comp_EVT), col = "red")

########################################################################
# Poisson Likelihood Ratio

system.time(pois2_LR <- para_bandit_sim_LR_poisson(data = data_list_pois2, 
                                                   rounds = 10000, 
                                                   tau = tau_pois2, 
                                                   epsilon = epsilon_pois2))

save(pois2_LR, file = paste0(current_path, "pois2_LR.Rda"))
pois2_comp_LR <- compare_to_ground_truth(mu_pois2, pois2_LR, 
                                         tau_pois2, epsilon_pois2)$mean
save(pois2_comp_LR, file = paste0(current_path, "pois2_comp_LR.Rda"))
rm(pois2_LR)
gc()

########################################################################
# Anytime Parameter-free Thresholding Algorithm (Locatelli et al., 2016)

system.time(pois2_APT <- para_bandit_sim_APT(data = data_list_pois2, 
                                             rounds = 10000, 
                                             tau = tau_pois2, 
                                             epsilon = epsilon_pois2))

save(pois2_APT, file = paste0(current_path, "pois2_APT.Rda"))
pois2_comp_APT <- compare_to_ground_truth(mu_pois2, pois2_APT, 
                                          tau_pois2, epsilon_pois2)$mean
save(pois2_comp_APT, file = paste0(current_path, "pois2_comp_APT.Rda"))
rm(pois2_APT)
gc()

########################################################################
# Empirical Variance Guided Algorithm (Zhong et al., 2017)

system.time(pois2_EVT <- para_bandit_sim_EVT(data = data_list_pois2, 
                                             rounds = 10000, 
                                             tau = tau_pois2, 
                                             epsilon = epsilon_pois2))

save(pois2_EVT, file = paste0(current_path, "pois2_EVT.Rda"))
pois2_comp_EVT <- compare_to_ground_truth(mu_pois2, pois2_EVT, 
                                          tau_pois2, epsilon_pois2)$mean
save(pois2_comp_EVT, file = paste0(current_path, "pois2_comp_EVT.Rda"))
gc()

########################################################################
# Uniform Sampling

system.time(pois2_UNIFORM <- para_bandit_sim_uniform(data = data_list_pois2, 
                                                     rounds = 10000))
save(pois2_UNIFORM, file = paste0(current_path, "pois2_UNIFORM.Rda"))
pois2_comp_UNIFORM <- compare_to_ground_truth(mu_pois2, 
                                              pois2_UNIFORM, 
                                              tau_pois2, 
                                              epsilon_pois2)$mean
save(pois2_comp_UNIFORM, file = paste0(current_path, "pois2_comp_UNIFORM.Rda"))

########################################################################
# Augmented-UCB (Mukherjee et al., 2017)

pois2_AugUCB <- para_bandit_sim_AugUCB(data = data_list_pois2, 
                                       rounds = 10000, 
                                       tau = tau_pois2)
save(pois2_AugUCB, file = paste0(current_path, "pois2_AugUCB.Rda"))
pois2_comp_AugUCB <- compare_to_ground_truth(mu_pois2, pois2_AugUCB, 
                                             tau_pois2, epsilon_pois2)$mean
save(pois2_comp_AugUCB, file = paste0(current_path, "pois2_comp_AugUCB.Rda"))

