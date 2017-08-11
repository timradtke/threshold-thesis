########################################################################
# load functions

current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/AugUCB/"
source(paste0(current_path, "run_from_data.R"))
source(paste0(current_path, "run_parallelized_from_data.R"))
current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/synthetic_experiments/exponential/"

########################################################################
# Create the data from exponential distributions

set.seed(93468734)
tau_exp2 <- 1
epsilon_exp2 <- 0
mu_exp2 <- rexp(20, 1)

plot(c(0,8), c(0,2), type = "n")
for(i in order(mu_exp2)) {
  lines(seq(0,8,0.001), dexp(seq(0,8,0.001), 1/mu_exp2[i]), col = rainbow(length(mu_exp2))[i])
}

data_list_exp2 <- list()
for(j in 1:5000) {
  curr_data <- data.frame(rep(NA, times = 4000))
  for(i in 1:length(mu_exp2)) {
    curr_data[[i]] <- as.numeric(rexp(4000, 1/mu_exp2[i]))
  }
  names(curr_data) <- paste0("V", rep(1:length(mu_exp2)))
  data_list_exp2[[j]] <- curr_data
}
gc()

########################################################################
load(paste0(current_path, "exp2_comp_LR.Rda"))
load(paste0(current_path, "exp2_comp_APT.Rda"))
load(paste0(current_path, "exp2_comp_EVT.Rda"))
load(paste0(current_path, "exp2_comp_UNIFORM.Rda"))

plot(c(0,4000), c(0, -8), type = "n")
abline(h=log(0.01), lty = 2)
lines(log(exp2_comp_UNIFORM), col = "black")
lines(log(exp2_comp_LR), col = rainbow(4)[1])
lines(log(exp2_comp_APT), col = rainbow(4)[2])
lines(log(exp2_comp_EVT), col = rainbow(4)[3])
lines(log(exp2_comp_AugUCB), col = rainbow(4)[4])

########################################################################
# Exponential-Distribution Likelihood Ratio
system.time(exp2_LR <- para_bandit_sim_LR_exponential(data = data_list_exp2, 
                                                      rounds = 4000, 
                                                      tau = tau_exp2, 
                                                      epsilon = epsilon_exp2))

save(exp2_LR, file = paste0(current_path, "exp2_LR.Rda"))
exp2_comp_LR <- compare_to_ground_truth(mu_exp2, exp2_LR, 
                                           tau_exp2, epsilon_exp2)$mean
save(exp2_comp_LR, file = paste0(current_path, "exp2_comp_LR.Rda"))
gc()

########################################################################
# Anytime Parameter Free (Locatelli et al., 2016)

system.time(exp2_APT <- para_bandit_sim_APT(data = data_list_exp2, 
                                               rounds = 4000, 
                                               tau = tau_exp2, 
                                               epsilon = epsilon_exp2))
save(exp2_APT, file = paste0(current_path, "exp2_APT.Rda"))
exp2_comp_APT <- compare_to_ground_truth(mu_exp2, exp2_APT, 
                                            tau_exp2, epsilon_exp2)$mean
save(exp2_comp_APT, file = paste0(current_path, "exp2_comp_APT.Rda"))
rm(exp2_APT)
gc()

########################################################################
# Empirical Variance Guided Algorithm (Zhong et al., 2017)

system.time(exp2_EVT <- para_bandit_sim_EVT(data = data_list_exp2, 
                                             rounds = 4000, 
                                             tau = tau_exp2, 
                                             epsilon = epsilon_exp2))
save(exp2_EVT, file = paste0(current_path, "exp2_EVT.Rda"))
exp2_comp_EVT <- compare_to_ground_truth(mu_exp2, exp2_EVT, 
                                            tau_exp2, epsilon_exp2)$mean
save(exp2_comp_EVT, file = paste0(current_path, "exp2_comp_EVT.Rda"))
rm(exp2_EVT)
gc()

########################################################################
# Uniform Sampling

system.time(exp2_UNIFORM <- para_bandit_sim_uniform(data = data_list_exp2, 
                                                    rounds = 4000))
save(exp2_UNIFORM, file = paste0(current_path, "exp2_UNIFORM.Rda"))
exp2_comp_UNIFORM <- compare_to_ground_truth(mu_exp2, 
                                               exp2_UNIFORM, 
                                               tau_exp2, 
                                               epsilon_exp2)$mean
save(exp2_comp_UNIFORM, file = paste0(current_path, "exp2_comp_UNIFORM.Rda"))
rm(exp2_UNIFORM)
gc()

########################################################################
# Augmented-UCB (Mukherjee et al., 2017)

#augtest <- AugUCB_from_tsdata(data = data_list_exp2[[1]], rounds = 2000, tau = tau_exp2,
#                   verbose = TRUE)

exp2_AugUCB <- para_bandit_sim_AugUCB(data = data_list_exp2, 
                                      rounds = 4000, 
                                      tau = tau_exp2)
save(exp2_AugUCB, file = paste0(current_path, "exp2_AugUCB.Rda"))
exp2_comp_AugUCB <- compare_to_ground_truth(mu_exp2, exp2_AugUCB, 
                                            tau_exp2, epsilon_exp2)$mean
save(exp2_comp_AugUCB, file = paste0(current_path, "exp2_comp_AugUCB.Rda"))
rm(exp2_AugUCB)
gc()
