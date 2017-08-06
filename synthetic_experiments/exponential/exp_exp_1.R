########################################################################
# load functions

current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/AugUCB/"
source(paste0(current_path, "run_from_data.R"))
source(paste0(current_path, "run_parallelized_from_data.R"))
current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/synthetic_experiments/exponential/"

########################################################################
# Create the data from exponential distributions

set.seed(58459230)
tau_expEasy <- 1
epsilon_expEasy <- 0
mu_expEasy <- c(0.01, 0.3, 0.6, 0.8, 0.8, 
                1.2, 1.2, 1.4, 2, 5)

plot(c(0,8), c(0,2), type = "n")
for(i in 1:length(mu_expEasy)) {
  lines(seq(0,8,0.001), dexp(seq(0,8,0.001), mu_expEasy[i]), col = rainbow(length(mu_expEasy))[i])
}

data_list_expEasy <- list()
for(j in 1:5000) {
  curr_data <- data.frame(rep(NA, times = 2000))
  for(i in 1:length(mu_expEasy)) {
    curr_data[[i]] <- as.numeric(rexp(2000, 1/mu_expEasy[i]))
  }
  names(curr_data) <- paste0("V", rep(1:length(mu_expEasy)))
  data_list_expEasy[[j]] <- curr_data
}
gc()

########################################################################

plot(c(0,2000), c(0, -8), type = "n")
abline(h=log(0.01), lty = 2)
lines(log(expEasy_comp_APT), col = "black")
lines(log(expEasy_comp_LR), col = "blue")
lines(log(expEasy_comp_EVT), col = "red")

########################################################################
# Gaussian Likelihood Ratio
system.time(expEasy_LR <- para_bandit_sim_LR_exponential(data = data_list_expEasy, 
                                                         rounds = 2000, 
                                                         tau = tau_expEasy, 
                                                         epsilon = epsilon_expEasy))



save(expEasy_LR, file = paste0(current_path, "expEasy_LR.Rda"))
expEasy_comp_LR <- compare_to_ground_truth(mu_expEasy, expEasy_LR, 
                                           tau_expEasy, epsilon_expEasy)$mean
save(expEasy_comp_LR, file = paste0(current_path, "expEasy_comp_LR.Rda"))
gc()

########################################################################

system.time(expEasy_APT <- para_bandit_sim_APT(data = data_list_expEasy, 
                                               rounds = 2000, 
                                               tau = tau_expEasy, 
                                               epsilon = epsilon_expEasy))
save(expEasy_APT, file = paste0(current_path, "expEasy_APT.Rda"))
expEasy_comp_APT <- compare_to_ground_truth(mu_expEasy, expEasy_APT, 
                                            tau_expEasy, epsilon_expEasy)$mean
save(expEasy_comp_APT, file = paste0(current_path, "expEasy_comp_APT.Rda"))
gc()

########################################################################

system.time(expEasy_EVT <- para_bandit_sim_EVT(data = data_list_expEasy, 
                                               rounds = 2000, 
                                               tau = tau_expEasy, 
                                               epsilon = epsilon_expEasy))
save(expEasy_EVT, file = paste0(current_path, "expEasy_EVT.Rda"))
expEasy_comp_EVT <- compare_to_ground_truth(mu_expEasy, expEasy_EVT, 
                                            tau_expEasy, epsilon_expEasy)$mean
save(expEasy_comp_EVT, file = paste0(current_path, "expEasy_comp_EVT.Rda"))
gc()

