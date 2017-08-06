########################################################################
# load functions

current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/AugUCB/"
source(paste0(current_path, "run_from_data.R"))
source(paste0(current_path, "run_parallelized_from_data.R"))
current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/synthetic_experiments/gaussian/"

########################################################################
# Create the data
# Setup as in Experiment 5 of Mukherjee et al.

set.seed(58459230)
tau_gaxEasy <- 0.5
epsilon_gaxEasy <- 0
mu_gaxEasy <- c(0.45, 0.55,
                0.45, 0.55,
                0.45, 0.55,
                0.45, 0.55)
sd_gaxEasy <- c(0.05, 0.1,
                0.1, 0.15,
                0.15, 0.2,
                0.05, 0.025)

plot(seq(-2,3,by=0.0001), dnorm(seq(-2,3,by=0.0001), mu_gaxEasy[1], sd_gaxEasy[1]), type = "l",
     ylim = c(0,5))
lines(seq(-2,3,by=0.0001), dnorm(seq(-2,3,by=0.0001), mu_gaxEasy[2], sd_gaxEasy[2]), col = "red")
abline(v = 0.5)

data_list_gaxEasy <- list()
for(j in 1:5000) {
  curr_data <- data.frame(rep(NA, times = 2000))
  for(i in 1:length(mu_gaxEasy)) {
    curr_data[[i]] <- as.numeric(rnorm(2000, mu_gaxEasy[i], sd_gaxEasy[i]))
  }
  names(curr_data) <- paste0("V", rep(1:length(mu_gaxEasy)))
  data_list_gaxEasy[[j]] <- curr_data
}
gc()

########################################################################

plot(c(0,2000), c(0, -8), type = "n")
lines(log(gaxEasy_comp_APT), col = "black")
lines(log(gaxEasy_comp_LR), col = "blue")

########################################################################
# Gaussian Likelihood Ratio
system.time(gaxEasy_LR <- para_bandit_sim_LR_gaussian(data = data_list_gaxEasy, 
                                                   rounds = 2000, 
                                                   tau = tau_gaxEasy, 
                                                   epsilon = epsilon_gaxEasy))



save(gaxEasy_LR, file = paste0(current_path, "gaxEasy_LR.Rda"))
gaxEasy_comp_LR <- compare_to_ground_truth(mu_gaxEasy, gaxEasy_LR, 
                                        tau_gaxEasy, epsilon_gaxEasy)$mean
save(gaxEasy_comp_LR, file = paste0(current_path, "gaxEasy_comp_LR.Rda"))
gc()

########################################################################

system.time(gaxEasy_APT <- para_bandit_sim_APT(data = data_list_gaxEasy, 
                                            rounds = 2000, 
                                            tau = tau_gaxEasy, 
                                            epsilon = epsilon_gaxEasy))
save(gaxEasy_APT, file = paste0(current_path, "gaxEasy_APT.Rda"))
gaxEasy_comp_APT <- compare_to_ground_truth(mu_gaxEasy, gaxEasy_APT, 
                                         tau_gaxEasy, epsilon_gaxEasy)$mean
save(gaxEasy_comp_APT, file = paste0(current_path, "gaxEasy_comp_APT.Rda"))
gc()

