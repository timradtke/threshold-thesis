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
tau_gax3 <- 0.5
epsilon_gax3 <- 0
mu_gax3 <- c(rep(0.4, 3),
             rep(0.6, 3),
             rep(0.3, 4))
sd_gax3 <- c(rep(sqrt(0.4), 3),
             rep(sqrt(1), 3),
             sqrt(runif(4,0.2,0.4)))

plot(seq(-2,3,by=0.0001), dnorm(seq(-2,3,by=0.0001), 0.4, sqrt(0.4)), type = "l",
     ylim = c(0,2))
lines(seq(-2,3,by=0.0001), dnorm(seq(-2,3,by=0.0001), 0.6, sqrt(1)), col = "red")
abline(v = 0.5)

data_list_gax3 <- list()
for(j in 1:1000) {
  curr_data <- data.frame(rep(NA, times = 10000))
  for(i in 1:length(mu_gax3)) {
    curr_data[[i]] <- as.numeric(rnorm(10000, mu_gax3[i], sd_gax3[i]))
  }
  names(curr_data) <- paste0("V", rep(1:length(mu_gax3)))
  data_list_gax3[[j]] <- curr_data
}
gc()

########################################################################

plot(c(0,10000), c(0, -8), type = "n")
lines(log(gax3_comp_APT), col = "black")
lines(log(gax3_comp_LR), col = "red")

########################################################################
# Gaussian Likelihood Ratio
system.time(gax3_LR <- para_bandit_sim_LR_gaussian(data = data_list_gax3, 
                                                   rounds = 10000, 
                                                   tau = tau_gax3, 
                                                   epsilon = epsilon_gax3))

save(gax3_LR, file = paste0(current_path, "gax3_LR.Rda"))
gax3_comp_LR <- compare_to_ground_truth(mu_gax3, gax3_LR, 
                                        tau_gax3, epsilon_gax3)$mean
save(gax3_comp_LR, file = paste0(current_path, "gax3_comp_LR.Rda"))
rm(gax3_LR)
gc()

########################################################################

system.time(gax3_APT <- para_bandit_sim_APT(data = data_list_gax3, 
                                            rounds = 10000, 
                                            tau = tau_gax3, 
                                            epsilon = epsilon_gax3,
                                            do_verbose = TRUE))
save(gax3_APT, file = paste0(current_path, "gax3_APT.Rda"))
gax3_comp_APT <- compare_to_ground_truth(mu_gax3, gax3_APT, 
                                         tau_gax3, epsilon_gax3)$mean
save(gax3_comp_APT, file = paste0(current_path, "gax3_comp_APT.Rda"))
gc()

