########################################################################
# load functions

current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/AugUCB/"
source(paste0(current_path, "run_from_data.R"))
source(paste0(current_path, "run_parallelized_from_data.R"))
current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/synthetic_experiments/gaussian/"

########################################################################
# Create the data

set.seed(45844069)
tau_gex <- 0
epsilon_gex <- 0
(mu_gex <- rnorm(20, 0, 0.3))
(sd_gex <- sqrt(1/rgamma(50, 2, 0.1)))

sd_gex <- sd_gex[order(mu_gex)]
mu_gex <- mu_gex[order(mu_gex)]

rbind(mu_gex, sd_gex)

plot(c(-2,2), c(0,2), type = "n")
for (i in 1:length(mu_gex)) {
  lines(seq(-2,2,0.0001), dnorm(seq(-2,2,0.0001), mu_gex[order(mu_gex)][i], 
                                sd_gex[order(mu_gex)][i]),
        col = rainbow(50)[i])
}

plot(c(-2,2), c(0,2), type = "n")
for (i in 1:length(mu_gex)) {
  lines(seq(-2,2,0.0001), dnorm(seq(-2,2,0.0001), mu_gex[order(sd_gex)][i], 
                                sd_gex[order(sd_gex)][i]),
        col = rainbow(50)[i])
}

data_gex <- list()
for(j in 1:5000) {
  curr_data <- data.frame(rep(NA, times = 10000))
  for(i in 1:length(mu_gex)) {
    curr_data[[i]] <- as.numeric(rnorm(10000, mu_gex[i], sd_gex[i]))
  }
  names(curr_data) <- paste0("V", rep(1:length(mu_gex)))
  data_gex[[j]] <- curr_data
}
gc()

########################################################################



########################################################################
# Gaussian Likelihood Ratio 2D-SLR

system.time(gex_LR <- para_bandit_sim_LR_gaussian(data = data_gex[1:1000], 
                                                   rounds = 8000, 
                                                   tau = tau_gex, 
                                                   epsilon = epsilon_gex))

#save(gex_LR, file = paste0(current_path, "gex_LR.Rda"))
gex_comp_LR <- compare_to_ground_truth(mu_gex, gex_LR, 
                                        tau_gex, epsilon_gex)$mean
#save(gex_comp_LR, file = paste0(current_path, "gex_comp_LR.Rda"))
rm(gex_LR)
gc()

########################################################################
# Anytime Parameter Free (Locatelli et. al, 2016)

system.time(gex_APT <- para_bandit_sim_APT(data = data_gex[1:1000], 
                                            rounds = 8000, 
                                            tau = tau_gex, 
                                            epsilon = epsilon_gex,
                                            do_verbose = TRUE))
#save(gex_APT, file = paste0(current_path, "gex_APT.Rda"))
gex_comp_APT <- compare_to_ground_truth(mu_gex, gex_APT, 
                                         tau_gex, epsilon_gex)$mean
#save(gex_comp_APT, file = paste0(current_path, "gex_comp_APT.Rda"))
gc()



plot(c(0,8000), c(-8,0), type = "n")
lines(log(gex_comp_LR))
lines(log(gex_comp_APT), col = "blue")
