########################################################################
# load functions

current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/AugUCB/"
source(paste0(current_path, "run_from_data.R"))
source(paste0(current_path, "run_parallelized_from_data.R"))
current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/synthetic_experiments/gaussian/"

########################################################################
# Create the data

set.seed(45844069)
tau_gexl <- 0.5
epsilon_gexl <- 0
(mu_gexl <- c(0.1,0.3,0.3,0.4,0.4,0.6,0.6,0.7,0.7,0.9))
#(sd_gexl <- c(1.5,1.5,1.5,1.5,1.5,4,4,4,4,4))
(sd_gexl <- c(0.8,0.8,0.8,0.8,0.8,2,2,2,2,2))

sd_gexl <- sd_gexl[order(mu_gexl)]
mu_gexl <- mu_gexl[order(mu_gexl)]

rbind(mu_gexl, sd_gexl)

plot(c(-2,2), c(0,2), type = "n")
for (i in 1:length(mu_gexl)) {
  lines(seq(-2,2,0.0001), dnorm(seq(-2,2,0.0001), mu_gexl[order(mu_gexl)][i], 
                                sd_gexl[order(mu_gexl)][i]),
        col = rainbow(50)[i])
}

plot(c(-2,2), c(0,2), type = "n")
for (i in 1:length(mu_gexl)) {
  lines(seq(-2,2,0.0001), dnorm(seq(-2,2,0.0001), mu_gexl[order(sd_gexl)][i], 
                                sd_gexl[order(sd_gexl)][i]),
        col = rainbow(50)[i])
}

data_gexl <- list()
for(j in 1:5000) {
  curr_data <- data.frame(rep(NA, times = 10000))
  for(i in 1:length(mu_gexl)) {
    curr_data[[i]] <- as.numeric(rnorm(10000, mu_gexl[i], sd_gexl[i]))
  }
  names(curr_data) <- paste0("V", rep(1:length(mu_gexl)))
  data_gexl[[j]] <- curr_data
}
gc()

########################################################################



########################################################################
# Gaussian Likelihood Ratio 2D-SLR

system.time(gexl_LR <- para_bandit_sim_LR_gaussian(data = data_gexl[1:500], 
                                                   rounds = 3000, 
                                                   tau = tau_gexl, 
                                                   epsilon = epsilon_gexl))

#save(gexl_LR, file = paste0(current_path, "gexl_LR.Rda"))
gexl_comp_LR <- compare_to_ground_truth(mu_gexl, gexl_LR, 
                                        tau_gexl, epsilon_gexl)$mean
#save(gexl_comp_LR, file = paste0(current_path, "gexl_comp_LR.Rda"))
#rm(gexl_LR)
gc()

########################################################################
# Anytime Parameter Free (Locatelli et. al, 2016)

system.time(gexl_APT <- para_bandit_sim_APT(data = data_gexl[1:500], 
                                            rounds = 3000, 
                                            tau = tau_gexl, 
                                            epsilon = epsilon_gexl,
                                            do_verbose = TRUE))
#save(gexl_APT, file = paste0(current_path, "gexl_APT.Rda"))
gexl_comp_APT <- compare_to_ground_truth(mu_gexl, gexl_APT, 
                                         tau_gexl, epsilon_gexl)$mean
#save(gexl_comp_APT, file = paste0(current_path, "gexl_comp_APT.Rda"))
gc()

########################################################################
# Empirical Variance Guided Algorithm (Zhong et al., 2017)

system.time(gexl_EVT <- para_bandit_sim_EVT(data = data_gexl[1:500], 
                                            rounds = 3000, 
                                            tau = tau_gexl, 
                                            epsilon = epsilon_gexl))
#save(gexl_EVT, file = paste0(current_path, "gexl_EVT.Rda"))
gexl_comp_EVT <- compare_to_ground_truth(mu_gexl, gexl_EVT,
                                         tau_gexl, epsilon_gexl)$mean
#save(gexl_comp_EVT, file = paste0(current_path, "gexl_comp_EVT.Rda"))
#rm(exp2_EVT)
gc()

########################################################################
# Uniform Sampling

system.time(gexl_UNIFORM <- para_bandit_sim_uniform(data = data_gexl[1:500], 
                                                    rounds = 6000))
#save(gexl_UNIFORM, file = paste0(current_path, "gexl_UNIFORM.Rda"))
gexl_comp_UNIFORM <- compare_to_ground_truth(mu_gexl, 
                                             gexl_UNIFORM, 
                                             tau_gexl, 
                                             epsilon_gexl)$mean
#save(gexl_comp_UNIFORM, file = paste0(current_path, "gexl_comp_UNIFORM.Rda"))
#rm(gexl_UNIFORM)
gc()



plot(c(0,8000), c(-6,0), type = "n")
lines(log(gexl_comp_LR))
lines(log(gexl_comp_APT), col = "blue")
lines(log(gexl_comp_EVT), col = "red")
lines(log(gexl_comp_UNIFORM), col = "green")
