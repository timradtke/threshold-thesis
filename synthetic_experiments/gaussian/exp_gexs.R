########################################################################
# load functions

current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/AugUCB/"
source(paste0(current_path, "run_from_data.R"))
source(paste0(current_path, "run_parallelized_from_data.R"))
current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/synthetic_experiments/gaussian/"

########################################################################
# Create the data

set.seed(45844069)
tau_gexs <- 0.5
epsilon_gexs <- 0
(mu_gexs <- c(0.1,0.3,0.4,0.4,0.45,0.55,0.6,0.6,0.7,0.9))
(sd_gexs <- c(0.9,0.9,0.9,0.2,0.9,0.45,0.1,0.45,0.45,0.45))

sd_gexs <- sd_gexs[order(mu_gexs)]
mu_gexs <- mu_gexs[order(mu_gexs)]

rbind(mu_gexs, sd_gexs)

plot(c(-2,2), c(0,2), type = "n")
for (i in 1:length(mu_gexs)) {
  lines(seq(-2,2,0.0001), dnorm(seq(-2,2,0.0001), mu_gexs[order(mu_gexs)][i], 
                                sd_gexs[order(mu_gexs)][i]),
        col = rainbow(50)[i])
}
abline(v = tau_gexs, lty = 2)

plot(c(-2,2), c(0,2), type = "n")
for (i in 1:length(mu_gexs)) {
  lines(seq(-2,2,0.0001), dnorm(seq(-2,2,0.0001), mu_gexs[order(sd_gexs)][i], 
                                sd_gexs[order(sd_gexs)][i]),
        col = rainbow(50)[i])
}
abline(v = tau_gexs, lty = 2)

data_gexs <- list()
for(j in 1:5000) {
  curr_data <- data.frame(rep(NA, times = 10000))
  for(i in 1:length(mu_gexs)) {
    curr_data[[i]] <- as.numeric(rnorm(10000, mu_gexs[i], sd_gexs[i]))
  }
  names(curr_data) <- paste0("V", rep(1:length(mu_gexs)))
  data_gexs[[j]] <- curr_data
}
gc()

########################################################################

load(file = paste0(current_path, "gexs_LR.Rda"))
load(file = paste0(current_path, "gexs_APT.Rda"))
rbind(table(gexs_LR[[6]]$arm_sequence), mu_gexs, sd_gexs)
rbind(table(gexs_APT[[2]]$arm_sequence), mu_gexs, sd_gexs)

########################################################################
# Gaussian Likelihood Ratio 2D-SLR

system.time(gexs_LR <- para_bandit_sim_LR_gaussian(data = data_gexs[1:1000], 
                                                  rounds = 6000, 
                                                  tau = tau_gexs, 
                                                  epsilon = epsilon_gexs))

save(gexs_LR, file = paste0(current_path, "gexs_LR.Rda"))
gexs_comp_LR <- compare_to_ground_truth(mu_gexs, gexs_LR, 
                                       tau_gexs, epsilon_gexs)$mean
save(gexs_comp_LR, file = paste0(current_path, "gexs_comp_LR.Rda"))
rm(gexs_LR)
gc()

########################################################################
# Anytime Parameter Free (Locatelli et. al, 2016)

system.time(gexs_APT <- para_bandit_sim_APT(data = data_gexs[1:1000], 
                                           rounds = 6000, 
                                           tau = tau_gexs, 
                                           epsilon = epsilon_gexs,
                                           do_verbose = TRUE))
save(gexs_APT, file = paste0(current_path, "gexs_APT.Rda"))
gexs_comp_APT <- compare_to_ground_truth(mu_gexs, gexs_APT, 
                                        tau_gexs, epsilon_gexs)$mean
save(gexs_comp_APT, file = paste0(current_path, "gexs_comp_APT.Rda"))
rm(gexs_APT)
gc()

########################################################################
# Uniform Sampling

system.time(gexs_UNIFORM <- para_bandit_sim_uniform(data = data_gexs[1:1000], 
                                                    rounds = 6000))
save(gexs_UNIFORM, file = paste0(current_path, "gexs_UNIFORM.Rda"))
gexs_comp_UNIFORM <- compare_to_ground_truth(mu_gexs, 
                                             gexs_UNIFORM, 
                                             tau_gexs, 
                                             epsilon_gexs)$mean
save(gexs_comp_UNIFORM, file = paste0(current_path, "gexs_comp_UNIFORM.Rda"))
rm(gexs_UNIFORM)
gc()



plot(c(0,8000), c(-8,0), type = "n")
lines(log(gexs_comp_LR))
lines(log(gexs_comp_APT), col = "blue")
lines(log(gexs_comp_UNIFORM), col = "red")
