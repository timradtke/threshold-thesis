# run a simulation based on a simple example that can hopefully serve to
# calibrate the algorithms
########################################################################
# load functions

current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/AugUCB/"
source(paste0(current_path, "run_from_data.R"))
source(paste0(current_path, "run_parallelized_from_data.R"))
current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/synthetic_experiments/"

########################################################################
# Create the data

mu_loc6g <- c(0.01, 0.06, 0.11, 0.16,
              0.21, 0.29,
              0.34, 0.39, 0.44, 0.49)
get_ber_var <- function(x) x*(1-x)
sd_loc6g <- sqrt(get_ber_var(mu_loc6g))
tau_loc6g <- 0.25
epsilon_loc6g <- 0

rbind(mu_loc6g, sd_loc6g)

plot(c(-2,2), c(0,2), type = "n")
for (i in 1:length(mu_loc6g)) {
  lines(seq(-2,2,0.0001), dnorm(seq(-2,2,0.0001), mu_loc6g[order(mu_loc6g)][i], 
                                sd_loc6g[order(mu_loc6g)][i]),
        col = rainbow(50)[i])
}
abline(v = tau_loc6g, lty = 2)

plot(c(-2,2), c(0,2), type = "n")
for (i in 1:length(mu_loc6g)) {
  lines(seq(-2,2,0.0001), dnorm(seq(-2,2,0.0001), mu_loc6g[order(sd_loc6g)][i], 
                                sd_loc6g[order(sd_loc6g)][i]),
        col = rainbow(50)[i])
}
abline(v = tau_loc6g, lty = 2)

data_loc6g <- list()
for(j in 1:5000) {
  curr_data <- data.frame(rep(NA, times = 10000))
  for(i in 1:length(mu_loc6g)) {
    curr_data[[i]] <- as.numeric(rnorm(10000, mu_loc6g[i], sd_loc6g[i]))
  }
  names(curr_data) <- paste0("V", rep(1:length(mu_loc6g)))
  data_loc6g[[j]] <- curr_data
}
gc()

########################################################################

load(file = paste0(current_path, "loc6g_LR.Rda"))
load(file = paste0(current_path, "loc6g_APT.Rda"))
rbind(table(loc6g_LR[[6]]$arm_sequence), mu_loc6g, sd_loc6g)
rbind(table(loc6g_APT[[2]]$arm_sequence), mu_loc6g, sd_loc6g)

########################################################################
# Gaussian Likelihood Ratio 2D-SLR

system.time(loc6g_LR <- para_bandit_sim_LR_gaussian(data = data_loc6g[1:500], 
                                                   rounds = 3000, 
                                                   tau = tau_loc6g, 
                                                   epsilon = epsilon_loc6g))

save(loc6g_LR, file = paste0(current_path, "loc6g_LR.Rda"))
loc6g_comp_LR <- compare_to_ground_truth(mu_loc6g, loc6g_LR, 
                                        tau_loc6g, epsilon_loc6g)$mean
save(loc6g_comp_LR, file = paste0(current_path, "loc6g_comp_LR.Rda"))
rm(loc6g_LR)
gc()

########################################################################
# Anytime Parameter Free (Locatelli et. al, 2016)

system.time(loc6g_APT <- para_bandit_sim_APT(data = data_loc6g[1:500], 
                                            rounds = 3000, 
                                            tau = tau_loc6g, 
                                            epsilon = epsilon_loc6g,
                                            do_verbose = TRUE))
save(loc6g_APT, file = paste0(current_path, "loc6g_APT.Rda"))
loc6g_comp_APT <- compare_to_ground_truth(mu_loc6g, loc6g_APT, 
                                         tau_loc6g, epsilon_loc6g)$mean
save(loc6g_comp_APT, file = paste0(current_path, "loc6g_comp_APT.Rda"))
rm(loc6g_APT)
gc()

########################################################################
# Uniform Sampling

system.time(loc6g_UNIFORM <- para_bandit_sim_uniform(data = data_loc6g[1:500], 
                                                    rounds = 3000))
save(loc6g_UNIFORM, file = paste0(current_path, "loc6g_UNIFORM.Rda"))
loc6g_comp_UNIFORM <- compare_to_ground_truth(mu_loc6g, 
                                             loc6g_UNIFORM, 
                                             tau_loc6g, 
                                             epsilon_loc6g)$mean
save(loc6g_comp_UNIFORM, file = paste0(current_path, "loc6g_comp_UNIFORM.Rda"))
rm(loc6g_UNIFORM)
gc()



plot(c(0,8000), c(-8,0), type = "n")
lines(log(loc6g_comp_LR))
lines(log(loc6g_comp_APT), col = "blue")
lines(log(loc6g_comp_UNIFORM), col = "red")
