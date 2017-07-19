sessions_threshold <- 0
set.seed(512)
selected_products <- sample(1:196, 15)
niter <- 5000
nrounds <- 1440 # Analyzing a data set of a day length
nahead <- 10080 # Computing CR a week ahead
source("/Users/timradtke/Desktop/thesis_data/get_dataset.R")

dim(pv_list[[1]])
pv_list[[1]]

pv_list[[1]]
pv_list[[5000]]
round(pv_next_means[[5000]],4)
round(pv_next_means[[1]],4)

#res <- matrix(ncol = 10, nrow = 5000)
#for(i in 1:5000) {
#  res[i,] <- round(pv_next_means[[i]] - pv_own_means[[i]], 4)
#}
#res

####################################################################

source("/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/AugUCB/run_from_data.R")
source("/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/AugUCB/run_parallelized_from_data.R")
current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/run_from_amo_data"

####################################################################

tau_amo1 <- 0.05
epsilon_amo1 <- 0.01
data_amo1 <- pv_list
rm(pv_list)
gc()

########################################################################
# Standard APT Algorithm
system.time(amo1_APT <- para_bandit_sim_APT(data = data_amo1, rounds = 1440, 
                                            tau = tau_amo1, epsilon = epsilon_amo1))
save(amo1_APT, file = paste0(current_path, "amo1_APT.Rda"))
amo1_comp_APT <- compare_to_cv_data(pv_next_means, amo1_APT, tau_amo1, epsilon_amo1)$mean
save(amo1_comp_APT, file = paste0(current_path, "amo1_comp_APT.Rda"))
rm(amo1_APT)
gc()

########################################################################

plot(c(0,1440), c(0, -4), type = "n")
lines(log(amo1_comp_APT), col = "red")
lines(log(loc2_comp_UNIFORM), col = "black")
lines(log(loc2_comp_PI), col = "darkgreen")
lines(log(loc2_comp_TTS), col = "lightgreen")
lines(log(loc2_comp_BUCB), col = "orange")
#lines(log(loc2_comp_KL), col = "lightblue")
lines(log(loc2_comp_KL_not_tau), col = "blue")
lines(log(loc2_comp_KL_horizon), col = "darkblue")
lines(log(loc2_comp_KL_horizon2), col = "darkred")
lines(log(loc2_comp_KL_horizon1), col = "red")
