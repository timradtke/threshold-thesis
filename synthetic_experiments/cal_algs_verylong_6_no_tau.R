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

mean_loc6nt <- c(0.001, 0.005, 0.01, 0.015,
                 0.04, 0.06,
                 0.085, 0.09, 0.095, 0.099)
tau_loc6nt <- 0.05
epsilon_loc6nt <- 0
H_loc6nt <- get_complexity(mean_loc6nt, tau_loc6nt, epsilon_loc6nt)

plot(mean_loc6nt, rep(1,10), main = paste0("Complexity of ", round(H_loc6nt,2)))
abline(v=tau_loc6nt)
abline(v=tau_loc6nt+epsilon_loc6nt, lty=2)
abline(v=tau_loc6nt-epsilon_loc6nt, lty=2)

data_list6 <- list()
set.seed(8247502)
for(j in 1:5000) {
  curr_data <- data.frame(rep(NA, times = 7000))
  for(i in 1:length(mean_loc6nt)) {
    curr_data[[i]] <- as.numeric(purrr::rbernoulli(7000, p  = mean_loc6nt[i]))
  }
  names(curr_data) <- paste0("V", rep(1:length(mean_loc6nt)))
  data_list6[[j]] <- curr_data
}

########################################################################

load(paste0(current_path, "/loc6nt_BUCB_horizon_7000.Rda"))
load(paste0(current_path, "/loc6nt_APT_7000.Rda"))
load(paste0(current_path, "/loc6nt_LR_7000.Rda"))
load(paste0(current_path, "/loc6nt_LRD.Rda"))

table(loc6nt_LR_7000[[1]]$arm_sequence)
table(loc6nt_LR_7000[[2]]$arm_sequence)
table(loc6nt_LR_7000[[3]]$arm_sequence)
table(loc6nt_LR_7000[[4]]$arm_sequence)
table(loc6nt_LR_7000[[5]]$arm_sequence)
tail(loc6nt_LR_7000[[5]]$mean_storage)

library(ggplot2)
library(ggjoy)
library(dplyr)
library(tidyr)

pulls_of_arm_6_LR <- vector(length = 5000)
for(i in 1:5000){
  pulls_of_arm_6_LR[i] <- table(loc6nt_LR_7000[[i]]$arm_sequence)[6]
}
summary(pulls_of_arm_6_LR)
quantile(pulls_of_arm_6_LR, 0.05)

pulls_of_arm_6_APT <- vector(length = 5000)
for(i in 1:5000){
  pulls_of_arm_6_APT[i] <- table(loc6nt_APT_7000[[i]]$arm_sequence)[6]
}
summary(pulls_of_arm_6_APT)
quantile(pulls_of_arm_6_APT, 0.05)

pulls_of_arm_6_BUCB <- vector(length = 5000)
for(i in 1:5000){
  pulls_of_arm_6_BUCB[i] <- table(loc6nt_BUCB_horizon_7000[[i]]$arm_sequence)[6]
}
summary(pulls_of_arm_6_BUCB)
quantile(pulls_of_arm_6_BUCB, 0.05)

pulls_of_arm_6_LRD <- vector(length = 5000)
for(i in 1:5000){
  pulls_of_arm_6_LRD[i] <- table(loc6nt_LRD[[i]]$arm_sequence)[6]
}

pulls_arm_6 <- data.frame(APT = pulls_of_arm_6_APT,
                          SLR = pulls_of_arm_6_LR,
                          SLRD = pulls_of_arm_6_LRD,
                          BUCB = pulls_of_arm_6_BUCB)

pulls_arm_6 %>% gather(key = "Algorithm", value = "Pulls") %>%
  ggplot(aes(x = Pulls, group = Algorithm)) +
  geom_histogram() +
  facet_grid(Algorithm~.)

pulls_of_arm_1_LR <- vector(length = 5000)
for(i in 1:5000){
  pulls_of_arm_1_LR[i] <- table(loc6nt_LR_7000[[i]]$arm_sequence)[1]
}

pulls_of_arm_1_LRD <- vector(length = 5000)
for(i in 1:5000){
  pulls_of_arm_1_LRD[i] <- table(loc6nt_LRD[[i]]$arm_sequence)[1]
}

pulls_of_arm_1_APT <- vector(length = 5000)
for(i in 1:5000){
  pulls_of_arm_1_APT[i] <- table(loc6nt_APT_7000[[i]]$arm_sequence)[1]
}

pulls_of_arm_1_BUCB <- vector(length = 5000)
for(i in 1:5000){
  pulls_of_arm_1_BUCB[i] <- table(loc6nt_BUCB_horizon_7000[[i]]$arm_sequence)[1]
}

pulls_arm_1 <- data.frame(APT = pulls_of_arm_1_APT,
                          SLR = pulls_of_arm_1_LR,
                          SLRD = pulls_of_arm_1_LRD,
                          BUCB = pulls_of_arm_1_BUCB)

pulls_arm_1 %>% gather(key = "Algorithm", value = "Pulls") %>%
  ggplot(aes(x = Pulls, group = Algorithm)) +
  geom_histogram() +
  facet_grid(Algorithm~.)

pulls_of_arm_10_LR <- vector(length = 5000)
for(i in 1:5000){
  pulls_of_arm_10_LR[i] <- table(loc6nt_LR_7000[[i]]$arm_sequence)[10]
}

pulls_of_arm_10_LRD <- vector(length = 5000)
for(i in 1:5000){
  pulls_of_arm_10_LRD[i] <- table(loc6nt_LRD[[i]]$arm_sequence)[10]
}

pulls_of_arm_10_APT <- vector(length = 5000)
for(i in 1:5000){
  pulls_of_arm_10_APT[i] <- table(loc6nt_APT_7000[[i]]$arm_sequence)[10]
}

pulls_of_arm_10_BUCB <- vector(length = 5000)
for(i in 1:5000){
  pulls_of_arm_10_BUCB[i] <- table(loc6nt_BUCB_horizon_7000[[i]]$arm_sequence)[10]
}

pulls_arm_10 <- data.frame(APT = pulls_of_arm_10_APT,
                          SLR = pulls_of_arm_10_LR,
                          SLRD = pulls_of_arm_10_LRD,
                          BUCB = pulls_of_arm_10_BUCB)

pulls_arm_10 %>% gather(key = "Algorithm", value = "Pulls") %>%
  ggplot(aes(x = Pulls, group = Algorithm)) +
  geom_histogram(binwidth = 50) +
  facet_grid(Algorithm~.)

pulls_arm_1_long <- pulls_arm_1 %>% gather(key = "Algorithm", value = "Pulls")
pulls_arm_6_long <- pulls_arm_6 %>% gather(key = "Algorithm", value = "Pulls")
pulls_arm_10_long <- pulls_arm_10 %>% gather(key = "Algorithm", value = "Pulls")

save(pulls_arm_1_long, pulls_arm_6_long, pulls_arm_10_long,
     file = paste0(current_path, "loc6nt_pulls_of_arms_vis.Rda"))
load(file = paste0(current_path, "loc6nt_pulls_of_arms_vis.Rda"))

pulls_arm_6_long %>%
  ggplot(aes(x = Pulls)) +
  geom_histogram(binwidth = 100) +
  facet_grid(.~Algorithm, scales = "free_y") +
  geom_vline(aes(xintercept = 7000/10), linetype = 2) +
  labs(x = "Pulls", y = "Count", 
       title = "Distribution of Pulls of Arm 6 Across Simulations",
       subtitle = "Binwidth = 100. Budget = 7000. 5000 Simulations.
Vertical line indicates pulls by uniform sampling given 10 arms in experiment.") +
  theme_bw()

cbind.data.frame(
  Arm = rep(c("01", "10"), each = 20000),
  rbind.data.frame(pulls_arm_1_long, pulls_arm_10_long),
  stringsAsFactors = FALSE
) %>% 
  ggplot(aes(x = Pulls, group = Arm)) +
  geom_histogram(binwidth = 50) +
  facet_grid(Arm~Algorithm, scales = "free_y")+
  geom_vline(aes(xintercept = 7000/10), linetype = 2) +
  labs(x = "Pulls", y = "Count", 
       title = "Distribution of Pulls of Arms 1 and 10 Across Simulations",
       subtitle = "Binwidth = 50. Budget = 7000. 5000 Simulations.
Vertical line indicates pulls by uniform sampling given 10 arms in experiment.") +
  theme_bw()

########################################################################
# Empirical Variance Guided Algorithm (Zhong et al., 2017)

loc6nt_EVT <- para_bandit_sim_EVT(data = data_list6, 
                                  rounds = 7000, 
                                  tau = tau_loc6nt, 
                                  epsilon = epsilon_loc6nt)
save(loc6nt_EVT, file = paste0(current_path, "loc6nt_EVT.Rda"))
loc6nt_comp_EVT <- compare_to_ground_truth(mean_loc6nt, loc6nt_EVT, 
                                           tau_loc6nt, epsilon_loc6nt)$mean
save(loc6nt_comp_EVT, file = paste0(current_path, "loc6nt_comp_EVT.Rda"))
rm(loc6nt_EVT)
gc()

########################################################################
# Plot the results

load(paste0(current_path, "/loc6nt_comp_BUCB_horizon_7000.Rda"))
load(paste0(current_path, "/loc6nt_comp_BUCB_horizon_c5_7000.Rda"))
load(paste0(current_path, "/loc6nt_comp_BUCB_horizon_c15_7000.Rda"))
load(paste0(current_path, "/loc6nt_comp_AugUCB_7000.Rda"))
load(paste0(current_path, "/loc6nt_comp_KLUCB370_5000.Rda"))
load(paste0(current_path, "/loc6nt_comp_UNIFORM_7000.Rda"))
load(paste0(current_path, "/loc6nt_comp_APT_7000.Rda"))
load(paste0(current_path, "/loc6nt_comp_LR_7000.Rda"))

plot(c(0,7000), c(0, -5), type = "n")
lines(log(loc6nt_comp_BUCB_7000), col = "black")
lines(log(loc6nt_comp_BUCB_horizon_7000), col = "blue")
lines(log(loc6nt_comp_BUCB_horizon_c5_7000), col = "lightblue")
lines(log(loc6nt_comp_BUCB_horizon_c15_7000), col = "darkblue")
#lines(log(loc6nt_comp_BUCB_horizon), col = "lightblue")
#lines(log(loc6nt_comp_APT), col = "pink")
lines(log(loc6nt_comp_LR_7000), col = "green")
lines(log(loc6nt_comp_LRD), col = "darkgreen")
lines(log(loc6nt_comp_APT_7000), col = "red")
lines(log(loc6nt_comp_UNIFORM_7000), col = "black")
abline(h=log(0.1))
lines(log(loc6nt_comp_PI), col = "orange")
#lines(log(loc6nt_comp_AugUCB), col = "grey")
lines(log(loc6nt_comp_AugUCB_7000), col = "grey")
#lines(log(loc6nt_comp_KLUCB7), col = "green")
#lines(log(loc6nt_comp_KLUCB30), col = "green")
#lines(log(loc6nt_comp_KLUCB370), col = "darkgreen")
#lines(log(loc6nt_comp_KLUCB370_5000), col = "darkgreen")

########################################################################
# Empirical Variance Guided Algorithm (Zhong et al., 2017)

loc6nt_EVT <- para_bandit_sim_EVT(data = data_list6, 
                                rounds = 7000, 
                                tau = tau_loc6nt, 
                                epsilon = epsilon_loc6nt)
save(loc6nt_EVT, file = paste0(current_path, "loc6nt_EVT.Rda"))
loc6nt_comp_EVT <- compare_to_ground_truth(mean_loc6nt, loc6nt_EVT, 
                                           tau_loc6nt, epsilon_loc6nt)$mean
save(loc6nt_comp_EVT, file = paste0(current_path, "loc6nt_comp_EVT.Rda"))
rm(loc6nt_EVT)
gc()

########################################################################
# Probability of Improvement

loc6nt_PI <- para_bandit_sim_PI(data = data_list6, 
                                rounds = 7000,
                                 tau = tau_loc6nt, 
                                 epsilon = epsilon_loc6nt, 
                                 alpha = tau_loc6nt, 
                                 beta = 1-tau_loc6nt,
                                tadj = FALSE)
save(loc6nt_PI, 
     file = paste0(current_path, "loc6nt_PI.Rda"))
loc6nt_comp_PI <- compare_to_ground_truth(mean_loc6nt, 
                                          loc6nt_PI,
                                          tau_loc6nt,
                                          epsilon_loc6nt)$mean
save(loc6nt_comp_PI, file = paste0(current_path, 
                                   "loc6nt_comp_PI.Rda"))


########################################################################

loc6nt_BUCB_horizon_7000 <- para_bandit_sim_bucb(data = data_list6, 
                                                 rounds = 7000, 
                                                 rate = "inverse_horizon",
                                                 tau = tau_loc6nt, 
                                                 epsilon = epsilon_loc6nt, 
                                                 alpha = tau_loc6nt, 
                                                 beta = 1-tau_loc6nt)
save(loc6nt_BUCB_horizon_7000, 
     file = paste0(current_path, "loc6nt_BUCB_horizon_7000.Rda"))
loc6nt_comp_BUCB_horizon_7000 <- compare_to_ground_truth(mean_loc6nt, 
                                                         loc6nt_BUCB_horizon_7000,
                                                         tau_loc6nt,
                                                         epsilon_loc6nt)$mean
save(loc6nt_comp_BUCB_horizon_7000, file = paste0(current_path, 
                                                  "loc6nt_comp_BUCB_horizon_7000.Rda"))

########################################################################
# Try a different exploration parameter for BUCB

loc6nt_BUCB_7000 <- para_bandit_sim_bucb(data = data_list6, 
                                          rounds = 7000, 
                                          rate = "inverse",
                                          tau = tau_loc6nt, 
                                          epsilon = epsilon_loc6nt, 
                                          alpha = tau_loc6nt, 
                                          beta = 1-tau_loc6nt)
save(loc6nt_BUCB_7000,
     file = paste0(current_path, "loc6nt_BUCB_7000.Rda"))
loc6nt_comp_BUCB_7000 <- compare_to_ground_truth(mean_loc6nt, 
                                                 loc6nt_BUCB_7000,
                                                  tau_loc6nt,
                                                  epsilon_loc6nt)$mean
save(loc6nt_comp_BUCB_7000, file = paste0(current_path, 
                                           "loc6nt_comp_BUCB_7000.Rda"))
rm(loc6nt_BUCB_7000)

########################################################################
# Try a different exploration parameter for BUCB

loc6nt_BUCB_horizon_c5_7000 <- para_bandit_sim_bucb(data = data_list6, 
                                                 rounds = 7000, 
                                                 rate = "inverse_horizon_c",
                                                 tau = tau_loc6nt, 
                                                 epsilon = epsilon_loc6nt, 
                                                 alpha = tau_loc6nt, 
                                                 beta = 1-tau_loc6nt,
                                                 const = 5)
save(loc6nt_BUCB_horizon_c5_7000, 
     file = paste0(current_path, "loc6nt_BUCB_horizon_c5_7000.Rda"))
loc6nt_comp_BUCB_horizon_c5_7000 <- compare_to_ground_truth(mean_loc6nt, 
                                                         loc6nt_BUCB_horizon_c5_7000,
                                                         tau_loc6nt,
                                                         epsilon_loc6nt)$mean
save(loc6nt_comp_BUCB_horizon_c5_7000, file = paste0(current_path, 
                                                  "loc6nt_comp_BUCB_horizon_c5_7000.Rda"))
rm(loc6nt_BUCB_horizon_c5_7000)

########################################################################
# Try a different exploration parameter for BUCB

loc6nt_BUCB_horizon_c15_7000 <- para_bandit_sim_bucb(data = data_list6, 
                                                    rounds = 7000, 
                                                    rate = "inverse_horizon_c",
                                                    tau = tau_loc6nt, 
                                                    epsilon = epsilon_loc6nt, 
                                                    alpha = tau_loc6nt, 
                                                    beta = 1-tau_loc6nt,
                                                    const = 1/5)
save(loc6nt_BUCB_horizon_c15_7000, 
     file = paste0(current_path, "loc6nt_BUCB_horizon_c15_7000.Rda"))
loc6nt_comp_BUCB_horizon_c15_7000 <- compare_to_ground_truth(mean_loc6nt, 
                                                            loc6nt_BUCB_horizon_c15_7000,
                                                            tau_loc6nt,
                                                            epsilon_loc6nt)$mean
save(loc6nt_comp_BUCB_horizon_c15_7000, file = paste0(current_path, 
                                                     "loc6nt_comp_BUCB_horizon_c15_7000.Rda"))
rm(loc6nt_BUCB_horizon_c15_7000)

########################################################################
# Standard AugUCB
loc6nt_AugUCB_7000 <- para_bandit_sim_AugUCB(data = data_list6, 
                                             rounds = 7000, 
                                             tau = tau_loc6nt)
save(loc6nt_AugUCB_7000, file = paste0(current_path, "loc6nt_AugUCB_7000.Rda"))
loc6nt_comp_AugUCB_7000 <- compare_to_ground_truth(mean_loc6nt, loc6nt_AugUCB_7000, 
                                                   tau_loc6nt, epsilon_loc6nt)$mean
save(loc6nt_comp_AugUCB_7000, file = paste0(current_path, 
                                            "loc6nt_comp_AugUCB_7000.Rda"))

########################################################################
# Standard APT
system.time(loc6nt_APT_7000 <- para_bandit_sim_APT(data = data_list6, 
                                                   rounds = 7000, 
                                                   tau = tau_loc6nt, 
                                                   epsilon = epsilon_loc6nt))

save(loc6nt_APT_7000, file = paste0(current_path, "loc6nt_APT_7000.Rda"))
loc6nt_comp_APT_7000 <- compare_to_ground_truth(mean_loc6nt, loc6nt_APT_7000, 
                                                tau_loc6nt, epsilon_loc6nt)$mean
save(loc6nt_comp_APT_7000, file = paste0(current_path, 
                                         "loc6nt_comp_APT_7000.Rda"))

########################################################################
# Standard Likelihood Ratio
system.time(loc6nt_LR_7000 <- para_bandit_sim_LR(data = data_list6, 
                                                 rounds = 7000, 
                                                 tau = tau_loc6nt, 
                                                 epsilon = epsilon_loc6nt))

save(loc6nt_LR_7000, file = paste0(current_path, "loc6nt_LR_7000.Rda"))
loc6nt_comp_LR_7000 <- compare_to_ground_truth(mean_loc6nt, loc6nt_LR_7000, 
                                                tau_loc6nt, epsilon_loc6nt)$mean
save(loc6nt_comp_LR_7000, file = paste0(current_path, 
                                        "loc6nt_comp_LR_7000.Rda"))

########################################################################
# Standard Likelihood Ratio with D-Tracking rule
system.time(loc6nt_LRD <- para_bandit_sim_LRD(data = data_list6, 
                                               rounds = 7000, 
                                               tau = tau_loc6nt, 
                                               epsilon = epsilon_loc6nt,
                                               do_verbose = TRUE))

save(loc6nt_LRD, file = paste0(current_path, "loc6nt_LRD.Rda"))
loc6nt_comp_LRD <- compare_to_ground_truth(mean_loc6nt, loc6nt_LRD, 
                                            tau_loc6nt, epsilon_loc6nt)$mean
save(loc6nt_comp_LRD, file = paste0(current_path, 
                                    "loc6nt_comp_LRD.Rda"))

########################################################################

system.time(loc6nt_UNIFORM_7000 <- para_bandit_sim_uniform(data = data_list6, 
                                                         rounds = 7000))
save(loc6nt_UNIFORM_7000, file = paste0(current_path, "loc6nt_UNIFORM_7000.Rda"))
loc6nt_comp_UNIFORM_7000 <- compare_to_ground_truth(mean_loc6nt, 
                                                    loc6nt_UNIFORM_7000, 
                                                    tau_loc6nt, 
                                                    epsilon_loc6nt)$mean
save(loc6nt_comp_UNIFORM_7000, file = paste0(current_path, 
                                             "loc6nt_comp_UNIFORM_7000.Rda"))

########################################################################

system.time(loc6nt_KLUCB370_5000 <- para_bandit_sim_KLUCB(data = data_list6, 
                                                          rounds = 5000, 
                                                          tau = tau_loc6nt, 
                                                          epsilon = epsilon_loc6nt,
                                                          horizon = 370000,
                                                          H = H_loc6nt))

save(loc6nt_KLUCB370_5000, file = paste0(current_path, 
                                         "loc6nt_KLUCB370_5000.Rda"))
loc6nt_comp_KLUCB370_5000 <- compare_to_ground_truth(mean_loc6nt, 
                                                     loc6nt_KLUCB370_5000, 
                                                     tau_loc6nt, epsilon_loc6nt)$mean
save(loc6nt_comp_KLUCB370_5000, file = paste0(current_path, 
                                              "loc6nt_comp_KLUCB370_5000.Rda"))

########################################################################

error_bound <- function(n, H, K) {
  (K*n+H)*exp(-n/H)
}

plot(1:500000, log(error_bound(1:500000, H_loc6nt, 10)), type = "l",
     ylim = c(-6,12))
abline(h=log(0.2))

data_list6_small <- data_list6[1:500]
rm(data_list6)
gc()

system.time(loc6nt_KLUCB30 <- para_bandit_sim_KLUCB(data = data_list6_small, 
                                                    rounds = 2000, 
                                                    tau = tau_loc6nt, 
                                                    epsilon = epsilon_loc6nt,
                                                    horizon = round(H_loc6nt),
                                                    H = H_loc6nt))

system.time(loc6nt_KLUCB370 <- para_bandit_sim_KLUCB(data = data_list6_small, 
                                                     rounds = 2000, 
                                                     tau = tau_loc6nt, 
                                                     epsilon = epsilon_loc6nt,
                                                     horizon = 370000,
                                                     H = H_loc6nt))

system.time(loc6nt_KLUCB7 <- para_bandit_sim_KLUCB(data = data_list6_small, 
                                                   rounds = 2000, 
                                                   tau = tau_loc6nt, 
                                                   epsilon = epsilon_loc6nt,
                                                   horizon = 7000,
                                                   H = H_loc6nt))

save(loc6nt_KLUCB30, file = paste0(current_path, "/loc6nt_KLUCB30.Rda"))
save(loc6nt_KLUCB370, file = paste0(current_path, "/loc6nt_KLUCB370.Rda"))
save(loc6nt_KLUCB7, file = paste0(current_path, "/loc6nt_KLUCB7.Rda"))

loc6nt_comp_KLUCB30 <- compare_to_ground_truth(mean_loc6nt, loc6nt_KLUCB30, 
                                               tau_loc6nt, 
                                               epsilon_loc6nt)$mean
loc6nt_comp_KLUCB370 <- compare_to_ground_truth(mean_loc6nt, loc6nt_KLUCB370, 
                                                tau_loc6nt, 
                                                epsilon_loc6nt)$mean
loc6nt_comp_KLUCB7 <- compare_to_ground_truth(mean_loc6nt, loc6nt_KLUCB7, 
                                              tau_loc6nt, 
                                              epsilon_loc6nt)$mean
save(loc6nt_comp_KLUCB30, file = paste0(current_path, "/loc6nt_comp_KLUCB30.Rda"))
save(loc6nt_comp_KLUCB370, file = paste0(current_path, "/loc6nt_comp_KLUCB370.Rda"))
save(loc6nt_comp_KLUCB7, file = paste0(current_path, "/loc6nt_comp_KLUCB7.Rda"))
