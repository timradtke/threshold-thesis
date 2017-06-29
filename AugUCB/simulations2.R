script_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/AugUCB/"

# load all required functions and bandit models
tau_001 <- 0.035
epsilon_001 <- 0.005

mean_001 <- c(0.0005, 0.0005, 0.001, 0.001, 
                  0.005, 0.005, 0.01, 0.01,
                  0.05, 0.05, 0.1, 0.1,
                  0.11, 0.12, 0.13, 0.14,
                  0.06, 0.07, 0.08, 0.08)
                  
length(mean_001)

####################################################################
# for each model, run 1000 repititions with a budget of 600 samples
source(paste0(script_path, "parallelized_simulations.R"))
####################################################################
source(paste0(script_path, "Uniform.R"))
# Standard uniform sampling (which doesn't need info about tau and epsilon!)
system.time(
  ber001_uniform <- para_bandit_sim_uniform(reps = 1000, means = mean_001,
                                            variances = NA, K = 20, rounds = 1000)
)
#  user  system elapsed 
# 0.889   0.275  36.618
save(ber001_uniform, file = paste0(script_path, "ber001_uniform.Rda"))
ber001_uniform_comp <- compare_to_ground_truth(mean_001, ber001_uniform, 
                                               tau = tau_001, epsilon = epsilon_001)
save(ber001_uniform_comp, file = paste0(script_path, "ber001_uniform_comp.Rda"))
lines(log(ber001_uniform_comp$mean), col = "black", xlim = c(0,1000))
####################################################################

source(paste0(script_path, "APT.R"))
# Standard APT
system.time(
  ber001_APT <- para_bandit_sim_APT(reps = 1000, means = mean_001, 
                                    variances = NA, K = 20, rounds = 1000,
                                    tau = tau_001, epsilon = epsilon_001)
)
# user  system elapsed 
#1.140   0.587 109.989 
save(ber001_APT, file = paste0(script_path, "ber001_APT.Rda"))
ber001_APT_comp <- compare_to_ground_truth(mean_001, ber001_APT, 
                                           tau = tau_001, epsilon = epsilon_001)
save(ber001_APT_comp, file = paste0(script_path, "ber001_APT_comp.Rda"))
lines(log(ber001_APT_comp$mean), col = "red", xlim = c(0,1000))

####################################################################

source(paste0(script_path, "KL.R"))
# Standard KL 
system.time(
  ber001_KL <- para_bandit_sim_KL(reps = 1000, means = mean_001, 
                                  variances = NA, K = 20, rounds = 1000,
                                  tau = tau_001, epsilon = epsilon_001)
)
# user  system elapsed 
#1.105   0.708 162.445 
save(ber001_KL, file = paste0(script_path, "ber001_KL.Rda"))
ber001_KL_comp <- compare_to_ground_truth(mean_001, ber001_KL, 
                                          tau = tau_001, epsilon = epsilon_001)
save(ber001_KL_comp, file = paste0(script_path, "ber001_KL_comp.Rda"))
lines(log(ber001_KL_comp$mean), col = "violet", xlim = c(0,1000))

####################################################################

source(paste0(script_path, "AugUCB.R"))
# Standard AugUCB which doesn't take an epsilon
system.time(
  ber001_AugUCB <- para_bandit_sim_AugUCB(reps = 1000, means = mean_001, 
                                          variances = NA, K = 20, rounds = 1000,
                                          tau = tau_001)
)
# user  system elapsed 
# 2.160   2.160 829.202
save(ber001_AugUCB, file = paste0(script_path, "ber001_AugUCB.Rda"))
ber001_AugUCB_comp <- compare_to_ground_truth(mean_001, ber001_AugUCB, 
                                              tau = tau_001, epsilon = epsilon_001)
save(ber001_AugUCB_comp, file = paste0(script_path, "ber001_AugUCB_comp.Rda"))
lines(log(ber001_AugUCB_comp$mean), col = "green", xlim = c(0,1000))

####################################################################

source(paste0(script_path, "PI.R"))
# Probability of Improvement with informative prior
system.time(
  ber001_PI_informative_prior <- para_bandit_sim_PI(reps = 1000, means = mean_001, 
                                                variances = NA, 
                                                K = 20, rounds = 1000,
                                                tau = tau_001, epsilon = epsilon_001,
                                                alpha = tau_001, beta = 1-tau_001,
                                                mean_prior = NA,
                                                n_prior = NA)
)

save(ber001_PI_informative_prior, file = paste0(script_path, 
                                                "ber001_PI_informative_prior.Rda"))
ber001_PI_informative_prior_comp <- compare_to_ground_truth(mean_001, 
                                                            ber001_PI_informative_prior, 
                                                            tau = tau_001, 
                                                            epsilon = epsilon_001)
save(ber001_PI_informative_prior_comp, file = paste0(script_path, "ber001_PI_informative_prior_comp.Rda"))
lines(log(ber001_PI_informative_prior_comp$mean), col = "blue", xlim = c(0,1000))

####################################################################

source(paste0(script_path, "TTS.R"))
# Threshold Thompson Sampling with informative prior but no samples
system.time(
  ber001_TTS_informative <- para_bandit_sim_tts(reps = 1000, 
                                                means = mean_001, variances = NA, 
                                                K = 20, rounds = 1000,
                                                tau = tau_001, epsilon = epsilon_001,
                                                mean_prior = NA, n_prior = NA,
                                                alpha = tau_001, beta = 1-tau_001)
)
#user  system elapsed 
#1.269   0.668 132.992 
save(ber001_TTS_informative, file = paste0(script_path, "ber001_TTS_informative.Rda"))
ber001_TTS_informative_comp <- compare_to_ground_truth(mean_001, 
                                                       ber001_TTS_informative, 
                                                       tau = tau_001, 
                                                       epsilon = epsilon_001)
save(ber001_TTS_informative_comp, 
     file = paste0(script_path, "ber001_TTS_informative_comp.Rda"))
lines(log(ber001_TTS_informative_comp$mean), col = "orange", xlim = c(0,1000))

####################################################################

source(paste0(script_path, "BayesUCB.R"))
system.time(
  ber001_BUCB_informative <- para_bandit_sim_bucb(reps = 1000,
                                                  means = mean_001,
                                                  K = 20, rounds = 1000,
                                                  rate = "inverse",
                                                 tau = tau_001, epsilon = epsilon_001,
                                                 alpha = tau_001, beta = 1-tau_001)
)
save(ber001_BUCB_informative, file = paste0(script_path, "ber001_BUCB_informative.Rda"))
ber001_BUCB_informative_comp <- compare_to_ground_truth(mean_001, 
                                                        ber001_BUCB_informative, 
                                                        tau = tau_001, 
                                                        epsilon = epsilon_001)
save(ber001_BUCB_informative_comp, 
     file = paste0(script_path, "ber001_BUCB_informative_comp.Rda"))
lines(log(ber001_BUCB_informative_comp$mean), col = "darkgreen", xlim = c(0,1000))
