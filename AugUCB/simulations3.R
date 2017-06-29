script_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/AugUCB/"

# load all required functions and bandit models
tau_002 <- 0.035
epsilon_002 <- 0.005

# symmetric distances
mean_002 <- c(0.033, 0.033, 0.037, 0.037,
              0.031, 0.031, 0.039, 0.039,
              0.027, 0.027, 0.043, 0.043,
              0.02, 0.02, 0.05, 0.05,
              0.005, 0.005, 0.065, 0.065)

length(mean_002)

####################################################################
# for each model, run 1000 repititions with a budget of 600 samples
source(paste0(script_path, "parallelized_simulations.R"))
####################################################################
source(paste0(script_path, "Uniform.R"))
# Standard uniform sampling (which doesn't need info about tau and epsilon!)
system.time(
  ber002_uniform <- para_bandit_sim_uniform(reps = 1000, means = mean_002,
                                            variances = NA, K = 20, rounds = 1000)
)
save(ber002_uniform, file = paste0(script_path, "ber002_uniform.Rda"))
ber002_uniform_comp <- compare_to_ground_truth(mean_002, ber002_uniform, 
                                               tau = tau_002, epsilon = epsilon_002)
save(ber002_uniform_comp, file = paste0(script_path, "ber002_uniform_comp.Rda"))
plot(log(ber002_uniform_comp$mean), type = "l", col = "black", xlim = c(0,1000), ylim = c(-1,0))
####################################################################

source(paste0(script_path, "APT.R"))
# Standard APT
system.time(
  ber002_APT <- para_bandit_sim_APT(reps = 1000, means = mean_002, 
                                    variances = NA, K = 20, rounds = 1000,
                                    tau = tau_002, epsilon = epsilon_002)
)
save(ber002_APT, file = paste0(script_path, "ber002_APT.Rda"))
ber002_APT_comp <- compare_to_ground_truth(mean_002, ber002_APT, 
                                           tau = tau_002, epsilon = epsilon_002)
save(ber002_APT_comp, file = paste0(script_path, "ber002_APT_comp.Rda"))
lines(log(ber002_APT_comp$mean), col = "red", xlim = c(0,1000))

####################################################################

source(paste0(script_path, "KL.R"))
# Standard KL 
system.time(
  ber002_KL <- para_bandit_sim_KL(reps = 1000, means = mean_002, 
                                  variances = NA, K = 20, rounds = 1000,
                                  tau = tau_002, epsilon = epsilon_002)
)
save(ber002_KL, file = paste0(script_path, "ber002_KL.Rda"))
ber002_KL_comp <- compare_to_ground_truth(mean_002, ber002_KL, 
                                          tau = tau_002, epsilon = epsilon_002)
save(ber002_KL_comp, file = paste0(script_path, "ber002_KL_comp.Rda"))
lines(log(ber002_KL_comp$mean), col = "violet", xlim = c(0,1000))

####################################################################

source(paste0(script_path, "KL.R"))
# KL not taking into account the epsilon
system.time(
  ber002_KL_at_tau <- para_bandit_sim_KL(reps = 1000, means = mean_002, 
                                  variances = NA, K = 20, rounds = 1000,
                                  tau = tau_002, epsilon = epsilon_002,
                                  at_tau = TRUE)
)
save(ber002_KL_at_tau, file = paste0(script_path, "ber002_KL_at_tau.Rda"))
ber002_KL_at_tau_comp <- compare_to_ground_truth(mean_002, ber002_KL_at_tau, 
                                          tau = tau_002, epsilon = epsilon_002)
save(ber002_KL_at_tau_comp, file = paste0(script_path, "ber002_KL_at_tau_comp.Rda"))
lines(log(ber002_KL_at_tau_comp$mean), col = "yellow", xlim = c(0,1000))

####################################################################

source(paste0(script_path, "AugUCB.R"))
# Standard AugUCB which doesn't take an epsilon
system.time(
  ber002_AugUCB <- para_bandit_sim_AugUCB(reps = 1000, means = mean_002, 
                                          variances = NA, K = 20, rounds = 1000,
                                          tau = tau_002)
)
save(ber002_AugUCB, file = paste0(script_path, "ber002_AugUCB.Rda"))
ber002_AugUCB_comp <- compare_to_ground_truth(mean_002, ber002_AugUCB, 
                                              tau = tau_002, epsilon = epsilon_002)
save(ber002_AugUCB_comp, file = paste0(script_path, "ber002_AugUCB_comp.Rda"))
lines(log(ber002_AugUCB_comp$mean), col = "green", xlim = c(0,1000))

####################################################################

source(paste0(script_path, "PI.R"))
# Probability of Improvement with informative prior
system.time(
  ber002_PI_informative_prior <- para_bandit_sim_PI(reps = 1000, means = mean_002, 
                                                    variances = NA, 
                                                    K = 20, rounds = 1000,
                                                    tau = tau_002, epsilon = epsilon_002,
                                                    alpha = tau_002, beta = 1-tau_002,
                                                    mean_prior = NA,
                                                    n_prior = NA)
)

save(ber002_PI_informative_prior, file = paste0(script_path, 
                                                "ber002_PI_informative_prior.Rda"))
ber002_PI_informative_prior_comp <- compare_to_ground_truth(mean_002, 
                                                            ber002_PI_informative_prior, 
                                                            tau = tau_002, 
                                                            epsilon = epsilon_002)
save(ber002_PI_informative_prior_comp, file = paste0(script_path, "ber002_PI_informative_prior_comp.Rda"))
lines(log(ber002_PI_informative_prior_comp$mean), col = "blue", xlim = c(0,1000))

####################################################################

source(paste0(script_path, "TTS.R"))
# Threshold Thompson Sampling with informative prior but no samples
system.time(
  ber002_TTS_informative <- para_bandit_sim_tts(reps = 1000, 
                                                means = mean_002, variances = NA, 
                                                K = 20, rounds = 1000,
                                                tau = tau_002, epsilon = epsilon_002,
                                                mean_prior = NA, n_prior = NA,
                                                alpha = tau_002, beta = 1-tau_002)
)
save(ber002_TTS_informative, file = paste0(script_path, "ber002_TTS_informative.Rda"))
ber002_TTS_informative_comp <- compare_to_ground_truth(mean_002, 
                                                       ber002_TTS_informative, 
                                                       tau = tau_002, 
                                                       epsilon = epsilon_002)
save(ber002_TTS_informative_comp, 
     file = paste0(script_path, "ber002_TTS_informative_comp.Rda"))
lines(log(ber002_TTS_informative_comp$mean), col = "orange", xlim = c(0,1000))

####################################################################

system.time(
  ber002_TTS_informative_samples <- para_bandit_sim_tts(reps = 1000, 
                                                means = mean_002, variances = NA, 
                                                K = 20, rounds = 1000,
                                                tau = tau_002, epsilon = epsilon_002,
                                                mean_prior = NA, n_prior = NA,
                                                alpha = tau_002*10, beta = (1-tau_002)*10)
)
save(ber002_TTS_informative_samples, file = paste0(script_path, "ber002_TTS_informative_samples.Rda"))
ber002_TTS_informative_samples_comp <- compare_to_ground_truth(mean_002, 
                                                       ber002_TTS_informative_samples, 
                                                       tau = tau_002, 
                                                       epsilon = epsilon_002)
save(ber002_TTS_informative_samples_comp, 
     file = paste0(script_path, "ber002_TTS_informative_samples_comp.Rda"))
lines(log(ber002_TTS_informative_samples_comp$mean), col = "orange", xlim = c(0,1000))

####################################################################

system.time(
  ber002_TTS_informative_moresamples <- para_bandit_sim_tts(reps = 1000, 
                                                        means = mean_002, variances = NA, 
                                                        K = 20, rounds = 1000,
                                                        tau = tau_002, epsilon = epsilon_002,
                                                        mean_prior = NA, n_prior = NA,
                                                        alpha = tau_002*25, beta = (1-tau_002)*25)
)
save(ber002_TTS_informative_moresamples, file = paste0(script_path, "ber002_TTS_informative_moresamples.Rda"))
ber002_TTS_informative_moresamples_comp <- compare_to_ground_truth(mean_002, 
                                                                   ber002_TTS_informative_moresamples, 
                                                               tau = tau_002, 
                                                               epsilon = epsilon_002)
save(ber002_TTS_informative_moresamples_comp, 
     file = paste0(script_path, "ber002_TTS_informative_moresamples_comp.Rda"))
lines(log(ber002_TTS_informative_moresamples_comp$mean), col = "orange", xlim = c(0,1000))

####################################################################

source(paste0(script_path, "BayesUCB.R"))
system.time(
  ber002_BUCB_informative <- para_bandit_sim_bucb(reps = 1000,
                                                  means = mean_002,
                                                  K = 20, rounds = 1000,
                                                  rate = "inverse",
                                                  tau = tau_002, epsilon = epsilon_002,
                                                  alpha = tau_002, beta = 1-tau_002)
)
save(ber002_BUCB_informative, file = paste0(script_path, "ber002_BUCB_informative.Rda"))
ber002_BUCB_informative_comp <- compare_to_ground_truth(mean_002, 
                                                        ber002_BUCB_informative, 
                                                        tau = tau_002, 
                                                        epsilon = epsilon_002)
save(ber002_BUCB_informative_comp, 
     file = paste0(script_path, "ber002_BUCB_informative_comp.Rda"))
lines(log(ber002_BUCB_informative_comp$mean), col = "darkgreen", xlim = c(0,1000))
