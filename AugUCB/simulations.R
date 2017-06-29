script_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/AugUCB/"

# load all required functions and bandit models

source(paste0(script_path, "Uniform.R"))
source(paste0(script_path, "APT.R"))
source(paste0(script_path, "KL.R"))
source(paste0(script_path, "AugUCB.R"))
source(paste0(script_path, "PI.R"))
source(paste0(script_path, "TTS.R"))
source(paste0(script_path, "JustCount.R"))
source(paste0(script_path, "BayesUCB.R"))
source(paste0(script_path, "BayesLUCB.R"))
source(paste0(script_path, "parallelized_simulations.R"))

# specify the problem that we want to look at

# 20 arms, low threshold of 0.04, epsilon of 0.01
# 5 arms very low, 5 arms very high
# 3 arms below 0.04, 3 arms above 0.04
# 2 arms below within epsilon, 2 arms above within epsilon

mean_twenty <- c(0.0005, 0.0005, 0.001, 0.001, 0.005,
                 0.3, 0.4, 0.5, 0.6, 0.7,
                 0.01, 0.02, 0.02,
                 0.06, 0.06, 0.07,
                 0.035, 0.035, 0.045, 0.045)
length(mean_twenty)

# for each model, run 1000 repititions with a budget of 5000 samples

# Standard uniform sampling (which doesn't need info about tau and epsilon!)
system.time(
  ber20_uniform <- para_bandit_sim_uniform(reps = 1000, means = mean_twenty,
                                           variances = NA, K = 20, rounds = 5000)
)
# user   system  elapsed 
#5.128    4.020 1013.515 

save(ber20_uniform, file = paste0(script_path, "ber20_uniform.Rda"))
ber20_uniform_comp <- compare_to_ground_truth(mean_twenty, ber20_uniform, 
                                              tau = 0.04, epsilon = 0.01)
save(ber20_uniform_comp, file = paste0(script_path, "ber20_uniform_comp.Rda"))
plot(log(ber20_uniform_comp$mean), type = "l", xlim = c(0,1000))

# Standard APT
system.time(
  ber20_APT <- para_bandit_sim_APT(reps = 1000, means = mean_twenty, 
                                   variances = NA, K = 20, rounds = 5000,
                                   tau = 0.04, epsilon = 0.01)
)
#  user   system  elapsed 
# 4.055    4.255 1169.258 
save(ber20_APT, file = paste0(script_path, "ber20_APT.Rda"))
ber20_APT_comp <- compare_to_ground_truth(mean_twenty, ber20_APT, 
                                          tau = 0.04, epsilon = 0.01)
save(ber20_APT_comp, file = paste0(script_path, "ber20_APT_comp.Rda"))
lines(log(ber20_APT_comp$mean), col = "red", xlim = c(0,1000))

# Standard KL 
system.time(
  ber20_KL <- para_bandit_sim_KL(reps = 1000, means = mean_twenty, 
                                   variances = NA, K = 20, rounds = 5000,
                                   tau = 0.04, epsilon = 0.01)
)

save(ber20_KL, file = paste0(script_path, "ber20_KL.Rda"))
ber20_KL_comp <- compare_to_ground_truth(mean_twenty, ber20_KL, 
                                          tau = 0.04, epsilon = 0.01)
save(ber20_KL_comp, file = paste0(script_path, "ber20_KL_comp.Rda"))
lines(log(ber20_KL_comp$mean), col = "violet", xlim = c(0,1000))

# Standard AugUCB which doesn't take an epsilon
system.time(
  ber20_AugUCB <- para_bandit_sim_AugUCB(reps = 1000, means = mean_twenty, 
                                         variances = NA, K = 20, rounds = 5000,
                                         tau = 0.04)
)
#  user   system  elapsed 
# 9.397   12.624 4947.836
save(ber20_AugUCB, file = paste0(script_path, "ber20_AugUCB.Rda"))
ber20_AugUCB_comp <- compare_to_ground_truth(mean_twenty, ber20_AugUCB, 
                                             tau = 0.04, epsilon = 0.01)
save(ber20_AugUCB_comp, file = paste0(script_path, "ber20_AugUCB_comp.Rda"))
lines(log(ber20_AugUCB_comp$mean), col = "green", xlim = c(0,1000))

# Probability of Improvement with (1,1) prior
system.time(
  ber20_PI_uniform_prior <- para_bandit_sim_PI(reps = 1000, means = mean_twenty, 
                                              variances = NA, 
                                              K = 20, rounds = 5000,
                                              tau = 0.04, epsilon = 0.01,
                                              alpha = 1, beta = 1,
                                              mean_prior = NA,
                                              n_prior = NA)
)
#  user   system  elapsed 
# 4.009    4.513 1314.920 
save(ber20_PI_uniform_prior, file = paste0(script_path, "ber20_PI_uniform_prior.Rda"))
ber20_PI_uniform_prior_comp <- compare_to_ground_truth(mean_twenty, 
                                                       ber20_PI_uniform_prior, 
                                                       tau = 0.04, epsilon = 0.01)
save(ber20_PI_uniform_prior_comp, file = paste0(script_path, "ber20_PI_uniform_prior_comp.Rda"))
lines(log(ber20_PI_uniform_prior_comp$mean), col = "blue", xlim = c(0,1000))

# Probability of Improvement with prior at tau (0.04*26, 0.96*26)
system.time(
  ber20_PI_informative <- para_bandit_sim_PI(reps = 1000, means = mean_twenty, 
                                             variances = NA, 
                                             K = 20, rounds = 5000,
                                             tau = 0.04, epsilon = 0.01,
                                             alpha = 0.04*26, beta = 0.96*26,
                                             mean_prior = NA,
                                             n_prior = NA)
)
#  user   system  elapsed 
# 4.036    4.462 1332.926
save(ber20_PI_informative, file = paste0(script_path, "ber20_PI_informative.Rda"))
ber20_PI_informative_comp <- compare_to_ground_truth(mean_twenty, 
                                                     ber20_PI_informative, 
                                                     tau = 0.04, epsilon = 0.01)
save(ber20_PI_informative_comp, file = paste0(script_path, "ber20_PI_informative_comp.Rda"))
lines(log(ber20_PI_informative_comp$mean), col = "blue", xlim = c(0,1000))

# Probability of Improvement with prior at tau (0.04, 0.96)
system.time(
  ber20_PI_informative_nosamples <- para_bandit_sim_PI(reps = 1000, means = mean_twenty, 
                                                       variances = NA, 
                                                       K = 20, rounds = 5000,
                                                       tau = 0.04, epsilon = 0.01,
                                                       alpha = 0.04, beta = 0.96)
)
#  user   system  elapsed 
# 4.783    5.328 1564.269
save(ber20_PI_informative_nosamples, file = paste0(script_path, "ber20_PI_informative_nosamples.Rda"))
ber20_PI_informative_nosamples_comp <- compare_to_ground_truth(mean_twenty, 
                                                               ber20_PI_informative_nosamples, 
                                                               tau = 0.04, epsilon = 0.01)
save(ber20_PI_informative_nosamples_comp, file = paste0(script_path, "ber20_PI_informative_nosamples_comp.Rda"))
lines(log(ber20_PI_informative_nosamples_comp$mean), col = "pink", xlim = c(0,1000))

# Threshold Thompson Sampling with uniform prior
system.time(
  ber20_TTS <- para_bandit_sim_tts(reps = 1000, means = mean_twenty, variances = NA, 
                                   K = 20, rounds = 5000,
                                   tau = 0.04, epsilon = 0.01,
                                   mean_prior = NA, n_prior = NA,
                                   alpha = 1, beta = 1)
)
#  user   system  elapsed 
# 4.708    5.439 1494.472 
save(ber20_TTS, file = paste0(script_path, "ber20_TTS.Rda"))
ber20_TTS_comp <- compare_to_ground_truth(mean_twenty, 
                                          ber20_TTS, 
                                          tau = 0.04, epsilon = 0.01)
save(ber20_TTS_comp, file = paste0(script_path, "ber20_TTS_comp.Rda"))
lines(log(ber20_TTS_comp$mean), col = "orange", xlim = c(0,1000))

# Threshold Thompson Sampling with informative prior but no samples
system.time(
  ber20_TTS_informative <- para_bandit_sim_tts(reps = 1000, 
                                               means = mean_twenty, variances = NA, 
                                               K = 20, rounds = 5000,
                                               tau = 0.04, epsilon = 0.01,
                                               mean_prior = NA, n_prior = NA,
                                               alpha = 0.04, beta = 0.96)
)
#  user   system  elapsed 
# 4.708    5.439 1494.472 
save(ber20_TTS_informative, file = paste0(script_path, "ber20_TTS.Rda"))
ber20_TTS_informative_comp <- compare_to_ground_truth(mean_twenty, 
                                                      ber20_TTS_informative, 
                                          tau = 0.04, epsilon = 0.01)
save(ber20_TTS_informative_comp, 
     file = paste0(script_path, "ber20_TTS_informative_comp.Rda"))
lines(log(ber20_TTS_informative_comp$mean), col = "orange", xlim = c(0,1000))



# Just Count Bandit with informative Prior
system.time(
  ber20_JC_informative <- para_bandit_sim_tts(reps = 1000, 
                                               means = mean_twenty, 
                                               K = 20, rounds = 2000,
                                               tau = 0.04, epsilon = 0.01,
                                               alpha = 0.04, beta = 0.96)
)
#  user   system  elapsed 
# 4.708    5.439 1494.472 
save(ber20_JC_informative, file = paste0(script_path, "ber20_JC_informative.Rda"))
ber20_JC_informative_comp <- compare_to_ground_truth(mean_twenty, 
                                                     ber20_JC_informative, 
                                                     tau = 0.04, epsilon = 0.01)
save(ber20_JC_informative_comp, 
     file = paste0(script_path, "ber20_JC_informative_comp.Rda"))



# Bayes UCB with informative prior
system.time(
  ber20_BUCB_informative <- para_bandit_sim_bucb(reps = 1000, 
                                                 means = mean_twenty, 
                                                 K = 20, rounds = 2000,
                                                 tau = 0.04, epsilon = 0.01,
                                                 alpha = 0.04, beta = 0.96)
)
save(ber20_BUCB_informative, file = paste0(script_path, "ber20_BUCB_informative.Rda"))
ber20_BUCB_informative_comp <- compare_to_ground_truth(mean_twenty, 
                                                       ber20_BUCB_informative, 
                                                       tau = 0.04, epsilon = 0.01)
save(ber20_BUCB_informative_comp, 
     file = paste0(script_path, "ber20_BUCB_informative_comp.Rda"))


# Bayes UCB with uninformative prior
system.time(
  ber20_BUCB <- para_bandit_sim_bucb(reps = 1000, 
                                     means = mean_twenty, 
                                     K = 20, rounds = 2000,
                                     tau = 0.04, epsilon = 0.01,
                                     alpha = 1, beta = 1)
)
#  user   system  elapsed 
# 4.708    5.439 1494.472 
save(ber20_BUCB, file = paste0(script_path, "ber20_BUCB.Rda"))
ber20_BUCB_comp <- compare_to_ground_truth(mean_twenty, 
                                           ber20_BUCB, 
                                           tau = 0.04, epsilon = 0.01)
save(ber20_BUCB_comp, 
     file = paste0(script_path, "ber20_BUCB_comp.Rda"))

# Get Bayes UCB with different rate for upper bounds, lower prior
system.time(
  ber20_BUCB_inverse_lowPrior <- para_bandit_sim_bucb(reps = 1000, 
                                                     means = mean_twenty, 
                                                     K = 20, rounds = 2000,
                                                     tau = 0.04, epsilon = 0.01,
                                                     alpha = 0.03, beta = 0.97,
                                                     rate = "inverse")
)

save(ber20_BUCB_inverse_lowPrior, file = paste0(script_path, 
                                                "ber20_BUCB_inverse_lowPrior.Rda"))
ber20_BUCB_inverse_lowPrior_comp <- compare_to_ground_truth(mean_twenty, 
                                                            ber20_BUCB_inverse_lowPrior, 
                                                            tau = 0.04, epsilon = 0.01)
save(ber20_BUCB_inverse_lowPrior_comp, 
     file = paste0(script_path, "ber20_BUCB_inverse_lowPrior_comp.Rda"))

# Get Bayes UCB with different rate for upper bounds
system.time(
  ber20_BUCB_inverse_horizon <- para_bandit_sim_bucb(reps = 1000, 
                                     means = mean_twenty, 
                                     K = 20, rounds = 2000,
                                     tau = 0.04, epsilon = 0.01,
                                     alpha = 0.04, beta = 0.96,
                                     rate = "inverse_horizon")
)
#  user   system  elapsed 
# 4.708    5.439 1494.472 
save(ber20_BUCB_inverse_horizon, file = paste0(script_path, "ber20_BUCB_inverse_horizon.Rda"))
ber20_BUCB_inverse_horizon_comp <- compare_to_ground_truth(mean_twenty, 
                                                           ber20_BUCB_inverse_horizon, 
                                           tau = 0.04, epsilon = 0.01)
save(ber20_BUCB_inverse_horizon_comp, 
     file = paste0(script_path, "ber20_BUCB_inverse_horizon_comp.Rda"))



# Bayes LUCB with informative prior
system.time(
  ber20_BLUCB <- para_bandit_sim_blucb(reps = 1000,
                                       means = mean_twenty, 
                                       K = 20, rounds = 2000,
                                       tau = 0.04, epsilon = 0.01,
                                       alpha = 0.04, beta = 0.96)
)
#  user   system  elapsed 
# 4.708    5.439 1494.472 
save(ber20_BLUCB, file = paste0(script_path, "ber20_BLUCB.Rda"))
ber20_BLUCB_comp <- compare_to_ground_truth(mean_twenty, 
                                            ber20_BLUCB, 
                                            tau = 0.04, epsilon = 0.01)
save(ber20_BLUCB_comp, 
     file = paste0(script_path, "ber20_BLUCB_comp.Rda"))
