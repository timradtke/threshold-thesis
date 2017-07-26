# a script that takes experiment results and transforms them into
# data frames for piping into ggplot

library(plyr)
current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/AugUCB/"
source(paste0(current_path, "run_from_data.R"))
source(paste0(current_path, "run_parallelized_from_data.R"))
current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/synthetic_experiments/"

######################################################################

#amo1

load(paste0("/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/real_experiments/amo1_APT.Rda"))
load(paste0("/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/real_experiments/amo1_BUCB.Rda"))

arm_seq_APT_amo1 <- data.frame(t(laply(amo1_APT, 
                                       function(x) x$arm_sequence)))
arm_seq_BUCB_amo1 <- data.frame(t(laply(amo1_BUCB, 
                                       function(x) x$arm_sequence)))

save(arm_seq_APT_amo1,
     file = "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/real_experiments/arm_seq_APT_amo1.Rda")
save(arm_seq_BUCB_amo1,
     file = "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/real_experiments/arm_seq_BUCB_amo1.Rda")

library(ggplot)
library(ggjoy)
library(dplyr)
library(tidyr)

arm_seq_APT_amo1 %>% tbl_df() %>% 
  mutate(index = 1:7580) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::filter(index > 10) %>%
  dplyr::group_by(arm, index) %>%
  dplyr::summarize(n = n()) %>%
  ungroup() %>%
  mutate(arm = as.factor(arm)) %>%
  ggplot(aes(x = index, y = arm, height = n)) +
  geom_joy(stat = "identity", scale = 1.2, size = 0.2) + theme_joy() +
  labs(title = "APT", y = "Arm", x = "Round",
       subtitle = "Pulls per arm and round after 5000 simulations")

arm_seq_APT_amo1 %>% tbl_df() %>% 
  mutate(index = 1:7580) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::filter(index > 10) %>%
  dplyr::group_by(arm, index) %>%
  dplyr::summarize(n = n()) %>%
  ungroup() %>%
  mutate(arm = as.factor(arm)) %>%
  ggplot(aes(x = index, y = n, group = arm, color = arm)) +
  geom_line(size = 0.1)

arm_seq_BUCB_amo1 %>% tbl_df() %>% 
  mutate(index = 1:7580) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::filter(index > 10) %>%
  dplyr::group_by(arm, index) %>%
  dplyr::summarize(n = n()) %>%
  ungroup() %>%
  mutate(arm = as.factor(arm)) %>%
  ggplot(aes(x = index, y = arm, height = n)) +
  geom_joy(stat = "identity", scale = 1.2, size = 0.2) + theme_joy() +
  labs(title = "Bayes-UCB", y = "Arm", x = "Round",
       subtitle = "Pulls per arm and round after 5000 simulations")
######################################################################

# loc2

mean_loc2 <- c(0.4-0.2^(1:4),
               0.45, 0.55,
               0.6+0.2^(5-(1:4)))
tau_loc2 <- 0.5
epsilon_loc2 <- 0.1

load(file = paste0(current_path, "loc2_BUCB_horizon_4000.Rda"))
load(file = paste0(current_path, "loc2_APT_2000.Rda"))
load(file = paste0(current_path, "loc2_UNIFORM_2000.Rda"))
load(file = paste0(current_path, "loc2_AugUCB_2000.Rda"))
load(file = paste0(current_path, "loc2_KLUCB80_long.Rda"))

arm_seq_APT_loc2 <- data.frame(t(laply(loc2_APT_2000, 
                                  function(x) x$arm_sequence)))
arm_seq_BUCB_loc2 <- data.frame(t(laply(loc2_BUCB_horizon_4000, 
                                   function(x) x$arm_sequence)))
arm_seq_AugUCB_loc2 <- data.frame(t(laply(loc2_AugUCB_2000, 
                                     function(x) x$arm_sequence)))
arm_seq_KLUCB_loc2 <- data.frame(t(laply(loc2_KLUCB80_long, 
                                    function(x) x$arm_sequence)))

save(arm_seq_APT_loc2, arm_seq_BUCB_loc2, 
     arm_seq_AugUCB_loc2, arm_seq_KLUCB_loc2,
     file = paste0(current_path, "arm_seq_loc2.Rda"))
rm(arm_seq_APT_loc2, arm_seq_BUCB_loc2, 
   arm_seq_AugUCB_loc2, arm_seq_KLUCB_loc2)
gc()

wrong_arms_UNIFORM_loc2 <- get_wrong_arms_per_round(mean_loc2, 
                                                    loc2_UNIFORM_2000,
                                                    tau_loc2, epsilon_loc2)
wrong_arms_APT_loc2 <- get_wrong_arms_per_round(mean_loc2, loc2_APT_2000,
                                                tau_loc2, epsilon_loc2)
wrong_arms_AugUCB_loc2 <- get_wrong_arms_per_round(mean_loc2, loc2_AugUCB_2000,
                                                tau_loc2, epsilon_loc2)
wrong_arms_BUCB_loc2 <- get_wrong_arms_per_round(mean_loc2, loc2_BUCB_horizon_4000,
                                                tau_loc2, epsilon_loc2)
wrong_arms_KLUCB_loc2 <- get_wrong_arms_per_round(mean_loc2, loc2_KLUCB80_long,
                                                tau_loc2, epsilon_loc2)

save(wrong_arms_BUCB_loc2,
     file = paste0(current_path, "wrong_arms_BUCB_loc2.Rda"))
save(wrong_arms_UNIFORM_loc2, wrong_arms_BUCB_loc2,
     wrong_arms_APT_loc2, wrong_arms_AugUCB_loc2, wrong_arms_KLUCB_loc2,
     file = paste0(current_path, "wrong_arms_loc2.Rda"))
rm(wrong_arms_UNIFORM_loc2, wrong_arms_BUCB_loc2,
   wrong_arms_APT_loc2, wrong_arms_AugUCB_loc2, wrong_arms_KLUCB_loc2,
   loc2_UNIFORM_2000, loc2_APT_2000, loc2_comp_AugUCB_2000,
   loc2_BUCB_horizon_4000, loc2_KLUCB80_long)
gc()

######################################################################

# loc3

mean_loc3 <- c(0.001, 0.005, 0.01, 0.015,
               0.0475, 0.0525,
               0.085, 0.09, 0.095, 0.099)
tau_loc3 <- 0.05
epsilon_loc3 <- 0.005

load(file = paste0(current_path, "loc3_BUCB_horizon_7000.Rda"))
load(file = paste0(current_path, "loc3_APT_7000.Rda"))
load(file = paste0(current_path, "loc3_UNIFORM_7000.Rda"))
load(file = paste0(current_path, "loc3_AugUCB_7000.Rda"))
load(file = paste0(current_path, "loc3_KLUCB370_5000.Rda"))

arm_seq_APT_loc3 <- data.frame(t(laply(loc3_APT_7000, 
                                       function(x) x$arm_sequence)))
arm_seq_BUCB_loc3 <- data.frame(t(laply(loc3_BUCB_horizon_7000, 
                                        function(x) x$arm_sequence)))
arm_seq_AugUCB_loc3 <- data.frame(t(laply(loc3_AugUCB_7000, 
                                          function(x) x$arm_sequence)))
arm_seq_KLUCB_loc3 <- data.frame(t(laply(loc3_KLUCB370_5000, 
                                         function(x) x$arm_sequence)))

save(arm_seq_APT_loc3, arm_seq_BUCB_loc3, 
     arm_seq_AugUCB_loc3, arm_seq_KLUCB_loc3,
     file = paste0(current_path, "arm_seq_loc3.Rda"))
rm(arm_seq_APT_loc3, arm_seq_BUCB_loc3, 
   arm_seq_AugUCB_loc3, arm_seq_KLUCB_loc3)
gc()

wrong_arms_UNIFORM_loc3 <- get_wrong_arms_per_round(mean_loc3, 
                                                    loc3_UNIFORM_7000,
                                                    tau_loc3, epsilon_loc3)
wrong_arms_APT_loc3 <- get_wrong_arms_per_round(mean_loc3, loc3_APT_7000,
                                                tau_loc3, epsilon_loc3)
wrong_arms_AugUCB_loc3 <- get_wrong_arms_per_round(mean_loc3, loc3_AugUCB_7000,
                                                   tau_loc3, epsilon_loc3)
wrong_arms_BUCB_loc3 <- get_wrong_arms_per_round(mean_loc3, loc3_BUCB_horizon_7000,
                                                 tau_loc3, epsilon_loc3)
wrong_arms_KLUCB_loc3 <- get_wrong_arms_per_round(mean_loc3, loc3_KLUCB370_5000,
                                                  tau_loc3, epsilon_loc3)

save(wrong_arms_UNIFORM_loc3, wrong_arms_BUCB_loc3,
     wrong_arms_APT_loc3, wrong_arms_AugUCB_loc3, wrong_arms_KLUCB_loc3,
     file = paste0(current_path, "wrong_arms_loc3.Rda"))
rm(wrong_arms_UNIFORM_loc3, wrong_arms_BUCB_loc3,
   wrong_arms_APT_loc3, wrong_arms_AugUCB_loc3, wrong_arms_KLUCB_loc3,
   loc3_UNIFORM_7000, loc3_APT_7000, loc3_AugUCB_7000, loc3_BUCB_horizon_7000,
   loc3_KLUCB370_5000)
gc()

######################################################################

# loc4

mean_loc4 <- c(0.0005, 0.001, 0.001, 0.005, 0.009,
               0.02, 0.0275, 0.0325, 0.04,
               0.08, 0.09)
tau_loc4 <- 0.03
epsilon_loc4 <- 0.005

load(paste0(current_path, "loc4_BUCB_horizon_8000.Rda"))
load(paste0(current_path, "loc4_AugUCB_8000.Rda"))
load(paste0(current_path, "loc4_APT_8000.Rda"))
load(paste0(current_path, "loc4_KLUCB90_long.Rda"))
load(paste0(current_path, "loc4_UNIFORM_8000.Rda"))

arm_seq_APT_loc4 <- data.frame(t(laply(loc4_APT_8000, 
                                       function(x) x$arm_sequence)))
arm_seq_BUCB_loc4 <- data.frame(t(laply(loc4_BUCB_horizon_8000, 
                                        function(x) x$arm_sequence)))
arm_seq_AugUCB_loc4 <- data.frame(t(laply(loc4_AugUCB_8000, 
                                          function(x) x$arm_sequence)))
arm_seq_KLUCB_loc4 <- data.frame(t(laply(loc4_KLUCB90_long, 
                                         function(x) x$arm_sequence)))

save(arm_seq_APT_loc4, arm_seq_BUCB_loc4, 
     arm_seq_AugUCB_loc4, arm_seq_KLUCB_loc4,
     file = paste0(current_path, "arm_seq_loc4.Rda"))
rm(arm_seq_APT_loc4, arm_seq_BUCB_loc4, 
   arm_seq_AugUCB_loc4, arm_seq_KLUCB_loc4)
gc()

wrong_arms_UNIFORM_loc4 <- get_wrong_arms_per_round(mean_loc4, 
                                                    loc4_UNIFORM_7000,
                                                    tau_loc4, epsilon_loc4)
wrong_arms_APT_loc4 <- get_wrong_arms_per_round(mean_loc4, loc4_APT_7000,
                                                tau_loc4, epsilon_loc4)
wrong_arms_AugUCB_loc4 <- get_wrong_arms_per_round(mean_loc4, loc4_AugUCB_7000,
                                                   tau_loc4, epsilon_loc4)
wrong_arms_BUCB_loc4 <- get_wrong_arms_per_round(mean_loc4, loc4_BUCB_horizon_7000,
                                                 tau_loc4, epsilon_loc4)
wrong_arms_KLUCB_loc4 <- get_wrong_arms_per_round(mean_loc4, loc4_KLUCB470_5000,
                                                  tau_loc4, epsilon_loc4)

save(wrong_arms_UNIFORM_loc4, wrong_arms_BUCB_loc4,
     wrong_arms_APT_loc4, wrong_arms_AugUCB_loc4, wrong_arms_KLUCB_loc4,
     file = paste0(current_path, "wrong_arms_loc4.Rda"))
rm(wrong_arms_UNIFORM_loc4, wrong_arms_BUCB_loc4,
   wrong_arms_APT_loc4, wrong_arms_AugUCB_loc4, wrong_arms_KLUCB_loc4,
   loc4_UNIFORM_7000, loc4_APT_7000, loc4_AugUCB_7000, loc4_BUCB_horizon_7000,
   loc4_KLUCB370_5000)
gc()

######################################################################

# loc5

mean_loc5 <- c(0.0005, 0.001,
               0.05, 0.055,
               0.065, 0.075,
               0.085, 0.09,
               0.1, 0.1)
tau_loc5 <- 0.07
epsilon_loc5 <- 0.01

load(paste0(current_path, "loc5_KLUCB72_long.Rda"))
load(paste0(current_path, "loc5_BUCB_horizon_long_8000.Rda"))
load(paste0(current_path, "loc5_AugUCB_8000.Rda"))
load(paste0(current_path, "loc5_APT_8000.Rda"))
load(paste0(current_path, "loc5_UNIFORM_8000.Rda"))

arm_seq_APT_loc5 <- data.frame(t(laply(loc5_APT_8000, 
                                       function(x) x$arm_sequence)))
arm_seq_BUCB_loc5 <- data.frame(t(laply(loc5_BUCB_horizon_long_8000, 
                                        function(x) x$arm_sequence)))
arm_seq_AugUCB_loc5 <- data.frame(t(laply(loc5_AugUCB_8000, 
                                          function(x) x$arm_sequence)))
arm_seq_KLUCB_loc5 <- data.frame(t(laply(loc5_KLUCB72_long, 
                                         function(x) x$arm_sequence)))

save(arm_seq_APT_loc5, arm_seq_BUCB_loc5, 
     arm_seq_AugUCB_loc5, arm_seq_KLUCB_loc5,
     file = paste0(current_path, "arm_seq_loc5.Rda"))
rm(arm_seq_APT_loc3, arm_seq_BUCB_loc3, 
   arm_seq_AugUCB_loc3, arm_seq_KLUCB_loc3)
gc()

wrong_arms_UNIFORM_loc5 <- get_wrong_arms_per_round(mean_loc5, 
                                                    loc5_UNIFORM_7000,
                                                    tau_loc5, epsilon_loc5)
wrong_arms_APT_loc5 <- get_wrong_arms_per_round(mean_loc5, loc5_APT_7000,
                                                tau_loc5, epsilon_loc5)
wrong_arms_AugUCB_loc5 <- get_wrong_arms_per_round(mean_loc5, loc5_AugUCB_7000,
                                                   tau_loc5, epsilon_loc5)
wrong_arms_BUCB_loc5 <- get_wrong_arms_per_round(mean_loc5, 
                                                 loc5_BUCB_horizon_long_8000,
                                                 tau_loc5, epsilon_loc5)
wrong_arms_KLUCB_loc5 <- get_wrong_arms_per_round(mean_loc5, loc5_KLUCB370_5000,
                                                  tau_loc5, epsilon_loc5)

save(wrong_arms_BUCB_loc5, 
     file = paste0(current_path, "wrong_arms_BUCB_loc5.Rda"))
save(wrong_arms_UNIFORM_loc5, wrong_arms_BUCB_loc5,
     wrong_arms_APT_loc5, wrong_arms_AugUCB_loc5, wrong_arms_KLUCB_loc5,
     file = paste0(current_path, "wrong_arms_loc5.Rda"))
rm(wrong_arms_UNIFORM_loc5, wrong_arms_BUCB_loc5,
   wrong_arms_APT_loc5, wrong_arms_AugUCB_loc5, wrong_arms_KLUCB_loc5,
   loc5_UNIFORM_7000, loc5_APT_7000, loc5_AugUCB_7000, loc5_BUCB_horizon_7000,
   loc5_KLUCB370_5000)
gc()