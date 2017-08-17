# Post processing of data for arm sequence visualization.
# This is for Simulation 2 on synthetic data as used in thesis.

##########################################################

# load the algorithm results from Experiment 2

library(plyr)
current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/AugUCB/"
source(paste0(current_path, "run_from_data.R"))
source(paste0(current_path, "run_parallelized_from_data.R"))
current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/synthetic_experiments/"

load(paste0(current_path, "loc6nt_BUCB_horizon_7000.Rda"))
loc6nt_BUCB_arm_seq <- data.frame(t(laply(loc6nt_BUCB_horizon_7000, 
                                        function(x) x$arm_sequence)))
save(loc6nt_BUCB_arm_seq, file = "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/Thesis/experiment_postprocessing/loc6nt_BUCB_arm_seq.Rda")
rm(loc6nt_BUCB_7000, loc6nt_BUCB_arm_seq)
gc()

load(paste0(current_path, "loc6nt_APT_7000.Rda"))
loc6nt_APT_arm_seq <- data.frame(t(laply(loc6nt_APT_7000, 
                                       function(x) x$arm_sequence)))
save(loc6nt_APT_arm_seq, file = "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/Thesis/experiment_postprocessing/loc6nt_APT_arm_seq.Rda")
rm(loc6nt_APT_7000, loc6nt_APT_arm_seq)
gc()

load(paste0(current_path, "loc6nt_LR_7000.Rda"))
loc6nt_LR_arm_seq <- data.frame(t(laply(loc6nt_LR_7000, 
                                      function(x) x$arm_sequence)))
save(loc6nt_LR_arm_seq, file = "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/Thesis/experiment_postprocessing/loc6nt_LR_arm_seq.Rda")
rm(loc6nt_LR_7000, loc6nt_LR_arm_seq)
gc()


##########################################################

# Now pipe the new objects into ggplot

library(ggplot)
library(ggjoy)
library(dplyr)
library(tidyr)

current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/synthetic_experiments/"
load(paste0(current_path, "data_amo4_mean_firsthalf.Rda"))

load("/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/Thesis/experiment_postprocessing/loc6nt_LR_arm_seq.Rda")

loc6nt_LR_arm_seq_summary <- loc6nt_LR_arm_seq %>% 
  tbl_df() %>% 
  mutate(index = 1:7000) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::filter(index > 10) %>%
  dplyr::group_by(arm, index) %>%
  dplyr::summarize(n = n()) %>%
  ungroup() %>% 
  spread(arm, n) %>%
  as.data.frame()

names(loc6nt_LR_arm_seq_summary)[-1] <- c(paste0("V0", 1:9), "V10")

loc6nt_LR_ma <- data.frame(algorithm = rep("SLR",6981),
                           index = 1:6981,
                           lapply(loc6nt_LR_arm_seq_summary[,-1], FUN = zoo::rollapply, width = 10, mean)) %>%
  gather(arm, mean, -index, -algorithm)

load("/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/Thesis/experiment_postprocessing/loc6nt_APT_arm_seq.Rda")

loc6nt_APT_arm_seq_summary <- loc6nt_APT_arm_seq %>% 
  tbl_df() %>% 
  mutate(index = 1:7000) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::filter(index > 10) %>%
  dplyr::group_by(arm, index) %>%
  dplyr::summarize(n = n()) %>%
  ungroup() %>% 
  spread(arm, n) %>%
  as.data.frame()

names(loc6nt_APT_arm_seq_summary)[-1] <- c(paste0("V0", 1:9), "V10")

loc6nt_APT_ma <- data.frame(algorithm = rep("APT",6981),
                          index = 1:6981,
                          lapply(loc6nt_APT_arm_seq_summary[,-1], FUN = zoo::rollapply, width = 10, mean)) %>%
  gather(arm, mean, -index, -algorithm)

rbind(loc6nt_APT_ma, loc6nt_LR_ma) %>%
  ggplot(aes(x = index, y = mean, group = arm, color = arm)) +
  geom_line(size = 0.3) +
  facet_grid(.~algorithm) +
  labs(x = "Index", y = "10 Rounds Moving Average",
       title = "Average Number of Pulls Allocated to Each Arm per Round",
       subtitle = "5000 iterations; a value of 1000 means that in 20% of 
the iterations the arm was pulled in that round.")

save(loc6nt_LR_ma, loc6nt_APT_ma, file = "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/Thesis/experiment_postprocessing/loc6nt_ma.Rda")
