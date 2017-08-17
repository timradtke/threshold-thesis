# Post processing of data for arm sequence visualization.
# This is for Experiment 2 on real data as used in thesis.

##########################################################

# load the algorithm results from Experiment 2

library(plyr)
current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/AugUCB/"
source(paste0(current_path, "run_from_data.R"))
source(paste0(current_path, "run_parallelized_from_data.R"))
current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/real_experiments/"

load(paste0(current_path, "amo4_BUCB.Rda"))
amo4_BUCB_arm_seq <- data.frame(t(laply(amo4_BUCB, 
                                       function(x) x$arm_sequence)))
save(amo4_BUCB_arm_seq, file = "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/Thesis/experiment_postprocessing/amo4_BUCB_arm_seq.Rda")
rm(amo4_BUCB, amo4_BUCB_arm_seq)
gc()

load(paste0(current_path, "amo4_APT.Rda"))
amo4_APT_arm_seq <- data.frame(t(laply(amo4_APT, 
                                       function(x) x$arm_sequence)))
save(amo4_APT_arm_seq, file = "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/Thesis/experiment_postprocessing/amo4_APT_arm_seq.Rda")
rm(amo4_APT, amo4_APT_arm_seq)
gc()

load(paste0(current_path, "amo4_AugUCB.Rda"))
amo4_AugUCB_arm_seq <- data.frame(t(laply(amo4_AugUCB, 
                                       function(x) x$arm_sequence)))
save(amo4_AugUCB_arm_seq, file = "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/Thesis/experiment_postprocessing/amo4_AugUCB_arm_seq.Rda")
rm(amo4_AugUCB, amo4_AugUCB_arm_seq)
gc()

load(paste0(current_path, "amo4_LR.Rda"))
amo4_LR_arm_seq <- data.frame(t(laply(amo4_LR, 
                                          function(x) x$arm_sequence)))
save(amo4_LR_arm_seq, file = "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/Thesis/experiment_postprocessing/amo4_LR_arm_seq.Rda")
rm(amo4_LR, amo4_LR_arm_seq)
gc()

load(paste0(current_path, "amo4_EVT.Rda"))
amo4_EVT_arm_seq <- data.frame(t(laply(amo4_EVT, 
                                      function(x) x$arm_sequence)))
save(amo4_EVT_arm_seq, file = "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/Thesis/experiment_postprocessing/amo4_EVT_arm_seq.Rda")
rm(amo4_EVT, amo4_EVT_arm_seq)
gc()

load(paste0(current_path, "amo4sim_BUCB.Rda"))
amo4sim_BUCB_arm_seq <- data.frame(t(laply(amo4sim_BUCB, 
                                        function(x) x$arm_sequence)))
save(amo4sim_BUCB_arm_seq, file = "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/Thesis/experiment_postprocessing/amo4sim_BUCB_arm_seq.Rda")
rm(amo4sim_BUCB, amo4sim_BUCB_arm_seq)
gc()

load(paste0(current_path, "amo4sim_APT.Rda"))
amo4sim_APT_arm_seq <- data.frame(t(laply(amo4sim_APT, 
                                       function(x) x$arm_sequence)))
save(amo4sim_APT_arm_seq, file = "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/Thesis/experiment_postprocessing/amo4sim_APT_arm_seq.Rda")
rm(amo4sim_APT, amo4sim_APT_arm_seq)
gc()

load(paste0(current_path, "amo4sim_AugUCB.Rda"))
amo4sim_AugUCB_arm_seq <- data.frame(t(laply(amo4sim_AugUCB, 
                                          function(x) x$arm_sequence)))
save(amo4sim_AugUCB_arm_seq, file = "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/Thesis/experiment_postprocessing/amo4sim_AugUCB_arm_seq.Rda")
rm(amo4sim_AugUCB, amo4sim_AugUCB_arm_seq)
gc()

load(paste0(current_path, "amo4sim_LR.Rda"))
amo4sim_LR_arm_seq <- data.frame(t(laply(amo4sim_LR, 
                                      function(x) x$arm_sequence)))
save(amo4sim_LR_arm_seq, file = "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/Thesis/experiment_postprocessing/amo4sim_LR_arm_seq.Rda")
rm(amo4sim_LR, amo4sim_LR_arm_seq)
gc()

load(paste0(current_path, "amo4sim_EVT.Rda"))
amo4sim_EVT_arm_seq <- data.frame(t(laply(amo4sim_EVT, 
                                             function(x) x$arm_sequence)))
save(amo4sim_EVT_arm_seq, file = "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/Thesis/experiment_postprocessing/amo4sim_EVT_arm_seq.Rda")
rm(amo4sim_EVT, amo4sim_EVT_arm_seq)
gc()

##########################################################

# Now pipe the new objects into ggplot

library(ggplot)
library(ggjoy)
library(dplyr)
library(tidyr)

load(paste0(current_path, "data_amo4_mean_firsthalf.Rda"))

load("/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/Thesis/experiment_postprocessing/amo4_AugUCB_arm_seq.Rda")
load("/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/Thesis/experiment_postprocessing/amo4sim_AugUCB_arm_seq.Rda")

amo4_AugUCB_arm_seq_summary <- amo4_AugUCB_arm_seq %>% 
  tbl_df() %>% 
  mutate(index = 1:10080) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::filter(index > 10) %>%
  dplyr::group_by(arm, index) %>%
  dplyr::summarize(n = n()) %>%
  ungroup() %>% 
  spread(arm, n) %>%
  as.data.frame()
amo4sim_AugUCB_arm_seq_summary <- amo4sim_AugUCB_arm_seq %>% 
  tbl_df() %>% 
  mutate(index = 1:10080) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::filter(index > 10) %>%
  dplyr::group_by(arm, index) %>%
  dplyr::summarize(n = n()) %>%
  ungroup() %>% 
  spread(arm, n) %>%
  as.data.frame()

names(amo4_AugUCB_arm_seq_summary)[-1] <- names(data_amo4_mean_firsthalf)
names(amo4sim_AugUCB_arm_seq_summary)[-1] <- names(data_amo4_mean_firsthalf)

amo4_AugUCB_ma <- data.frame(algorithm = rep("AugUCB",10061),
                             source = rep("Real", 10061),
                        index = 1:10061,
                        lapply(amo4_AugUCB_arm_seq_summary[,-1], FUN = zoo::rollapply, width = 10, mean)) %>%
  gather(arm, mean, -index, -source, -algorithm)
amo4sim_AugUCB_ma <- data.frame(algorithm = rep("AugUCB",10061),
                                source = rep("Synthetic", 10061),
                             index = 1:10061,
                             lapply(amo4sim_AugUCB_arm_seq_summary[,-1], FUN = zoo::rollapply, width = 10, mean)) %>%
  gather(arm, mean, -index, -source, -algorithm)
  
rbind(amo4_AugUCB_ma, amo4sim_AugUCB_ma) %>%
  ggplot(aes(x = index, y = mean, group = arm, color = arm)) +
  geom_line(size = 0.3) +
  facet_grid(.~source) +
  labs(x = "Index", y = "10 Rounds Moving Average",
       title = "Average Number of Pulls Allocated to Each Arm per Round",
       subtitle = "Augmented-UCB")

##########################################################

# Now pipe the new objects into ggplot

library(ggplot)
library(ggjoy)
library(dplyr)
library(tidyr)

load(paste0(current_path, "data_amo4_mean_firsthalf.Rda"))

load("/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/Thesis/experiment_postprocessing/amo4_LR_arm_seq.Rda")
load("/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/Thesis/experiment_postprocessing/amo4sim_LR_arm_seq.Rda")

amo4_LR_arm_seq_summary <- amo4_LR_arm_seq %>% 
  tbl_df() %>% 
  mutate(index = 1:10080) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::filter(index > 10) %>%
  dplyr::group_by(arm, index) %>%
  dplyr::summarize(n = n()) %>%
  ungroup() %>% 
  spread(arm, n) %>%
  as.data.frame()
amo4sim_LR_arm_seq_summary <- amo4sim_LR_arm_seq %>% 
  tbl_df() %>% 
  mutate(index = 1:10080) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::filter(index > 10) %>%
  dplyr::group_by(arm, index) %>%
  dplyr::summarize(n = n()) %>%
  ungroup() %>% 
  spread(arm, n) %>%
  as.data.frame()

names(amo4_LR_arm_seq_summary)[-1] <- names(data_amo4_mean_firsthalf)
names(amo4sim_LR_arm_seq_summary)[-1] <- names(data_amo4_mean_firsthalf)

amo4_LR_ma <- data.frame(algorithm = rep("SLR",10061),
                         source = rep("Real", 10061),
                         index = 1:10061,
                         lapply(amo4_LR_arm_seq_summary[,-1], FUN = zoo::rollapply, width = 10, mean)) %>%
  gather(arm, mean, -index, -source, -algorithm)
amo4sim_LR_ma <- data.frame(algorithm = rep("SLR",10061),
                            source = rep("Synthetic", 10061),
                            index = 1:10061,
                            lapply(amo4sim_LR_arm_seq_summary[,-1], FUN = zoo::rollapply, width = 10, mean)) %>%
  gather(arm, mean, -index, -source, -algorithm)

rbind(amo4_LR_ma, amo4sim_LR_ma) %>%
  ggplot(aes(x = index, y = mean, group = arm, color = arm)) +
  geom_line(size = 0.3) +
  facet_grid(.~source) +
  labs(x = "Index", y = "10 Rounds Moving Average",
       title = "Average Number of Pulls Allocated to Each Arm per Round",
       subtitle = "Simple Likelihood Ratio")

##########################################################

# Now pipe the new objects into ggplot

library(ggplot)
library(ggjoy)
library(dplyr)
library(tidyr)

load(paste0(current_path, "data_amo4_mean_firsthalf.Rda"))

load("/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/Thesis/experiment_postprocessing/amo4_APT_arm_seq.Rda")
load("/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/Thesis/experiment_postprocessing/amo4sim_APT_arm_seq.Rda")

amo4_APT_arm_seq_summary <- amo4_APT_arm_seq %>% 
  tbl_df() %>% 
  mutate(index = 1:10080) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::filter(index > 10) %>%
  dplyr::group_by(arm, index) %>%
  dplyr::summarize(n = n()) %>%
  ungroup() %>% 
  spread(arm, n) %>%
  as.data.frame()
amo4sim_APT_arm_seq_summary <- amo4sim_APT_arm_seq %>% 
  tbl_df() %>% 
  mutate(index = 1:10080) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::filter(index > 10) %>%
  dplyr::group_by(arm, index) %>%
  dplyr::summarize(n = n()) %>%
  ungroup() %>% 
  spread(arm, n) %>%
  as.data.frame()

names(amo4_APT_arm_seq_summary)[-1] <- names(data_amo4_mean_firsthalf)
names(amo4sim_APT_arm_seq_summary)[-1] <- names(data_amo4_mean_firsthalf)

amo4_APT_ma <- data.frame(algorithm = rep("APT",10061),
                          source = rep("Real", 10061),
                          index = 1:10061,
                          lapply(amo4_APT_arm_seq_summary[,-1], FUN = zoo::rollapply, width = 10, mean)) %>%
  gather(arm, mean, -index, -source, -algorithm)
amo4sim_APT_ma <- data.frame(algorithm = rep("APT",10061),
                             source = rep("Synthetic", 10061),
                             index = 1:10061,
                             lapply(amo4sim_APT_arm_seq_summary[,-1], FUN = zoo::rollapply, width = 10, mean)) %>%
  gather(arm, mean, -index, -source, -algorithm)

rbind(amo4_APT_ma, amo4sim_APT_ma) %>%
  ggplot(aes(x = index, y = mean, group = arm, color = arm)) +
  geom_line(size = 0.3) +
  facet_grid(.~source) +
  labs(x = "Index", y = "10 Rounds Moving Average",
       title = "Average Number of Pulls Allocated to Each Arm per Round",
       subtitle = "Anytime Parameter Free")

##########################################################

# Now pipe the new objects into ggplot

library(ggplot)
library(ggjoy)
library(dplyr)
library(tidyr)

load(paste0(current_path, "data_amo4_mean_firsthalf.Rda"))

load("/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/Thesis/experiment_postprocessing/amo4_BUCB_arm_seq.Rda")
load("/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/Thesis/experiment_postprocessing/amo4sim_BUCB_arm_seq.Rda")

amo4_BUCB_arm_seq_summary <- amo4_BUCB_arm_seq %>% 
  tbl_df() %>% 
  mutate(index = 1:10080) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::filter(index > 10) %>%
  dplyr::group_by(arm, index) %>%
  dplyr::summarize(n = n()) %>%
  ungroup() %>% 
  spread(arm, n) %>%
  as.data.frame()
amo4sim_BUCB_arm_seq_summary <- amo4sim_BUCB_arm_seq %>% 
  tbl_df() %>% 
  mutate(index = 1:10080) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::filter(index > 10) %>%
  dplyr::group_by(arm, index) %>%
  dplyr::summarize(n = n()) %>%
  ungroup() %>% 
  spread(arm, n) %>%
  as.data.frame()

names(amo4_BUCB_arm_seq_summary)[-1] <- names(data_amo4_mean_firsthalf)
names(amo4sim_BUCB_arm_seq_summary)[-1] <- names(data_amo4_mean_firsthalf)

amo4_BUCB_ma <- data.frame(algorithm = rep("BUCB",10061),
                          source = rep("Real", 10061),
                          index = 1:10061,
                          lapply(amo4_BUCB_arm_seq_summary[,-1], FUN = zoo::rollapply, width = 10, mean)) %>%
  gather(arm, mean, -index, -source, -algorithm)
amo4sim_BUCB_ma <- data.frame(algorithm = rep("BUCB",10061),
                             source = rep("Synthetic", 10061),
                             index = 1:10061,
                             lapply(amo4sim_BUCB_arm_seq_summary[,-1], FUN = zoo::rollapply, width = 10, mean)) %>%
  gather(arm, mean, -index, -source, -algorithm)

rbind(amo4_BUCB_ma, amo4sim_BUCB_ma) %>%
  ggplot(aes(x = index, y = mean, group = arm, color = arm)) +
  geom_line(size = 0.3) +
  facet_grid(.~source) +
  labs(x = "Index", y = "10 Rounds Moving Average",
       title = "Average Number of Pulls Allocated to Each Arm per Round",
       subtitle = "Anytime Parameter Free")

##########################################################

# Now pipe the new objects into ggplot

library(ggplot)
library(ggjoy)
library(dplyr)
library(tidyr)

load(paste0(current_path, "data_amo4_mean_firsthalf.Rda"))

load("/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/Thesis/experiment_postprocessing/amo4_EVT_arm_seq.Rda")
load("/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/Thesis/experiment_postprocessing/amo4sim_EVT_arm_seq.Rda")

amo4_EVT_arm_seq_summary <- amo4_EVT_arm_seq %>% 
  tbl_df() %>% 
  mutate(index = 1:10080) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::filter(index > 20) %>%
  dplyr::group_by(arm, index) %>%
  dplyr::summarize(n = n()) %>%
  ungroup() %>% 
  spread(arm, n) %>%
  as.data.frame()
amo4sim_EVT_arm_seq_summary <- amo4sim_EVT_arm_seq %>% 
  tbl_df() %>% 
  mutate(index = 1:10080) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::filter(index > 20) %>%
  dplyr::group_by(arm, index) %>%
  dplyr::summarize(n = n()) %>%
  ungroup() %>% 
  spread(arm, n) %>%
  as.data.frame()

names(amo4_EVT_arm_seq_summary)[-1] <- names(data_amo4_mean_firsthalf)
names(amo4sim_EVT_arm_seq_summary)[-1] <- names(data_amo4_mean_firsthalf)

amo4_EVT_ma <- data.frame(algorithm = rep("EVT",10051),
                           source = rep("Real", 10051),
                           index = 11:10061,
                           lapply(amo4_EVT_arm_seq_summary[,-1], FUN = zoo::rollapply, width = 10, mean)) %>%
  gather(arm, mean, -index, -source, -algorithm)
amo4sim_EVT_ma <- data.frame(algorithm = rep("EVT",10051),
                              source = rep("Synthetic", 10051),
                              index = 11:10061,
                              lapply(amo4sim_EVT_arm_seq_summary[,-1], FUN = zoo::rollapply, width = 10, mean)) %>%
  gather(arm, mean, -index, -source, -algorithm)

rbind(amo4_EVT_ma, amo4sim_EVT_ma) %>%
  ggplot(aes(x = index, y = mean, group = arm, color = arm)) +
  geom_line(size = 0.3) +
  facet_grid(.~source) +
  labs(x = "Index", y = "10 Rounds Moving Average",
       title = "Average Number of Pulls Allocated to Each Arm per Round",
       subtitle = "EVT")

##########################################################

save(amo4_LR_ma, amo4sim_LR_ma,
     amo4_APT_ma, amo4sim_APT_ma,
     amo4_BUCB_ma, amo4sim_BUCB_ma,
     amo4_AugUCB_ma, amo4sim_AugUCB_ma,
     amo4_EVT_ma, amo4sim_EVT_ma,
     file = "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/Thesis/experiment_postprocessing/amo4_arm_seq_ma.Rda")

load(file = "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/Thesis/experiment_postprocessing/amo4_arm_seq_ma.Rda")

rbind(amo4_LR_ma, amo4sim_LR_ma, amo4_APT_ma, amo4sim_APT_ma) %>%
  ggplot(aes(x = index, y = mean, group = arm, color = arm)) +
  geom_line(size = 0.3) +
  facet_grid(source~algorithm) +
  labs(x = "Index", y = "10 Rounds Moving Average",
       title = "Average Number of Pulls Allocated to Each Arm per Round")

rbind(amo4_BUCB_ma, amo4sim_BUCB_ma, amo4_AugUCB_ma, amo4sim_AugUCB_ma) %>%
  ggplot(aes(x = index, y = mean, group = arm, color = arm)) +
  geom_line(size = 0.3) +
  facet_grid(source~algorithm) +
  labs(x = "Index", y = "10 Rounds Moving Average",
       title = "Average Number of Pulls Allocated to Each Arm per Round")


rbind(amo4_APT_ma, amo4sim_APT_ma, amo4_LR_ma, amo4sim_LR_ma, amo4_EVT_ma, amo4sim_EVT_ma, amo4_BUCB_ma, amo4sim_BUCB_ma) %>%
  ggplot(aes(x = index, y = mean, group = arm, color = arm)) +
  geom_line(size = 0.3) +
  facet_grid(source~algorithm) +
  labs(x = "Index", y = "10 Rounds Moving Average",
       title = "Share of Iterations Pulling an Arm at a Given Round") +
  theme(legend.position = "bottom") +
  scale_color_discrete(name = "Arms")
