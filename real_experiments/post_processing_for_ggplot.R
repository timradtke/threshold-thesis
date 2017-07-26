# a script that takes experiment results and transforms them into
# data frames for piping into ggplot

library(plyr)
current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/AugUCB/"
source(paste0(current_path, "run_from_data.R"))
source(paste0(current_path, "run_parallelized_from_data.R"))
current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/real_experiments/"

######################################################################

#amo1

load(paste0(current_path, "amo1_APT.Rda"))
arm_seq_APT_amo1 <- data.frame(t(laply(amo1_APT, 
                                       function(x) x$arm_sequence)))
save(arm_seq_APT_amo1,
     file = "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/real_experiments/arm_seq_APT_amo1.Rda")
rm(amo1_APT)
gc()

load(paste0(current_path, "amo1_BUCB.Rda"))
arm_seq_BUCB_amo1 <- data.frame(t(laply(amo1_BUCB, 
                                        function(x) x$arm_sequence)))
save(arm_seq_BUCB_amo1,
     file = "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/real_experiments/arm_seq_BUCB_amo1.Rda")
rm(amo1_BUCB)
gc()

load(paste0(current_path, "amo1_KLUCB.Rda"))
arm_seq_KLUCB_amo1 <- data.frame(t(laply(amo1_KLUCB, 
                                        function(x) x$arm_sequence)))
save(arm_seq_KLUCB_amo1,
     file = "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/real_experiments/arm_seq_KLUCB_amo1.Rda")
rm(amo1_KLUCB)
gc()

load(paste0(current_path, "amo1_AugUCB.Rda"))
arm_seq_AugUCB_amo1 <- data.frame(t(laply(amo1_AugUCB, 
                                        function(x) x$arm_sequence)))
save(arm_seq_AugUCB_amo1,
     file = "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/real_experiments/arm_seq_AugUCB_amo1.Rda")
rm(amo1_AugUCB)
gc()

load(paste0(current_path, "amo1_UNIFORM.Rda"))
arm_seq_UNIFORM_amo1 <- data.frame(t(laply(amo1_UNIFORM, 
                                        function(x) x$arm_sequence)))
save(arm_seq_UNIFORM_amo1,
     file = "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/real_experiments/arm_seq_UNIFORM_amo1.Rda")
rm(amo1_UNIFORM)
gc()

library(ggplot)
library(ggjoy)
library(dplyr)
library(tidyr)
load(paste0(current_path, "arm_seq_UNIFORM_amo1.Rda"))
load(paste0(current_path, "arm_seq_APT_amo1.Rda"))
load(paste0(current_path, "arm_seq_AugUCB_amo1.Rda"))
load(paste0(current_path, "arm_seq_BUCB_amo1.Rda"))
load(paste0(current_path, "arm_seq_KLUCB_amo1.Rda"))

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

arm_seq_KLUCB_amo1 %>% tbl_df() %>% 
  mutate(index = 1:7580) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::filter(index > 10) %>%
  dplyr::group_by(arm, index) %>%
  dplyr::summarize(n = n()) %>%
  ungroup() %>%
  mutate(arm = as.factor(arm)) %>%
  ggplot(aes(x = index, y = arm, height = n)) +
  geom_joy(stat = "identity", scale = 1.2, size = 0.2) + theme_joy() +
  labs(title = "KL-UCB", y = "Arm", x = "Round",
       subtitle = "Pulls per arm and round after 5000 simulations")

arm_seq_AugUCB_amo1 %>% tbl_df() %>% 
  mutate(index = 1:7580) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::filter(index > 10) %>%
  dplyr::group_by(arm, index) %>%
  dplyr::summarize(n = n()) %>%
  ungroup() %>%
  mutate(arm = as.factor(arm)) %>%
  ggplot(aes(x = index, y = arm, height = n)) +
  geom_joy(stat = "identity", scale = 1.2, size = 0.2) + theme_joy() +
  labs(title = "Augmented-UCB", y = "Arm", x = "Round",
       subtitle = "Pulls per arm and round after 5000 simulations")

arm_seq_UNIFORM_amo1 %>% tbl_df() %>% 
  mutate(index = 1:7580) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::filter(index > 10) %>%
  dplyr::group_by(arm, index) %>%
  dplyr::summarize(n = n()) %>%
  ungroup() %>%
  mutate(arm = as.factor(arm)) %>%
  ggplot(aes(x = index, y = arm, height = n)) +
  geom_joy(stat = "identity", scale = 1.2, size = 0.2) + theme_joy() +
  labs(title = "UNIFORM", y = "Arm", x = "Round",
       subtitle = "Pulls per arm and round after 5000 simulations")

######################################################################

amo1_tidy_UNIFORM <- arm_seq_UNIFORM_amo1 %>% tbl_df() %>% 
  mutate(index = 1:7580) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::filter(index > 10) %>%
  dplyr::group_by(arm, index) %>%
  dplyr::summarize(n = n()) %>%
  ungroup()

amo1_tidy_AugUCB <- arm_seq_AugUCB_amo1 %>% tbl_df() %>% 
  mutate(index = 1:7580) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::filter(index > 10) %>%
  dplyr::group_by(arm, index) %>%
  dplyr::summarize(n = n()) %>%
  ungroup()

amo1_tidy_KLUCB <- arm_seq_KLUCB_amo1 %>% tbl_df() %>% 
  mutate(index = 1:7580) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::filter(index > 10) %>%
  dplyr::group_by(arm, index) %>%
  dplyr::summarize(n = n()) %>%
  ungroup()

amo1_tidy_BUCB <- arm_seq_BUCB_amo1 %>% tbl_df() %>% 
  mutate(index = 1:7580) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::filter(index > 10) %>%
  dplyr::group_by(arm, index) %>%
  dplyr::summarize(n = n()) %>%
  ungroup()

amo1_tidy_APT <- arm_seq_APT_amo1 %>% tbl_df() %>% 
  mutate(index = 1:7580) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::filter(index > 10) %>%
  dplyr::group_by(arm, index) %>%
  dplyr::summarize(n = n()) %>%
  ungroup()

amo1_tidy_all <- rbind(cbind(amo1_tidy_UNIFORM, model = rep("UNIFORM", times = 75700)),
                       cbind(amo1_tidy_APT, model = rep("APT", times = 75700)),
                       cbind(amo1_tidy_AugUCB, model = rep("AugUCB", times = 74605)),
                       cbind(amo1_tidy_BUCB, model = rep("BUCB", times = 75700)),
                       cbind(amo1_tidy_KLUCB, model = rep("KLUCB", times = 75700)))

save(amo1_tidy_APT, amo1_tidy_BUCB, amo1_tidy_UNIFORM,
     amo1_tidy_KLUCB, amo1_tidy_AugUCB, amo1_tidy_all,
     file = paste0(current_path, "amo1_tidy_data.Rda"))
load(file = paste0(current_path, "amo1_tidy_data.Rda"))

head(amo1_tidy_all)


amo1_tidy_all %>%
  filter(arm == 1) %>%
  ggplot(aes(x = index, y = n, group = model, color = model)) +
  geom_line()

amo1_tidy_all %>%
  filter(arm == 2) %>%
  ggplot(aes(x = index, y = n, group = model, color = model)) +
  geom_line()

amo1_tidy_all %>%
  filter(arm == 3) %>%
  ggplot(aes(x = index, y = n, group = model, color = model)) +
  geom_line()

amo1_tidy_all %>%
  filter(arm == 4) %>%
  ggplot(aes(x = index, y = n, group = model, color = model)) +
  geom_line()

amo1_tidy_all %>%
  filter(arm == 5) %>%
  ggplot(aes(x = index, y = n, group = model, color = model)) +
  geom_line()

amo1_tidy_all %>%
  filter(arm == 6) %>%
  ggplot(aes(x = index, y = n, group = model, color = model)) +
  geom_line()

amo1_tidy_all %>%
  filter(arm == 7) %>%
  ggplot(aes(x = index, y = n, group = model, color = model)) +
  geom_line()

amo1_tidy_all %>%
  filter(arm == 8) %>%
  ggplot(aes(x = index, y = n, group = model, color = model)) +
  geom_line()

amo1_tidy_all %>%
  dplyr::filter(arm == 9) %>%
  ggplot(aes(x = index, y = n, group = model, color = model)) +
  geom_line()

amo1_tidy_all %>%
  dplyr::filter(arm == 10) %>%
  ggplot(aes(x = index, y = n, group = model, color = model)) +
  geom_line()

amo1_boxplot_UNIFORM <- arm_seq_UNIFORM_amo1 %>% tbl_df() %>% 
  mutate(index = 1:7580) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::group_by(iter, arm) %>%
  dplyr::summarize(n = n()) %>%
  ungroup()

amo1_boxplot_AugUCB <- arm_seq_AugUCB_amo1 %>% dplyr::tbl_df() %>% 
  mutate(index = 1:7580) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::group_by(iter, arm) %>%
  dplyr::summarize(n = n()) %>%
  ungroup()

amo1_boxplot_KLUCB <- arm_seq_KLUCB_amo1 %>% tbl_df() %>% 
  mutate(index = 1:7580) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::group_by(iter, arm) %>%
  dplyr::summarize(n = n()) %>%
  ungroup()

amo1_boxplot_BUCB <- arm_seq_BUCB_amo1 %>% tbl_df() %>% 
  mutate(index = 1:7580) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::group_by(iter, arm) %>%
  dplyr::summarize(n = n()) %>%
  ungroup()

amo1_boxplot_APT <- arm_seq_APT_amo1 %>% tbl_df() %>% 
  mutate(index = 1:7580) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::group_by(iter, arm) %>%
  dplyr::summarize(n = n()) %>%
  ungroup()

amo1_boxplot_all <- rbind(cbind(amo1_boxplot_UNIFORM, model = rep("UNIFORM", times = 50000)),
                       cbind(amo1_boxplot_APT, model = rep("APT", times = 50000)),
                       cbind(amo1_boxplot_AugUCB, model = rep("AugUCB", times = 50000)),
                       cbind(amo1_boxplot_BUCB, model = rep("BUCB", times = 50000)),
                       cbind(amo1_boxplot_KLUCB, model = rep("KLUCB", times = 50000)))

amo1_boxplot_all %>%
  dplyr::filter(arm == 1) %>%
  ggplot(aes(x = model, y = n, group = model)) +
  geom_boxplot()

amo1_boxplot_all %>%
  dplyr::filter(arm == 2) %>%
  ggplot(aes(x = model, y = n, group = model)) +
  geom_boxplot()

amo1_boxplot_all %>%
  dplyr::filter(arm == 3) %>%
  ggplot(aes(x = model, y = n, group = model)) +
  geom_boxplot()

amo1_boxplot_all %>%
  dplyr::filter(arm == 4) %>%
  ggplot(aes(x = model, y = n, group = model)) +
  geom_boxplot()

amo1_boxplot_all %>%
  dplyr::filter(arm == 5) %>%
  ggplot(aes(x = model, y = n, group = model)) +
  geom_boxplot()

amo1_boxplot_all %>%
  dplyr::filter(arm == 6) %>%
  ggplot(aes(x = model, y = n, group = model)) +
  geom_boxplot()

amo1_boxplot_all %>%
  dplyr::filter(arm == 7) %>%
  ggplot(aes(x = model, y = n, group = model)) +
  geom_boxplot()

amo1_boxplot_all %>%
  dplyr::filter(arm == 8) %>%
  ggplot(aes(x = model, y = n, group = model)) +
  geom_boxplot()

amo1_boxplot_all %>%
  dplyr::filter(arm == 9) %>%
  ggplot(aes(x = model, y = n, group = model)) +
  geom_boxplot()

amo1_boxplot_all %>%
  dplyr::filter(arm == 10) %>%
  ggplot(aes(x = model, y = n, group = model)) +
  geom_boxplot()

######################################################################

library(zoo)
library(tidyr)
library(ggplot2)
load(paste0(current_path, "data_amo1_wide.Rda"))

names(data_amo1_wide) <- paste0(1:10)

amo1_hourly_window <- data.frame(index = 1:20101,
           lapply(data_amo1_wide, FUN = rollapply, 
                  width = 60, mean)) %>%
  gather(arm, mean, -index)

amo1_12hour_window <- data.frame(index = 1:19441,
                                 lapply(data_amo1_wide, FUN = rollapply, 
                                        width = 12*60, mean)) %>%
  gather(arm, mean, -index)

amo1_daily_window <- data.frame(index = 1:18721,
                                lapply(data_amo1_wide, FUN = rollapply, 
                                       width = 24*60, mean)) %>%
  gather(arm, mean, -index)

amo1_weekly_window <- data.frame(index = 1:10081,
                                lapply(data_amo1_wide, FUN = rollapply, 
                                       width = 24*60*7, mean)) %>%
  gather(arm, mean, -index)

amo1_7580_window <- data.frame(index = 1:12581,
                                lapply(data_amo1_wide, FUN = rollapply, 
                                       width = 7580, mean)) %>%
  dplyr::filter(index < 5000) %>%
  gather(arm, mean, -index)

amo1_hourly_window %>%
  ggplot(aes(x = index, y = mean, group = arm, color = arm)) +
  geom_line()

amo1_12hour_window %>%
  ggplot(aes(x = index, y = mean, group = arm, color = arm)) +
  geom_line()

amo1_daily_window %>%
  ggplot(aes(x = index, y = mean, group = arm, color = arm)) +
  geom_line()

amo1_weekly_window %>%
  ggplot(aes(x = index, y = mean, group = arm, color = arm)) +
  geom_line()

amo1_7580_window %>%
  ggplot(aes(x = index, y = mean, group = arm, color = arm)) +
  geom_line() 

######################################################################

library(zoo)
library(tidyr)
library(ggplot2)
load(paste0(current_path, "data_amo2_wide.Rda"))

names(data_amo2_wide) <- paste0(1:10)

amo2_hourly_window <- data.frame(index = 1:20101,
                                 lapply(data_amo2_wide, FUN = rollapply, 
                                        width = 60, mean)) %>%
  gather(arm, mean, -index)

amo2_12hour_window <- data.frame(index = 1:19441,
                                 lapply(data_amo2_wide, FUN = rollapply, 
                                        width = 12*60, mean)) %>%
  gather(arm, mean, -index)

amo2_daily_window <- data.frame(index = 1:18721,
                                lapply(data_amo2_wide, FUN = rollapply, 
                                       width = 24*60, mean)) %>%
  gather(arm, mean, -index)

amo2_weekly_window <- data.frame(index = 1:10081,
                                 lapply(data_amo2_wide, FUN = rollapply, 
                                        width = 24*60*7, mean)) %>%
  gather(arm, mean, -index)

amo2_7580_window <- data.frame(index = 1:12581,
                               lapply(data_amo2_wide, FUN = rollapply, 
                                      width = 7580, mean)) %>%
  dplyr::filter(index < 5000) %>%
  gather(arm, mean, -index)

amo2_hourly_window %>%
  ggplot(aes(x = index, y = mean, group = arm, color = arm)) +
  geom_line()

amo2_12hour_window %>%
  ggplot(aes(x = index, y = mean, group = arm, color = arm)) +
  geom_line()

amo2_daily_window %>%
  ggplot(aes(x = index, y = mean, group = arm, color = arm)) +
  geom_line()

amo2_weekly_window %>%
  ggplot(aes(x = index, y = mean, group = arm, color = arm)) +
  geom_line()

amo2_7580_window %>%
  ggplot(aes(x = index, y = mean, group = arm, color = arm)) +
  geom_line() 
