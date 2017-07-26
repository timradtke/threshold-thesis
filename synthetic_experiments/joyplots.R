install.packages("devtools")
library(devtools)
install_github("clauswilke/ggjoy")
library(ggplot2)
library(ggjoy)

ggplot(diamonds, aes(x = price, y = cut)) +
  geom_joy(scale = 4) + theme_joy() +
  scale_y_discrete(expand = c(0.01, 0)) +   # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0))      # for both axes to remove unneeded padding
#> Picking joint bandwidth of 458

count_data <- data.frame(group = rep(letters[1:5], each = 20),
                         mean = rep(1:5, each = 20))
count_data$count <- rpois(nrow(count_data), count_data$mean)
ggplot(count_data, aes(x = count, y = group, fill = group)) +
  geom_joy(stat = "binline", binwidth = 1, scale = 0.9) +
  theme_joy()

count_data2 <- data.frame(group = rep(letters[1:5], each = 10),
                          index = rep(1:10, 5),
                          mean = rep(c(1,2,3,4,5,4,3,2,1,1),5))
count_data2$count <- rpois(nrow(count_data2), count_data2$mean)
ggplot(count_data2, aes(x = index, y = group, height = count, fill = group)) +
  geom_joy(stat = "identity", scale = 0.7)

current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/synthetic_experiments/"
load(file = paste0(current_path, "loc2_BUCB_horizon_4000.Rda"))
load(file = paste0(current_path, "loc2_APT_2000.Rda"))
load(file = paste0(current_path, "loc2_UNIFORM_2000.Rda"))
load(file = paste0(current_path, "loc2_AugUCB_2000.Rda"))
load(file = paste0(current_path, "loc2_KLUCB80_long.Rda"))



library(plyr)

arm_seq_APT <- data.frame(t(laply(loc2_APT_2000, function(x) x$arm_sequence)))
arm_seq_BUCB <- data.frame(t(laply(loc2_BUCB_horizon_4000, function(x) x$arm_sequence)))
arm_seq_AugUCB <- data.frame(t(laply(loc2_AugUCB_2000, function(x) x$arm_sequence)))
arm_seq_KLUCB <- data.frame(t(laply(loc2_KLUCB80_long, function(x) x$arm_sequence)))

save(arm_seq_APT, arm_seq_BUCB, arm_seq_AugUCB, arm_seq_KLUCB,
     file = paste0(current_path, "arm_seq_loc2.Rda"))
load(paste0(current_path, "arm_seq_loc2.Rda"))
library(tidyr)
library(dplyr)

arm_seq_APT %>% tbl_df() %>% 
  mutate(index = 1:2000) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::filter(index > 10) %>%
  dplyr::group_by(arm, index) %>%
  dplyr::summarize(n = n()) %>%
  ungroup() %>%
  mutate(arm = as.factor(arm)) %>%
  ggplot(aes(x = index, y = arm, height = n)) +
  geom_joy(stat = "identity", scale = 1.5, size = 0.2) + theme_joy() +
  labs(title = "APT")

arm_seq_BUCB %>% tbl_df() %>%
  mutate(index = 1:4000) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::filter(index > 10, index <= 2000) %>%
  dplyr::group_by(arm, index) %>%
  dplyr::summarize(n = n()) %>%
  ungroup() %>%
  mutate(arm = as.factor(arm)) %>%
  ggplot(aes(x = index, y = arm, height = n)) +
  geom_joy(stat = "identity", scale = 1.5, size = 0.2) + theme_joy() +
  labs(title = "Bayes-UCB")

arm_seq_KLUCB %>% tbl_df() %>%
  mutate(index = 1:1500) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::filter(index > 10, index <= 2000) %>%
  dplyr::group_by(arm, index) %>%
  dplyr::summarize(n = n()) %>%
  ungroup() %>%
  mutate(arm = as.factor(arm)) %>%
  ggplot(aes(x = index, y = arm, height = n)) +
  geom_joy(stat = "identity", scale = 1.5, size = 0.2) + theme_joy() +
  labs(title = "KL-UCB")

arm_seq_AugUCB %>% tbl_df() %>%
  mutate(index = 1:2000) %>%
  gather(key = iter, value = arm, -index) %>%
  dplyr::filter(index > 10, index <= 2000) %>%
  dplyr::group_by(arm, index) %>%
  dplyr::summarize(n = n()) %>%
  ungroup() %>%
  mutate(arm = as.factor(arm)) %>%
  ggplot(aes(x = index, y = arm, height = n)) +
  geom_joy(stat = "identity", scale = 1.5, size = 0.2) + theme_joy()+
  labs(title = "Augmented-UCB")
  
arm_seq_APT %>% tbl_df() %>% 
  mutate(index = 1:2000) %>%
  gather(key = iter, value = arm, -index) %>%
  mutate(arm = as.factor(arm)) %>%
  ggplot(aes(x = index, y = arm)) +
  geom_joy(scale = 2) + theme_joy() +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "APT")

arm_seq_BUCB %>% tbl_df() %>% 
  mutate(index = 1:4000) %>%
  gather(key = iter, value = arm, -index) %>%
  filter(index <= 2000) %>%
  mutate(arm = as.factor(arm)) %>%
  ggplot(aes(x = index, y = arm)) +
  geom_joy(scale = 2) + theme_joy() +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Bayes-UCB")

arm_seq_AugUCB %>% tbl_df() %>% 
  mutate(index = 1:2000) %>%
  gather(key = iter, value = arm, -index) %>%
  filter(index <= 2000) %>%
  mutate(arm = as.factor(arm)) %>%
  ggplot(aes(x = index, y = arm)) +
  geom_joy(scale = 2) + theme_joy() +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Augmented-UCB")


loc2_fullcomp_BUCB_horizon_4000 <- compare_to_ground_truth(mean_loc2, 
                                                      loc2_BUCB_horizon_4000,
                                                      tau_loc2,
                                                      epsilon_loc2)
head(loc2_fullcomp_BUCB_horizon_4000)


get_wrong_arms_per_round <- function(true_means, sim_res, tau, epsilon) {
  n <- length(true_means)
  message(paste0("Comparing ", length(sim_res), " simulations."))
  true_classification_up <- which(true_means > tau+epsilon)
  true_classification_down <- which(true_means < tau-epsilon)
  
  get_wrong_arms_in_iter <- function(x, true_class_up, true_class_down) {
    arms <- sort(c(true_class_up[!(true_class_up %in% which(x >= tau))],
           true_class_down[!(true_class_down %in% which(x < tau))]))
    ifelse(length(arms) == 0, NA, arms)
  }
  
  get_wrong_arms <- function(res_df, ...) {
    larms <- alply(res_df$mean_storage, .margins = 1, get_wrong_arms_in_iter, ...,
                   .expand = TRUE)
    dfarms <- data.frame(round = rep(0, n), arm = 1:n)
    for(i in 1:length(larms)) {
      dfarms <- rbind(dfarms,
                      data.frame(round = rep(i, length(larms[[i]])),
                                 arm = larms[[i]]))
    }
    return(dfarms)
  }
  
  comp_list <- ldply(sim_res, get_wrong_arms, 
                     true_class_up = true_classification_up,
                     true_class_down = true_classification_down,
                     .progress = "text")
  
  return(comp_list)
}

loc2_BUCB_wrong_arms_per_round <- get_wrong_arms_per_round(mean_loc2, 
                                                      loc2_BUCB_horizon_4000,
                                                      tau_loc2,
                                                      epsilon_loc2)
