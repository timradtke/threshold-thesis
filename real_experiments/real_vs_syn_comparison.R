current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/AugUCB/"
source(paste0(current_path, "run_from_data.R"))
source(paste0(current_path, "run_parallelized_from_data.R"))
current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/real_experiments/"

load(paste0(current_path, "data_amo1_mean_firsthalf"))
load(paste0(current_path, "data_amo1_mean_secondhalf"))
load(paste0(current_path, "data_amo1.Rda"))
tau_amo1 <- 3/60
epsilon_amo1 <- 1/60

data <- data_amo1[[1]]
rm(data_amo1)
gc()

library(dplyr)
library(tidyr)
library(ggplot2)

data.frame(index = 1:7521,
           lapply(data, FUN = zoo::rollapply, width = 60, mean)) %>%
  gather(arm, mean, -index) %>% 
  ggplot(aes(x = index, y = mean, group = arm, color = arm)) +
  geom_line() + geom_hline(aes(yintercept = tau_amo1))

data.frame(index = 1:6141,
           lapply(data, FUN = zoo::rollapply, width = 24*60, mean)) %>%
  gather(arm, mean, -index) %>% 
  ggplot(aes(x = index, y = mean, group = arm, color = arm)) +
  geom_line() + geom_hline(aes(yintercept = tau_amo1))

data.frame(index = 1:4701,
           lapply(data, FUN = zoo::rollapply, width = 2*24*60, mean)) %>%
  gather(arm, mean, -index) %>% 
  ggplot(aes(x = index, y = mean, group = arm, color = arm)) +
  geom_line() + geom_hline(aes(yintercept = tau_amo1))


real_mean <- unlist(lapply(data, mean))
H_real <- get_complexity(real_mean, tau_amo1, epsilon_amo1)


syn_data <- data.frame(rep(NA, times = 7580))
set.seed(124632728)
for(i in 1:length(real_mean)) {
  syn_data[[i]] <- as.numeric(purrr::rbernoulli(7580, p  = real_mean[i]))
}
names(syn_data) <- names(real_mean)

data.frame(index = 1:7521,
           lapply(syn_data, FUN = zoo::rollapply, width = 60, mean)) %>%
  gather(arm, mean, -index) %>% 
  ggplot(aes(x = index, y = mean, group = arm, color = arm)) +
  geom_line() + geom_hline(aes(yintercept = tau_amo1))

data.frame(index = 1:6141,
           lapply(syn_data, FUN = zoo::rollapply, width = 24*60, mean)) %>%
  gather(arm, mean, -index) %>% 
  ggplot(aes(x = index, y = mean, group = arm, color = arm)) +
  geom_line() + geom_hline(aes(yintercept = tau_amo1))

data.frame(index = 1:4701,
           lapply(syn_data, FUN = zoo::rollapply, width = 2*24*60, mean)) %>%
  gather(arm, mean, -index) %>% 
  ggplot(aes(x = index, y = mean, group = arm, color = arm)) +
  geom_line() + geom_hline(aes(yintercept = tau_amo1))



  ###############################################
real_daily <- data.frame(source = rep("Real", 6141),
                         index = 1:6141,
                        lapply(data, FUN = zoo::rollapply, width = 24*60, mean)) %>%
  gather(arm, mean, -index, -source)
syn_daily <- data.frame(source = rep("Synthetic", 6141),
                        index = 1:6141,
                        lapply(syn_data, FUN = zoo::rollapply, width = 24*60, mean)) %>%
  gather(arm, mean, -index, -source)

rbind(real_daily, syn_daily) %>%
  ggplot(aes(x = index, y = mean, group = arm, color = arm)) +
  geom_line() + geom_hline(aes(yintercept = tau_amo1)) +
  facet_grid(.~source)

real_2daily <- data.frame(source = rep("Real", 4701),
                         index = 1:4701,
                         lapply(data, FUN = zoo::rollapply, width = 2*24*60, mean)) %>%
  gather(arm, mean, -index, -source)
syn_2daily <- data.frame(source = rep("Synthetic", 4701),
                        index = 1:4701,
                        lapply(syn_data, FUN = zoo::rollapply, width = 2*24*60, mean)) %>%
  gather(arm, mean, -index, -source)

rbind(real_2daily, syn_2daily) %>%
  ggplot(aes(x = index, y = mean, group = arm, color = arm)) +
  geom_line() + geom_hline(aes(yintercept = tau_amo1)) +
  facet_grid(.~source)
  


  