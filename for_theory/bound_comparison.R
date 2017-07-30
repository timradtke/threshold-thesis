current_path <- "/Users/timradtke/Dropbox/1Master/Master Thesis/threshold-thesis/AugUCB/"
source(paste0(current_path, "run_from_data.R"))
source(paste0(current_path, "run_parallelized_from_data.R"))


get_chernoff_info <- function(x,y) {
  if(x < y) {
    S <- seq(x,y,by=0.00001)
  } else {
    S <- seq(y,x,by=0.00001)
  }
  klx <- kl_ber(S,x)
  kly <- kl_ber(S,y)
  kl_ber(S[which.min(abs(klx-kly))],y)
}

get_chernoff_info(0.2,0.4)

test_grid <- expand.grid(p = seq(0.001,0.999,0.001),
                         tau = seq(0.001,0.999,0.001))
test_grid <- test_grid[test_grid$p != test_grid$tau,]
test_grid$chernoff <- NA
test_grid$DeltaSquared <- (test_grid$p - test_grid$tau)^2
dim(test_grid)
#for(i in 1:10) {
#  test_grid[i,3] <- get_chernoff_info(test_grid[i,1], test_grid[i,2])
#}

library(dplyr)
test_grid %>% filter(p > 0.48, p<0.52)
test_grid %>% filter(chernoff < DeltaSquared)
dim(test_grid[is.na(test_grid$chernoff),])


#############################################

get_complexity
get_chernoff_complexity <- function(means, tau, epsilon) {
  res <- vector(length = length(means))
  for(i in 1:length(means)) {
    res[i] <- max(get_chernoff_info(means[i], tau), epsilon^2/2)
  }
  sum(1/res)
}
get_APT_complexity <- function(means, tau, epsilon) {
  res <- vector(length = length(means))
  for(i in 1:length(means)) {
    res[i] <- (abs(means[i] - tau) + epsilon)^2
  }
  sum(1/res)
}

mean_loc7nt <- c(10^-4, 10^-3.5, 10^-3.25, 10^-3, 
                 10^-2.75, 10^-2.5, 10^-2.25, 10^-1.85, 10^-1.75,
                 10^-1.5, 10^-1)
tau_loc7nt <- 10^-2
epsilon_loc7nt <- 0

mean_loc7nt <- c(10^-4, 10^-3.5, 10^-3.25, 10^-3, 
                 10^-2.75, 10^-2.5, 10^-2.25, 10^-1.85, 10^-1.75,
                 10^-1.5, 10^-1)
tau_loc7nt <- 10^-2
epsilon_loc7nt <- 0

get_complexity(mean_loc7nt, tau_loc7nt, epsilon_loc7nt)
get_APT_complexity(mean_loc7nt, tau_loc7nt, epsilon_loc7nt)/64
get_chernoff_complexity(mean_loc7nt, tau_loc7nt, epsilon_loc7nt)

APT_bound <- function(budget, complexity, K, R = 1) {
  exp(-1/64/(R^2)*budget/complexity + 2*log((log(budget)+1)*K))
}

chernoff_bound <- function(budget, complexity, K, R = 1) {
  exp(-1/2*budget/complexity)*K*complexity
}



APT_bound(10000, get_APT_complexity(mean_loc7nt, tau_loc7nt, epsilon_loc7nt), 11)
chernoff_bound(10000, get_chernoff_complexity(mean_loc7nt, tau_loc7nt, epsilon_loc7nt), 11)

APT_bound(20000, get_APT_complexity(mean_loc7nt, tau_loc7nt, epsilon_loc7nt), 11)
chernoff_bound(20000, get_chernoff_complexity(mean_loc7nt, tau_loc7nt, epsilon_loc7nt), 11)

APT_bound(50000, get_APT_complexity(mean_loc7nt, tau_loc7nt, epsilon_loc7nt), 11)
chernoff_bound(50000, get_chernoff_complexity(mean_loc7nt, tau_loc7nt, epsilon_loc7nt), 11)

APT_bound(100000, get_APT_complexity(mean_loc7nt, tau_loc7nt, epsilon_loc7nt), 11)
chernoff_bound(100000, get_chernoff_complexity(mean_loc7nt, tau_loc7nt, epsilon_loc7nt), 11)

APT_bound(200000, get_APT_complexity(mean_loc7nt, tau_loc7nt, epsilon_loc7nt), 11)
chernoff_bound(200000, get_chernoff_complexity(mean_loc7nt, tau_loc7nt, epsilon_loc7nt), 11)



mean_loc7nt <- c(0.4-0.2^(1:4),
               0.45, 0.55,
               0.6+0.2^(5-(1:4)))
tau_loc7nt <- 0.5
epsilon_loc7nt <- 0.1

k <- 10
APT_bound(10000, get_APT_complexity(mean_loc7nt, tau_loc7nt, epsilon_loc7nt), k)
chernoff_bound(10000, get_chernoff_complexity(mean_loc7nt, tau_loc7nt, epsilon_loc7nt), k)

APT_bound(20000, get_APT_complexity(mean_loc7nt, tau_loc7nt, epsilon_loc7nt), k)
chernoff_bound(20000, get_chernoff_complexity(mean_loc7nt, tau_loc7nt, epsilon_loc7nt), k)

APT_bound(50000, get_APT_complexity(mean_loc7nt, tau_loc7nt, epsilon_loc7nt), k)
chernoff_bound(50000, get_chernoff_complexity(mean_loc7nt, tau_loc7nt, epsilon_loc7nt), k)

APT_bound(100000, get_APT_complexity(mean_loc7nt, tau_loc7nt, epsilon_loc7nt), k)
chernoff_bound(100000, get_chernoff_complexity(mean_loc7nt, tau_loc7nt, epsilon_loc7nt), k)

APT_bound(200000, get_APT_complexity(mean_loc7nt, tau_loc7nt, epsilon_loc7nt), k)
chernoff_bound(200000, get_chernoff_complexity(mean_loc7nt, tau_loc7nt, epsilon_loc7nt), k)
