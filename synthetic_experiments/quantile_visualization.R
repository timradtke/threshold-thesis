# After four rounds, only 0 observations, alpha = 0.04 and beta = 0.96

alpha_prime <- 0.04
beta_prime <- 0.96+4
(mean_post <- alpha_prime/(alpha_prime+beta_prime))
mean <- 0
# tau = 0.05
Ti <- 4
N <- 2500
log(N)

plot(seq(0,1,by = 0.001),dbeta(seq(0,1,by = 0.001), alpha_prime, beta_prime), type = "l")
plot(seq(0,1,by = 0.001),pbeta(seq(0,1,by = 0.001), alpha_prime, beta_prime), type = "l")
round(qbeta(1 - 1/Ti, alpha_prime, beta_prime),5)
round(qbeta(1 - 1/sqrt(Ti), alpha_prime, beta_prime),5)
round(qbeta(1 - 1/(Ti)^2, alpha_prime, beta_prime),5)
round(qbeta(1 - 1/(Ti)^5, alpha_prime, beta_prime),5)
round(qbeta(1 - 1/((Ti)*log(N)), alpha_prime, beta_prime),5)
round(qbeta(1 - 1/((Ti)*log(N)^5), alpha_prime, beta_prime),5)
round(qbeta(1 - 1/((Ti)*N), alpha_prime, beta_prime),5)
