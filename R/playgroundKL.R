x <- seq(-10, 10, by = 0.01)
y <- (1+x^2)^0.5 # + x^(2/3)
y

plot(x,y,type = "l")
(-0.5)^(2/3)

p <- 0.25
q <- seq(0.01,0.99,by = 0.01)


dKL <- log(p/q)-log((1-p)/(1-q))
KL <- p*log(p/q) + (1-p)*log((1-p)/(1-q))
plot(q,KL, type = "l")
plot(q,dKL, type = "l")

bernoulli_KL <- function(p,q) {
  p*log(p/q) + (1-p)*log((1-p)/(1-q))
}

q <- 0.2
p <- seq(0.01,0.99,by = 0.01)
d <- abs(p-q)
B <- d*bernoulli_KL(p,q)
plot(p,B, type = "l")

bernoulli_KL(0.3,0.4)
bernoulli_KL(0.1,0.9)
bernoulli_KL(0.01,0.99)

beta_KL <- function(a1, b1, a2, b2) {
  log(beta(a2,b2)/beta(a1,b1)) + (a1-a2)*digamma(a1) +
    (b1-b2)*digamma(b1) + (a2-a1+b2-b1)*digamma(a1+b1)
}

beta_KL(1,1,3,3)
beta_KL(3,3,1,1)
beta_KL(5,5,3,3)
beta_KL(3,3,9,9)
beta_KL(9,9,3,3)
beta_KL(7,26,100,100)
beta_KL(25,26,100,100)
beta_KL(56,56,100,100)

beta_KL(100,100,7,26)
beta_KL(100,100,25,26)
beta_KL(100,100,56,56)

plot(seq(0,1,by=0.001), dbeta(seq(0,1,by=0.001), 100,100), type = "l")
lines(seq(0,1,by=0.001), dbeta(seq(0,1,by=0.001), 3,10), type = "l")
lines(seq(0,1,by=0.001), dbeta(seq(0,1,by=0.001), 56,56), type = "l")
lines(seq(0,1,by=0.001), dbeta(seq(0,1,by=0.001), 20,26), type = "l")

set.seed(512)
A <- rbernoulli(200, p = 0.2)
B <- rbernoulli(200, p = 0.5)
C <- rbernoulli(200, p = 0.3)
tau <- rep(1,200)

Acum <- cumsum(A)
Bcum <- cumsum(B)
Ccum <- cumsum(C)
taucum <- cumsum(tau)
A_tau_index <- seq(1,199,by = 2)
B_tau_index <- seq(2,200,by = 2)


beta_KL(Acum+1,(1:200)+10,3,10)
beta_KL(Bcum+1,(1:200)+10,3,10)
beta_KL(Ccum+1,(1:200)+10,3,10)
plot(1:200, beta_KL(Acum+1,(1:200)+10,3,10), type = "l", col = "blue")
lines(1:200, beta_KL(Bcum+1,(1:200)+10,3,10), col = "red")
lines(1:200, beta_KL(Ccum+1,(1:200)+10,3,10), col = "black", lty = 2)

plot(seq(0,1,by=0.001), dbeta(seq(0,1,by=0.001), 3,10), type = "l")
lines(seq(0,1,by=0.001), dbeta(seq(0,1,by=0.001), 3+Acum[200],10+200), col = "blue")
lines(seq(0,1,by=0.001), dbeta(seq(0,1,by=0.001), 3+Bcum[200],210), col = "red")
plot(seq(0,1,by=0.001), dbeta(seq(0,1,by=0.001), 3+Ccum[200],210),
      lty = 2, type = "l")

######################################################################

interval_prob <- function(tau, epsilon, a, b) {
  pbeta(tau+epsilon, a, b) - pbeta(tau-epsilon, a, b)
}

set.seed(1024)
A <- rbernoulli(1000, p = 0.2)
B <- rbernoulli(1000, p = 0.5)
C <- rbernoulli(1000, p = 0.3)
D <- rbernoulli(1000, p = 0.325)
Acum <- cumsum(A)
Bcum <- cumsum(B)
Ccum <- cumsum(C)
Dcum <- cumsum(D)
A_probs <- interval_prob(0.3, 0.0125, Acum+3, (1:1000)+10-Acum)
B_probs <- interval_prob(0.3, 0.0125, Bcum+3, (1:1000)+10-Bcum)
C_probs <- interval_prob(0.3, 0.0125, Ccum+3, (1:1000)+10-Ccum)
D_probs <- interval_prob(0.3, 0.0125, Dcum+3, (1:1000)+10-Dcum)
plot(c(0,1000),c(0,1), type = "n")
lines(1:1000, A_probs/(A_probs+B_probs+C_probs+D_probs), col = "blue")
lines(1:1000, B_probs/(A_probs+B_probs+C_probs+D_probs), col = "red")
lines(1:1000, C_probs/(A_probs+B_probs+C_probs+D_probs), lty = 2)
lines(1:1000, D_probs/(A_probs+B_probs+C_probs+D_probs), lty = 1)

plot(c(0,1),c(0,1), type = "n")
lines(seq(0,1,by=0.001), dbeta(seq(0,1,by=0.001), 3+Acum[1000],1010-Acum), col = "blue")
lines(seq(0,1,by=0.001), dbeta(seq(0,1,by=0.001), 3+Bcum[1000],1010-Bcum), col = "red")
lines(seq(0,1,by=0.001), dbeta(seq(0,1,by=0.001), 3+Ccum[1000],1010-Ccum), lty = 2)
lines(seq(0,1,by=0.001), dbeta(seq(0,1,by=0.001), 3+Dcum[1000],1010-Dcum), lty = 1)
abline(v=0.3-0.0125)
abline(v=0.3+0.0125)
