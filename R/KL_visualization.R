# Visualization of KL Divergence for different Bernoullis

KL_ber <- function(p,q) p*log(p/q) + (1-p)*log((1-p)/(1-q))

KL_ber_approx <- function(p,q) (p-q)^2

plot(seq(0.01,0.99,by = 0.01), KL_ber(seq(0.01,0.99,by = 0.01),0.001), type = "l")
lines(seq(0.01,0.99,by = 0.01), KL_ber_approx(seq(0.01,0.99,by = 0.01), 0.001))

lines(seq(0.01,0.99,by = 0.01), KL_ber(seq(0.01,0.99,by = 0.01), 0.0075),
      col = "blue")
lines(seq(0.01,0.99,by = 0.01), KL_ber(seq(0.01,0.99,by = 0.01), 0.01),
      col = "green")
lines(seq(0.01,0.99,by = 0.01), KL_ber_approx(seq(0.01,0.99,by = 0.01), 0.01),
      col = "red")
lines(seq(0.01,0.99,by = 0.01), KL_ber_approx(seq(0.01,0.99,by = 0.01), 0.01),
      col = "red")

plot(seq(0.01,0.99,by = 0.01), KL_ber(seq(0.01,0.99,by = 0.01),0.001), type = "l")
for (i in c(0.001,0.01,0.1,0.5)) {
  lines(seq(0.01,0.99,by = 0.01), KL_ber(seq(0.01,0.99,by = 0.01), i),)
  lines(seq(0.01,0.99,by = 0.01), KL_ber_approx(seq(0.01,0.99,by = 0.01), i),
        col = "red")
}

lines(seq(0.01,0.99,by = 0.01), 8.5*seq(0.01,0.99,by = 0.01),
      lty=2)

#### Show asymmetry of KL when moving close to 0 or 1
#### (the divergence is different on the two sides of the threshold even though
#### the absolute difference between theta and tau is the same on both sides)
KL_comp <- function(theta1, theta2, tau) round(KL_ber(theta1,tau) - KL_ber(theta2,tau),5)

# compare theta=0.45 and 0.55 to tau=0.5 (should be similar)
KL_comp(0.45,0.55,0.5)  # 0

# compare theta=0.35 and 0.45 to tau=0.4
KL_comp(0.35,0.45,0.4)  # 0.00015

KL_comp(0.25,0.35, 0.3) # 0.00038

KL_comp(0.15,0.25, 0.2) # 0.001

KL_comp(0.05,0.15, 0.1) # 0.00447

KL_comp(0.01,0.11, 0.06) # 0.01536
KL_comp(0.01,0.05, 0.03) # 0.00347
KL_ber(0.01,0.03) # 0.00921866
KL_ber(0.05,0.03) # 0.005748899

for(i in 1:length(seq(0.001,0.)))
