diarmid <- function(x, tau) {
  n <- length(x)
  mu <- ifelse(mean(x) == 1,
               1-1/(n)^2, mean(x))
  mu <- ifelse(mu == 0,
               1/(n)^2, mu)
  abs(
    mu*log(mu/(mu-1/n)) +
      (1-mu)*log((1-mu)/(1-mu+1/n)) +
      1/n*log((1-tau)/tau) +
      1/n*log((mu-1/n)/(1-mu+1/n))
  )
}

y <- c(1,purrr::rbernoulli(10000,0.05))

plot(c(1,10000), c(-15,0), type = "n")
for(i in 1:10000) {
  points(i, log(diarmid(y[1:i], 0.1)), pch = 19, cex = 0.1)
}

y <- purrr::rbernoulli(10000,0.01)

plot(c(1,10000), c(1,20), type = "n")
for(i in 1:10000) {
  points(i, diarmid(y[1:i], 0.1), pch = 19, cex = 0.5)
}

y <- purrr::rbernoulli(10000,0.001)

plot(c(1,10000), c(1,20), type = "n")
for(i in 1:10000) {
  points(i, diarmid(y[1:i], 0.1), pch = 19, cex = 0.5)
}

y <- purrr::rbernoulli(10000,0.01)

plot(c(1,10000), c(1,20), type = "n")
for(i in 1:10000) {
  points(i, diarmid(y[1:i], 0.4), pch = 19, cex = 0.5)
}



plot(c(1,100000), c(0,0.2), type = "n")
for(j in seq(0.001,0.999, by = 0.05)) {
  y <- purrr::rbernoulli(100000,j)
  for(i in seq(10,100000, by=50)) {
    points(i, diarmid(y[1:i], 0.5), pch = 19, cex = 0.1)
  }
}


plot(3:100, log(-log(1-0.5/(3:100))), type = "l")
lines(3:100, log(0.5/(3:100)), col = "red")
