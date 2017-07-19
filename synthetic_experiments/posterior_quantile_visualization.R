plot(qbeta(1-1/(rep(6000,times=100)),0.04,0.96+1:100),type = "l")
abline(h = 0.04, lty = 2)
lines(qbeta(1-1/(1:100), 0.04, 0.96+1:100), lty = 2)
lines(qbeta(1-1/(6000*1:100),0.04,0.96+1:100),col="red")
lines(qbeta(1-1/(log(6000)*1:100),0.04,0.96+1:100),col="blue")
lines(qbeta(1-1/(1:100)^2,0.04,0.96+1:100),col="green")
#lines(qbeta(1-1/(log(2000)*(1:100)^2),0.04,0.96+1:100),col="green")
lines(qbeta(1-1/(log(2000^2)*1:100),0.04,0.96+1:100),col="blue")
lines(qbeta(1-1/(log(2000^5)*1:100),0.04,0.96+1:100),col="blue")

plot(qbeta(1/(rep(6000,times=100)),0.04+seq(0.5,50,by = 0.5),0.96+1:100),type = "l")
abline(h = 0.04, lty = 2)
lines(qbeta(1/(log(6000)*1:100),0.04+seq(0.5,50,by = 0.5),0.96+1:100),col="blue")
lines(qbeta(1-1/(1:100), 0.04, 0.96+1:100), lty = 2)
lines(qbeta(1-1/(6000*1:100),0.04,0.96+1:100),col="red")
lines(qbeta(1-1/(1:100)^2,0.04,0.96+1:100),col="green")
#lines(qbeta(1-1/(log(2000)*(1:100)^2),0.04,0.96+1:100),col="green")
lines(qbeta(1-1/(log(2000^2)*1:100),0.04,0.96+1:100),col="blue")
lines(qbeta(1-1/(log(2000^5)*1:100),0.04,0.96+1:100),col="blue")

plot(0.04-qbeta(1/(rep(6000,times=1000)),0.04+seq(0.06,60,by = 0.06),0.96+1:1000),type = "l",
     ylim = c(-0.04,0.08))
lines(0.04-qbeta(1/(log(6000)*1:1000),0.04+seq(0.06,60,by = 0.06),0.96+1:1000),col="blue")
lines(0.04-qbeta(1/(1:100)^2,0.04+seq(0.06,6,by = 0.06),0.96+1:100),col="red")
lines(qbeta(1-1/(rep(6000,times=1000)), 0.04, 0.96+1:1000)-0.04, lty = 2)
lines(qbeta(1-1/(log(6000)*1:100), 0.04, 0.96+1:100)-0.04, lty = 2, col = "blue")
lines(qbeta(1-1/(1:100)^2, 0.04, 0.96+1:100)-0.04, lty = 2, col = "red")



plot(qbeta(1-1/(rep(2000,times=100)),1,1+1:100),type = "l")
lines(qbeta(1-1/(2000*1:100),1,1+1:100),col="red")
lines(qbeta(1-1/(log(2000)*1:100),1,1+1:100),col="blue")
lines(qbeta(1-1/(1:100)^2,1,1+1:100),col="green")
lines(qbeta(1-1/(log(2000^2)*1:100),1,1+1:100),col="blue")
lines(qbeta(1-1/(log(2000^5)*1:100),1,1+1:100),col="blue")

1/0.04

plot(qbeta(1-1/(rep(2000,times=100)),1,24+1:100),type = "l")
lines(qbeta(1-1/(2000*1:100),1,24+1:100),col="red")
lines(qbeta(1-1/(log(2000)*1:100),1,24+1:100),col="blue")
lines(qbeta(1-1/(1:100)^2,1,24+1:100),col="green")
lines(qbeta(1-1/(log(2000^2)*1:100),1,24+1:100),col="blue")
lines(qbeta(1-1/(log(2000^5)*1:100),1,24+1:100),col="blue")



x1 <- purrr::rbernoulli(5000, 0.1)
x1_pulled <- x1[seq(10,5000,by=10)]
x1_pulled
x1_mean <- c()
for(i in 1:length(x1_pulled)) {
  if(mean(x1_pulled[1:i]) > 0.5) {
    x1_mean[i] <- mean(x1_pulled[1:i])-5000/i-
      sqrt(5000*mean(x1_pulled[1:i])*(1-mean(x1_pulled[1:i]))/i)
  } else {
    x1_mean[i] <- mean(x1_pulled[1:i])+5000/i+
      sqrt(5000*mean(x1_pulled[1:i])*(1-mean(x1_pulled[1:i]))/i)
  }
}
