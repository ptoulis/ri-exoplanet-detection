#' Code from Example 1, Section 1 in Toulis, P. and Bean, J. (2021). "Randomization Inference of Periodicity in Unequally Spaced Time Series with Application to Exoplanet Detection."
#' https://www.ptoulis.com/s/astro_main.pdf
#' 
#' Panos Toulis, ptoulis@chicagobooth.edu, 03/2021
#' 
rm(list=ls())
source("rv_lib.R")

set.seed(1) # Seed the dataset
n = 100 # datapoints
var_days = 0.05 # daily variation
days = seq(1, n) + var_days*runif(n, min=-1, max=1) # observation times
theta0 = sqrt(2) # true periodicity

y = 1.5 * cos(2*pi*days/theta0) + rnorm(n)  # synthetic RV data.
rv = data.frame(val=y, days=days, se=rep(1, n))  # RV data frame.
all_P = 10^seq(-1, log(50, base=10), length.out=10000)
#
out = lombe_scragle_fast(rv, all_P, v=T) # LS periodogram
out$Phat # Periodogram peak.

# Generate sampling distribution.
sampl = replicate(100, {
  days = seq(1, n) + var_days*runif(n, min=-1, max=1) # sample observation times T^n
  y = 1.5 * cos(2*pi*days/theta0) + rnorm(n) # Sample data $Y^n
  rv = data.frame(val=y, days=days, se=rep(1, n)) # Create periodogram
  out = lombe_scragle_fast(rv, all_P, v=F) # Calculate LS periodogram.
  out$Phat # Periodogram peak.
})


par(mfrow=c(1, 2))
par(mar=rep(4.5, 4))

#plot(days, y, pch=20, xlab="time", ylab="radial velocity", cex.lab=1.5)
# lines(days, y, col=rgb(0, 0, 1, alpha=0.2))
out = lombe_scragle_fast(rv, all_P, v=F)
plot(all_P, out$Pf, log="x", type="l", xlab="period", ylab="power", 
     cex.lab=1.5)
periods = c(theta0, out$Phat)
Pff = sapply(periods, function(p) out$Pf[which.min(abs(p-all_P))])
points(c(theta0, out$Phat), Pff, pch=20, cex=2, col=c("red", "blue"))
legend("bottomright", legend = c("true period", "peak"), col = c("red", "blue"),pch = c(20, 20))

# Generates
hist(sampl, breaks=50, 
     main="Sampling distribution of periodogram peak", 
     cex.lab=1.5, xlab="value", freq=F, 
     cex.main=1.5)

