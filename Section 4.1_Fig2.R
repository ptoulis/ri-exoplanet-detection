#' Data from Raup and Sepkoski 1986, Science -- Extinction data.
#' https://science.sciencemag.org/content/231/4740/833
#' 
#' Code from Section 4.1 (and Figure 2) in Toulis, P. and Bean, J. (2021). "Randomization Inference of Periodicity in Unequally Spaced Time Series with Application to Exoplanet Detection."
#' https://www.ptoulis.com/s/astro_main.pdf
#' 
#' Panos Toulis, ptoulis@chicagobooth.edu, 03/2021
#' 
rm(list=ls())
#' Derived the following data manually from Fig 2.
#' The following are cm measured from 0 of the corresponding axis.
Ti_plot = c(6.32, 6.20, 6.04, 5.90, 5.80, 5.70, 5.52, 5.40, 5.25, 5.12,
           4.90, 4.80, 4.63, 4.50, 4.35, 4.15, 4.00, 3.88, 3.70, 3.60,
           3.40, 3.28, 3.13, 3.00, 2.85, 2.70, 2.56, 2.40, 2.31, 2.20,
           2.05, 2.00, 1.74, 1.60, 1.40, 1.30, 1.20, 1.02, 0.90, 0.78,
           0.60, 0.48, 0.35, 0.25)

Yi_plot = c(1.54, 1.65, 4.40, 4.26, 3.29, 2.10, 1.67, 2.85, 2.00, 3.28,
            0.48, 1.00, 1.90, 1.50, 0.50, 1.58, 1.20, 1.39, 1.30, 1.13,
            2.11, 0.48, 0.98, 0.75, 0.82, 1.30, 0.52, 0.64, 0.80, 1.90,
            0.56, 0.80, 1.00, 3.48, 0.60, 0.45, 0.50, 0.92, 1.19, 0.49, 
            0.31, 0.30, 0.26, 0.80)

## Calibration to match data to figure in (R&S, 1986)
## 2.4cm in the plot => 100m
fctr = 248 / 251.667   # (calibrate wrt to paper.)
Ti = round(fctr * round(100*(Ti_plot / 2.4), 3), 3)  # time observations
## 4.3cm => 60% extinction
Yi = round(60 * (Yi_plot/4.3), 3)

# Data are (Ti, Yi)= (times, %extinctions.)
par(mar=rep(5, 4))
plot(-Ti, Yi, pch=20, type="b", labels=FALSE,  ylab="% extinction",
     xlab="time (million years from today)", cex.lab=1.8)
axis(1, at = c(-250, -200, -150, -100, -50, 0), 
     labels = c(250, 200, 150, 100, 50, 0), tick = T)


# Load radial velocity library.
source("rv_lib.R")

rv = data.frame(val=Yi, days=Ti, se=1) # SE does not matter here.
# all_P = seq(0.01, 62.5, length.out=1000000)
all_P = seq(12.5, 62.5, length.out=20000) # as in (R&S, 1986)

# LS periodogram for extinciotn data
out = lombe_scragle_fast(rv, all_P, v=T)

out$Phat  # periodogram peak

# Build Confidence Set.
Theta_hat = Build_ConfidenceSet(rv, all_P, null_samples = 250, time_budget_mins = 1)
# Confidence set is Theta_hat and includes all peaks.

#' Non-parametric test of Raup and Sepkoski (1986)
#' 
fit_g = function(C) {
  D = Inf
  t0 = 0
  val = NA
  while(D > 1e-4) {
    Ji = sapply(Pi, function(p) {
      cand_i = as.integer(C * (2*p - t0) / (2*C^2))
      f1 = (p - t0 - C*cand_i)^2
      f2 = (p - t0 - C*(cand_i+1))^2
      ifelse(f1 < f2, cand_i, cand_i+1)
    })
    t0_new = mean(Pi) - mean(Ji) * C
    D = abs(t0_new - t0)
    t0 = t0_new
    Pi_hat = t0 + Ji*C
    val = mean((Pi_hat - Pi)^2)
  }
  return(val)
}

all_C = seq(12.5, 60, length.out=50)  # periods to test.
fs = sapply(all_C, fit_g) # 
plot(all_C, fs, type="b", pch=20)


print(paste("Min C=", all_C[which.min(y)]))
