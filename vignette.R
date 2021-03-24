#' This illustrates the use of code for Toulis, P. and Bean, J. (2021). "Randomization Inference of Periodicity in Unequally Spaced Time Series with Application to Exoplanet Detection."
#' https://www.ptoulis.com/s/astro_main.pdf
#' 
#' Panos Toulis, ptoulis@chicagobooth.edu, 03/2021
#' 
rm(list=ls())
# Example 1. Create radial velocity data frame.
set.seed(1)
n = 100 # sample size
var_days = 0.05 # daily variation
Tn = seq(1, n) + var_days*runif(n, min=-1, max=1) # observation times
theta_true = sqrt(2) # true periodicity

Yn = 1.5 * cos(2*pi*Tn/theta_true) + rnorm(n)  # synthetic RV measurements
rv = data.frame(val=Yn, days=Tn, se=rep(1, n))  # RV data frame.

#' Example 2. Generate periodogram for the RV data.
source("rv_lib.R") # Loads the main library.
all_P = 10^seq(-1, log(50, base=10), length.out=10000) # universe of periods (denoted Θ in paper)
ls0 = lombe_scragle_fast(rv, all_P, v=T) # LS periodogram
ls0$Phat # Periodogram peak should be at 3.14 days.

#' Example 3. Test whether the true period could be equal to sqrt(2) (i.e., the ground truth)
test = Test_Period_Theta0(rv, all_P, theta0=theta_true, null_samples=100)
test$pval # should be 0.38. We cannot reject that the true period is sqrt(2) even though the peak is at 3.14

#' Example 4. Build 99% confidence set.
ci = Build_ConfidenceSet(rv, all_P, time_budget_mins = 2)
# ci is an approximate 99% Confidence Set. 
#' In a full application we need to increase num_samples and the time budget.
#' See paper for details.

# Example 5. Real exoplanet detection problem: 51 Pegasi b
source("rv_data.R")
rvPeg = load_Dataset("51Pegb") # load RV dataset.
all_P = 10^seq(-1, 2.5, length.out=25000) # universe of periods (denoted Θ in paper)

ls0 = lombe_scragle_fast(rvPeg, all_P, v=T) # LS periodogram of 51Pegasi B
ls0$Phat # Periodogram peak should be at 3.14 days.

# Build 99% Confidence Set for 51 Pegasi b period. Takes 5 minute.
ci = Build_ConfidenceSet(rv, all_P, time_budget_mins = 5)

