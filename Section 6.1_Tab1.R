#' Analysis of 4 exoplanets in Section 6 in Toulis, P. and Bean, J. (2021). "Randomization Inference of Periodicity in Unequally Spaced Time Series with Application to Exoplanet Detection."
#' https://www.ptoulis.com/s/astro_main.pdf
#' 
#' Panos Toulis, ptoulis@chicagobooth.edu, 03/2021
#' 
#' Here, we show that improved designs can resolve the identification issues.
#' For the α Centauri B example, define D_index=4.
#' For the Proxima Centauri example, define D_index=5.
#' The improved designs will be automatically defined according to Table 1, see parameters DELTA and NUM_NEW.
#' The code will produce results in new folders showing that the confidence set is sharp under the improved designs.
#' 
rm(list=ls())
source("rv_data.R") # loads lib.R

NUM_SAMPLES =  499
D_index = 5 # ESPRESSO (4= α Centauri B)
USE_DESIGN_A = T # A = uses delta (variation, fixed n); B = uses (n', new datapoints)

DELTA = c(0.18, 0.06)[D_index - 3]
NUM_NEW = c(137, 17)[D_index - 3]

# D = load_Dataset_Index(5)  # ESPRESSO data.
D = load_Dataset_Index(D_index)

rvData = D$RV
all_P = D$Periods
load("ALL_RESULTS_n100000.rdata")
ALL_RESULTS = as.data.frame(ALL_RESULTS)
A  = subset(ALL_RESULTS, dataset==D_index & (!reject))
rownames(A) = seq(1, nrow(A))


ls0 = lombe_scragle_fast(rvData, all_P, v=T)

#' Plots observed LS and null LS under signal Periodic(theta0).
#' 
#' Also calculates pvalue of P^-P(theta0) and saves file.
plot_overlap = function(theta0, suffix, RV, num_samples=NUM_SAMPLES) {
  ##
  print("---==-=-=-=====-=--=--=-=====-=-=-=--=-====-=-=-=-====-=")
  set.seed(123)
  print(paste("> Overlap for P0=", theta0, "Design=", suffix, "with", num_samples,"samples."))
  synth = Null_Periodograms_fast(RV, all_P, theta0=theta0, null_samples = num_samples)
  # Find overlap with observed GLS periodogram.
  dname =  sprintf("./DATA=%d", D_index)
  if(!dir.exists(dname)) {
    dir.create(dname)
  }
  fname = sprintf("%s/plot_%s(θ0=%.4f).png", dname, suffix, theta0)
  png(file=fname, width = 2000, height=800)
  ls_obs = lombe_scragle_fast(RV, all_P, v=T)
  for(j in 1:10) {
    lines(all_P, synth$all_Pf[j,], lwd=2, col=rgb(1, 0, 0, 0.25))
  }
  dev.off()
  id0= which.min(abs(all_P - theta0))
  tsamples = apply(synth$all_Pf, 1, function(pf) {
    max(pf) - pf[id0]
  })
  tobs = max(ls_obs$Pf) - ls_obs$Pf[id0]
  
  pval =  mean(tsamples >= tobs)
  print(paste(">> pval=", pval))
  fname = sprintf("%s/plot_%s(θ0=%.4f)_tstat.png", dname, suffix, theta0)
  
  png(file=fname, width = 2000, height=800)
  hist(tsamples, breaks=50, main=paste(suffix,"pval=", pval))
  abline(v=tobs, col="red")
  dev.off()
}

## Remove the problem with better design?
create_signal = function(ts_new, true_signal) {
  ts_new  = round(sort(ts_new), 3)
  y = rvData$val
  ei = y - ls0$lm_fit$yhat
  bhat = ls0$lm_fit$bhat
  
  P0 = true_signal
  signal = bhat[1] + bhat[2] * ts_new +  bhat[3] * cos(2*pi*ts_new/P0) +  bhat[4] * sin(2*pi*ts_new/P0) 
  # plot(ts_new, signal, pch=20, type="l")
  n  = length(ts_new)
  std_e = rvData$se
  if(n > length(ei)) {
    warning("Length of ts new > length(ei)")
    ei = sample(ei, size=n, replace=T)
    std_e = sample(rvData$se, size=n, replace=T)
  }
  rvData_new = data.frame(days=ts_new, 
                          val=signal + ei*sample(c(-1, 1), size=n, replace=T), 
                          se=std_e)
  # lombe_scragle_fast(rvData_new, all_P, v=T)
  return(rvData_new)
}

ts = rvData$days
n = length(ts)
tt = tail(order(A$`pval(Pf^-Pf(0))`), 2)
theta1d = A[tt[1],]$theta0
print(paste(">> PROBLEMATIC SIGNAL = ", theta1d))
print(A[tt[1],])

# plot_overlap(theta1d, suffix = "D0", RV=rvData, num_samples = NUM_SAMPLES)
print(paste("===>"))
print(paste(""))

if(USE_DESIGN_A) {
  # Design 1B: t' = t + 0.1 * Unif.
  print(paste("> Using Design (A)"))
  SUFFIX=paste("delta=", DELTA)
  ts_new = ts + DELTA*runif(n, min=-1, max=1)
  rvData_new = create_signal(ts_new, true_signal = ls0$Phat)
  plot_overlap(theta1d, suffix = SUFFIX, RV=rvData_new)
  new_theta = setdiff(Get_CandidatePeriods(rvData_new, all_P, threshold_power = .5), ls0$Phat)
  print(paste(">> Consider additional theta0 that exceed 50% of peak..Total",length(new_theta), "found."))

  for(th in new_theta) {
    plot_overlap(th, suffix =SUFFIX, RV=rvData_new)
  }
} else {
  print(paste("> Using Design (B)"))
  
  added_ts = sort(sample(seq(1, 100), size=NUM_NEW, replace=T) + 0.05*runif(NUM_NEW, min=-1, max=1))
  ts_new = c(ts, max(ts) + added_ts)
  print(paste("Old obs=", length(ts), 
              "Total new obs=", length(ts_new) - length(ts)))
  rvData_new = create_signal(ts_new, true_signal = ls0$Phat)
  plot_overlap(theta1d, suffix = SUFFIX, RV=rvData_new)
  new_theta = setdiff(Get_CandidatePeriods(rvData_new, all_P, threshold_power = .5), ls0$Phat)
  print(paste(">> Consider additional theta0 that exceed 50% of peak..Total",length(new_theta), "found."))
  
  for(th in new_theta) {
    plot_overlap(th, suffix =SUFFIX, RV=rvData_new)
  }
}



