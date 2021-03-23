#' Basic library of functions for the paper 
#' "Toulis, P. and Bean, J. (2021). Randomization Inference of Periodicity in Unequally Spaced Time Series with Application to Exoplanet Detection."
#' https://www.ptoulis.com/s/astro_main.pdf
#' 
print(">> LOADING rv_lib.R")
#' Checks whether the argument has the correct structure as a radial velocity data frame.
#' We define radial velocity as data.frame(days, val, se)
#' where days= observation times in JB days, val= radial velocity measurements, se=standard errors.
#' 
CHECK_rv = function(rvData) {
  stopifnot(class(rvData)=="data.frame")
  stopifnot(all(c("days", "val", "se") %in% names(rvData)))
  stopifnot(all(rvData$se >= 0))
}

#' Converts an RV data frame into a matrix of covariates. Our convention is rvData=data.frame(days, val, se, X1, X2, ...)
#'  where X1, X2, ... correspond to potentially additional covariates (e.g., stellar signals)
#' @return A nx(p-1) Matrix = [1, t, X1, X2, ...] where n=number of datapoints, p=number of columns in rvData.
get_X = function(rvData) {
  ts  = rvData$days
  X0 = cbind(b0=1, ts)  # default covariates.
  p = ncol(rvData) - 3 # no of covariate signals.
  if(p > 0) {
    X0 = cbind(X0, rvData[, seq(4, 3+p)])
  }
  return(as.matrix(X0))
}

#' Fit RV data.
#' @param rvData RV data (days, val, se), days=time of obervation, val=RV value, se=standard error.
#' @param add_P Vector (potentially empty) of periods which determines what periodic terms to include.
#' @return A LIST containing (sigma, bhat, yhat, logL, lm_fit) with information about the best weighted least squares fit (see below.)
#' 
fit_weighted = function(rvData, add_P=c()) {
  
  t0_start = proc.time()[3]
  
  CHECK_rv(rvData) # check whether rvData has correct structure.
  y = rvData$val; se = rvData$se; ts = rvData$days
  X0 = get_X(rvData) # baseline covariates.
  
  # Do we need to add periodic terms?
  if(length(add_P) > 0) {
    for(j in 1:length(add_P)) {
      X0 = cbind(X0, cos(2*pi*ts / add_P[j]), sin(2*pi*ts/ add_P[j]))
    }
  }
  
  # Quadratic loss function.
  loss = function(par, ret_logl=TRUE) {
    s  = exp(par[1])  # sigma parameter.
    v = s^2 + se^2
    fit = lm.wfit(X0, y, w = 1/v) # faster version than lm().
    # bhat = coef(fit)
    # res = y - arg_X %*% bhat
    res = fit$residuals
    df = nrow(X0) - ncol(X0)
    logL = sum(dnorm(res, sd=sqrt(v), log=T))
    # 
    if(!ret_logl) {
      yhat = y - res
      return(list(sigma=sqrt(v), bhat=as.numeric(fit$coefficients),
                  yhat=yhat,logL=logL, lm_fit=fit, 
                  R2=cor(y, yhat)^2,  df=df, 
                  dt=proc.time()[3] - t0_start))
    }
    return(-logL)
  }
  
  best = optim(par=c(log(mean(se))), fn=loss, method="L-BFGS-B", lower=-10, upper=10, control=list(factr=1e9))
  # Return optimized parameters.
  return(loss(best$par, ret_logl = FALSE))
}

#' Fit RV data using simple OLS ignoring the measurement standard error.
#' @param rvData RV data (days, val, se), days=time of obervation, val=RV value, se=standard error.
#' @param add_P Vector (potentially empty) of periods which determines what periodic terms to include.
#' @return A LIST containing (sigma, bhat, yhat, logL, lm_fit) with information about the best weighted least squares fit (see below.)
#' 
fit_OLS = function(rvData, add_P=c()) {
  #'
  t0_start = proc.time()[3]
  CHECK_rv(rvData)
  y = rvData$val; ts = rvData$days
  X0 = get_X(rvData)
  if(length(add_P) > 0) {
    for(j in 1:length(add_P)) {
      X0 = cbind(X0, cos(2*pi*ts / add_P[j]), sin(2*pi*ts/ add_P[j]))
    }
  }
  
  df = nrow(X0) - ncol(X0)
  fit = .lm.fit(X0, y) # lm(y ~ arg_X + 0)
  
  #' Save
  sigma = sqrt(sum(fit$residuals^2)/df) # summary(fit)$sigma
  bhat = as.numeric(coef(fit))
  yhat = y - fit$residuals
  logL = sum(dnorm(fit$residuals, sd=sigma, log=TRUE))
  return(list(sigma=sigma, bhat=bhat, 
              yhat=yhat, R2=cor(y, yhat)^2, df=df,
              logL=logL, lm_fit=fit, dt=proc.time()[3]-t0_start))
}

process_method = function(method) {
  fit_fn = NA  # fit_fn(y, X)
  if(method=="weighted_OLS") {
    fit_fn= fit_weighted
  } else if(method=="OLS") {
    fit_fn = fit_OLS # Speed things up.
  } else {
    stop(sprintf("Not supported method = %s", method))
  }
  return(fit_fn)
}

#' Lombe-Scragle periodogram.
#' Fits model
#'  y(t) = β0 + β1*(1,t,X) + β2*cos(2πt/θ) + β3*sin(2πt/θ) + ε.
#' 
#' @param rvData RV data to be analyzed.
#' @param all_P Vector of periods P to be searched (Θ in paper).
#' @param method Method to be used. One of {full, OLS}.
#' @return A LIST that includes (Pf, Phat) where Pf = vector of periodogram values; Phat = periodogram peak.
#' 
lombe_scragle = function(rvData, all_P, verbose=F, method="OLS", fap_level=0.01) { 
  t0_start = proc.time()[3]
  # 1a. Load (t, y)
  CHECK_rv(rvData)
  rv = rvData; ts = rv$days
  # X's
  
  # Fitting method.
  fit_fn = process_method(method)
  
  # Fit y  ~ 1.
  chi0 = fit_fn(rv)$logL
  
  # Scan all periods -- Periodogram.
  Pf = sapply(all_P, function(P) {
    chi1 = fit_fn(rv, add_P=P)$logL  # add a periodic signal
    # 1-chi1/chi0
    (chi1 - chi0) # LR statistic.
  })
  
  # find peak
  Phat = all_P[which.max(Pf)]
  # Get residual signal.
  best_fit = fit_fn(rv, add_P=Phat)$lm_fit
  Residual_rv = rv$val - best_fit$yhat

  nbreaks = length(all_P)
  #false-alarm
  a = seq(1, 100, length.out=10000)
  z = exp(nbreaks * log(pexp(a, rate=1/2)))  # F_n(1), F_n(1.01), ....F_n(100)
  z_q = a[which.min(abs(z-0.99))] / 2 # because we consider Δ logL as our statistic.
  
  if(verbose) {
    par(mfrow=c(1, 1))
    print(sprintf(">> Removed X signal. Method=%s", method))
    print(">> Plotting Periodogram.")
    options(scipen=5)
    plot(all_P, Pf, col="blue", type="l", lwd=2, main=sprintf("LS Periodogram (%d periods)", length(all_P)), 
         xlab="period", ylab="d(log L)", log="x", ylim=c(0, 1.3*max(Pf)))
    abline(v=Phat, col="red", lwd=0.5)
    text(Phat + 5, 1.8*max(Pf), sprintf("P^=%.2f", Phat), col="red", cex=.8)
    # false-alarm
    print(sprintf(">> %.1f%% percentile = %.1f", 100*(1-fap_level), z_q))
    abline(h=z_q, lty=3, col="gray")
    
    for(k in 1:10) {
      abline(v=Phat/k, col="red", lty=6, lwd=.5)
      abline(v=k*Phat, col="red", lty=6, lwd=.5)
      
    }
  }
  
  return(list(all_P=all_P, Pf=Pf, Phat=Phat, 
              exceed_fap=(max(Pf) >= z_q),
              res_rv=Residual_rv, lm_fit=best_fit,
              dt=proc.time()[3]-t0_start))
}

#' Fast implementation of LS periodogram. Based on Zechmeister and M. Kurster (2009).
#' https://www.aanda.org/articles/aa/abs/2009/11/aa11296-08/aa11296-08.html
#' See Equations (4)-(18)
#' @param rvData RV data to be analyzed.
#' @param all_P Vector of periods P to be searched (Θ in paper).
#' @return A LIST that includes (Pf, Phat) where Pf = vector of periodogram values; Phat = periodogram peak.
#' 
lombe_scragle_fast = function(rvData, all_P, verbose=F, fap_level=0.01) { 
  t0_start = proc.time()[3]
  # 1a. Load times and data.
  CHECK_rv(rvData)
  rv = rvData; 
  y = rv$val; ts = rv$days; se = rv$se
  y2 = y^2
  # Covariates X.
  w = (1/se^2); w = w / sum(w)
  Y = sum(w*y)
  
  Pf = sapply(all_P, function(P) {
    omega = 2*pi/P
    cosP = cos(omega*ts)
    sinP = sin(omega*ts)
    
    C = sum(w*cosP)
    S = sum(w*sinP)
    
    YY = sum(w*y2) - Y^2
    YC = sum(w*y*cosP) - Y*C
    YS =  sum(w*y*sinP) - Y*S
    CC = sum(w*cosP^2) - C^2
    SS = sum(w*sinP^2) - S^2
    CS = sum(w*sinP*cosP) - C*S
    D = CC*SS - CS^2
    
    ftr = (1/(YY*D))
    return(ftr*(SS*YC^2 + CC*YS^2 - 2*CS*YC*YS))
  })
  
  Phat = all_P[which.max(Pf)]
  fit_best = fit_weighted(rvData, add_P=Phat)
  
  # Get residual signal.
  nbreaks = length(all_P)
  #false-alarm
  a = seq(1, 100, length.out=10000)
  z = exp(nbreaks * log(pexp(a, rate=1/2)))  # F_n(1), F_n(1.01), ....F_n(100)
  z_q = a[which.min(abs(z-0.99))] / 2 # because we consider Δ logL as our statistic.
  
  if(verbose) {
    par(mfrow=c(1, 1))
    print(">> Plotting Periodogram.")
    options(scipen=5)
    par(mar=rep(5.5, 4))
    plot(all_P, Pf, col="blue", type="l", lwd=2, 
         main=sprintf("GLS Periodogram (%d periods)", length(all_P)), 
         xlab="period", ylab=expression(hat(A)[n] * "(" * theta * ")"),
                                   log="x", ylim=c(0, 1.3*max(Pf)),
         cex.lab=1.8, cex.axis=1.6)
    abline(v=Phat, col="red", lwd=0.5)
    text(Phat + 5, 1.8*max(Pf), sprintf("P^=%.2f", Phat), col="red", cex=.8)
    # false-alarm
    print(sprintf(">> %.1f%% percentile = %.1f", 100*(1-fap_level), z_q))
    print(sprintf(">> Best P^ =  %.2f", Phat))
    
    abline(h=z_q, lty=3, col="gray")
    
    for(k in 1:2) {
      #abline(v=Phat/k, col="black", lty=6, lwd=.5)
      #abline(v=k*Phat, col="black", lty=6, lwd=.5)
    }
  }
  
  return(list(all_P=all_P, Pf=Pf, Phat=Phat, 
              lm_fit=fit_best,
              exceed_fap=(max(Pf) >= z_q),
              dt=proc.time()[3]-t0_start))
}

#' Samples periodograms under the null that θ* = θ0.
#' @param rvData RV data to be analyzed.
#' @param all_P Vector of periods P to be searched (Θ in paper).
Null_Periodograms_fast = function(rvData, all_P, theta0, null_samples=100, verbose=F) {
  ##
  t0_start = proc.time()[3] # to measure execution time
  Yobs = rvData$val
  ts = rvData$days
  n = length(Yobs)
  
  # 1.  Fit y(t)  ~ X + Per(θ0, t) + ε(t).
  fit = fit_weighted(rvData, add_P=theta0)
  Yhat0 = fit$yhat
  # Residuals under H0  
  e0 = Yobs - Yhat0
  signs = c(rep(-1, n/2), rep(1, n/2))
  if(length(signs) < n) signs = c(signs, rep(1, n - length(signs)))
  
  # 2. Lombe-Scragle  periodogram.
  ls0 = function(y1) {
    rv_new = rvData
    rv_new$val = y1
    # rv_new = data.frame(val=y1, se=rvData$se, days=rvData$days)
    out = lombe_scragle_fast(rv_new, all_P)
    return(out)
  }
  
  all_Pf <<- matrix(0, nrow=null_samples, ncol=length(all_P))
  # Only assume sign symmetric errors.
  for(j in 1:null_samples) {
    e = sample(signs) * e0  # Sample under H0
    # e = sigma * rnorm(n)
    y_new = Yhat0 + e
    out = ls0(y_new)
    all_Pf[j, ] <<- out$Pf
  }
  
  out = list(all_Pf=all_Pf, all_P=all_P, theta0=theta0, 
             dt=proc.time()[3] - t0_start)
  return(out)
}

# Test decision at 0.1% level.
decide_test = function(tsamples, tobs, alpha=0.1/100) {
  R = length(tsamples)
  k = floor(R*(1-alpha))
  iOrd = order(tsamples)
  T_k = tsamples[iOrd[k]]
  R_plus = sum(tsamples > T_k)
  R_0 = sum(tsamples == T_k)
  threshold = (R*alpha - R_plus)/R_0
  if(tobs > T_k) { return(TRUE) 
  } else if(tobs == T_k) {
    return(runif(1) <= threshold)
  } else {
    return(FALSE)
  }
}
CHECK_decide = function() {
  
  sample_data = function(n) {
    B = rbinom(n, size=1, prob=0.6)
    X = 6*rbeta(n, shape1 = 2, shape2=3)# sample(u, size=n, replace=T)
    B * rep(0, n) + (1-B)*X # 0-inflation
  }
  X = sample_data(5000)
  
  y = replicate(20000, {
    Xobs = sample(X, 1)# sample_data(1) 
    decide_test(X, Xobs, alpha=0.1/100) 
  }) 
  100*mean(y)
}


#' Gets candidate periods to check.
#' @param rvData RV data to be analyzed.
#' @param all_P Vector of periods P to be searched (Θ in paper).
#' @description This method will return all periods in the periodogram that exceed "threshold_power"% of the peak.
#' @return Vector of periods.
Get_CandidatePeriods = function(rvData, all_P, 
                                threshold_power=0.2) {
  # Compute periodogram
  ls0 = lombe_scragle_fast(rvData, all_P)
  # Check which Period is peak.
  isPeak = sapply(seq(2, length(ls0$Pf)-1), function(i) {
    (ls0$Pf[i] > ls0$Pf[i-1]) & (ls0$Pf[i+1] < ls0$Pf[i])
  })
  isPeak = c(F, isPeak, F)
  
  # Filter out those peaks that are not high enough.
  J = intersect(which(isPeak), which(ls0$Pf > threshold_power*max(ls0$Pf)))
  
  return(all_P[J])
}


#' Tests H0: theta*=theta0
Test_Period_Theta0 = function(rvData, all_P, theta0, null_samples, ls0=NULL) {
  time_start  = proc.time()[3]
  
  if(is.null(ls0)) {
    ls0 = lombe_scragle_fast(rvData, all_P, v=F)
  }
  
  id0 = which.min(abs(all_P - theta0))
  
  Pf_Null = Null_Periodograms_fast(rvData, all_P, theta0, null_samples = null_samples)
  # Test statistic
  tstat = function(Pf) { 
    # take diff in peaks?
    best = tail(order(Pf), 2)
    max(Pf) - Pf[id0] 
  }
  
  # 4a. Observed test statistic.
  sobs = tstat(ls0$Pf)
  # 4b. Null distribution,
  S_null = apply(Pf_Null$all_Pf, 1, function(pf) { tstat(pf) })
  
  list(sobs=sobs, svals=S_null, dt=proc.time()[3] - time_start, 
       pval=mean(S_null >= sobs))
}

#' Main method that builds confidence set Θ_{1-α} in Equations (7) and (12) of paper.
#' @param rvData RV data to be analyzed.
#' @param all_P Vector of periods P to be searched (Θ in paper).
#' @param alpha Confidence level.
#' @param null_samples Number of randomization samples (variable R in paper).
#' @param time_budget_mins How many mins to spend in computation. 
#' 
#' @description The method automatically calculates how many periods to consider based on the desired time budget.
#' 
Build_ConfidenceSet = function(rvData, all_P, alpha=0.01, 
                               null_samples=100, 
                               time_budget_mins=1, retPval=FALSE) {

  print("[ SET IDENTIFICATION OF HIDDEN PERIODICITY. ]")
  print("[ -=-====-----=-=-=-===--------=-==-====--= ]")
  print(paste(">> Time budget=", time_budget_mins, " mins."))
  
  # Return result.
  pvals_m = matrix(0, nrow=0, ncol=4)
  colnames(pvals_m) = c("theta0", "<=sobs", ">=sobs", "==sobs")
  
  # 1. Get candidate thetas
  # How many thetas to consider given the budget. 
  # 1a. Calculate Time/periodogram (should be ~0.2sec on conventional laptop)
  print(paste(">> Estimating time/periodogram calculation rate."))
  t0 = proc.time()[3]
  out = replicate(10, { lombe_scragle_fast(rvData, all_P, v=F) })
  RATE = (proc.time()[3]-t0) / 10   # secs/periodogram
  
  all_threshold = seq(0.1, 0.8, length.out=50)
  z = sapply(all_threshold, function(i) length(Get_CandidatePeriods(rvData , all_P, threshold_power = i)))
  best_threshold = all_threshold[which.min(abs(z*null_samples*RATE/60 - time_budget_mins))]  # time in mins
  
  thetas = Get_CandidatePeriods(rvData , all_P, threshold_power = best_threshold)
  print(paste(">> Will test", length(thetas), "candidate periods."))
  if(null_samples < 100) {
    warning("Number of null samples should be larger than 100. e.g., 1,000")
  }
  # 1b. Observed periodogram
  ls0 = lombe_scragle_fast(rvData, all_P, v=F)
  cnt = 0
  time_start = proc.time()[3]
  for(theta0 in thetas) {
    cnt = cnt + 1
    # H0: theta*=theta0
    id0 = which.min(abs(all_P - theta0))
    if(length(id0)==0) {
      stop(paste(theta0, " is not in the Theta (universe set of periods)"))
    }
    # 2. For each theta0, sample from the null.
    print(paste(">> Sampling", null_samples, "periodograms. Please wait.."))
    Pf_Null = Null_Periodograms_fast(rvData, all_P, theta0, null_samples = null_samples)
    
    # dim(Pf_Null$all_Pf)   #100 x 10,000 should be.
    # 3. Test statistic.
    tstat = function(Pf) { 
      max(Pf) - Pf[id0] 
      # Pf[id0] / mean(Pf)
      # max(Pf) / mean(Pf)
    }
    
    # 4a. Observed test statistic.
    sobs = tstat(ls0$Pf)
    # 4b. Null distribution,
    S_null = apply(Pf_Null$all_Pf, 1, function(pf) { tstat(pf) })
    
    print(paste(">> Sampling took ", round(Pf_Null$dt, 2), "seconds. Now saving p-values.."))
    par(mar=rep(2.5, 4))
    par(mfrow=c(2, 1))
    
    hist(S_null, breaks=30, main=paste("theta0=", round(theta0, 2), "sobs=",round(sobs, 2)))
    abline(v=sobs, col="red", lwd=2)
    
    plot(all_P, ls0$Pf, log="x", col="blue", type="l", lwd=1.5)
    for(j in 1:10) {
      lines(all_P, Pf_Null$all_Pf[j, ], col=rgb(1, 0, 0, alpha=0.1))
    }
    
    
    pvals_m = rbind(pvals_m, c(theta0, mean(S_null <= sobs), 
                               mean(S_null >= sobs),
                               mean(S_null == sobs)))
    print(pvals_m)
    print(paste("> [", cnt,"/", length(thetas)," ", round(proc.time()[3]-time_start, 2), "secs.] Last test, H0: theta*=", theta0))
    
  }
  # Return p-values.
  if(retPval) { return(pvals_m) 
  } else {
    incl = which(pvals_m[,3] > alpha)
    return(as.numeric(pvals_m[incl, 1]))
  }
}

#' Non-parametric test of periodicity.
#' See Section 6.2 of the paper.
Test_Period_np = function(rvData, all_P, theta0, null_samples=100, verbose=F) {
  #  theta0 = 12;   num_splits = 100
  print(paste("Testing theta0=", theta0))
  ts = rvData$days # time stamps.
  yobs = rvData$val
  
  id0 = which.min(abs(all_P - theta0))
  
  test_stat = function(y) {
    rv_new = rvData
    rv_new$val = y
    out = lombe_scragle_fast(rv_new, all_P, v=F)
    max(out$Pf) - out$Pf[id0]
  }
  
  # observed test statistic.
  tobs = test_stat(yobs)
  
  # Resampling.
  rem = theta0 * ((ts/theta0) - floor(ts/theta0))  # remainders.
  # Split [0, theta0] in num_splits.
  ## Choose num_splits
  print(paste(">> Choosing num_splits..."))
  all_splits = seq(1e-3, 0.05, length.out=100) 
  w = sapply(all_splits, function(s) {
    num_splits = floor(diff(range(ts)) / s)
    idx = as.numeric(cut(rem, breaks = num_splits))
    J = as.numeric(table(idx))
    log(prod(J^J))
  })
  I = which(w > log(100))
  e_split = NA
  WARNING=FALSE
  if(length(I) ==0) { 
    e_split = max(all_splits)
    WARNING=TRUE
    warning("Cannot test reliably.")
  } else {
    e_split = min(all_splits[I])
  }

  num_splits = floor(diff(range(ts)) / e_split)  # t = t'  if  |t-t'| < range(T) / num_splits = ε.
  print(paste(">> Delta(t, t')=",round(24*e_split, 2)," hours. Will split observation times into num_splits=",num_splits))
  
  idx = as.numeric(cut(rem, breaks = num_splits))
  # We can permute within the same class.
  # 1    2    3    4    5
  # 1   17  74     82   82    # swap 4 and 5.
  cl = unique(idx)
  time0 = proc.time()[3]
  tvals = replicate(null_samples, {
    # Sample new data: ynew
    ynew = yobs
    for(icl in cl) {
      idx_icl = which(idx==icl)  # which indices correspond to class icl.
      # permute within class.
      if(length(idx_icl) > 1) {
        ynew[idx_icl] = ynew[sample(idx_icl, replace=T)]
      }
    }
    test_stat(ynew)
  })
  
  if(verbose) {
   print(paste("Time elapsed...", round(proc.time()[3] - time0, 2), "secs."))
   hist(tvals, breaks=40)
   abline(v=tobs, col="red", lwd=3)
  }
  pval1 = mean(tvals >= tobs)
  pval2 = mean(tvals <= tobs)
  pval0 = mean(tvals==tobs)
  return(list(left=pval2, right=pval1, on=pval0, 
              obs=tobs, tvals=tvals, num_tvals=length(tvals),
              hasWarning=WARNING, 
              dt= round(proc.time()[3] - time0, 2)))
}


Build_ConfidenceSet_np = function(rvData, all_P, 
                                      null_samples=100, 
                                      time_budget_mins=1) {
  
  print("[ NONPARAMETRIC SET IDENTIFICATION OF HIDDEN PERIODICITY. ]")
  print("[ -=-====-----=-=-=-===--------=-==-====--= ]")
  print(paste(">> Time budget=", time_budget_mins, " mins."))
  
  # Return result.
  cols = c("theta0", "<=sobs", ">=sobs", "==sobs", "warning")

  pvals_m = matrix(0, nrow=0, ncol=length(cols))
  colnames(pvals_m) = cols
  
  # 1. Get candidate thetas
  # How many thetas to consider given the budget. 
  # 1a. Calculate Time/periodogram (should be ~0.2sec on conventional laptop)
  print(paste(">> Estimating time/periodogram calculation rate."))
  t0 = proc.time()[3]
  out = replicate(10, { lombe_scragle_fast(rvData, all_P, v=F) })
  RATE = (proc.time()[3]-t0) / 10   # secs/periodogram
  
  all_threshold = seq(0.1, 0.8, length.out=50)
  z = sapply(all_threshold, function(i) length(Get_CandidatePeriods(rvData , all_P, threshold_power = i)))
  best_threshold = all_threshold[which.min(abs(z*null_samples*RATE/60 - time_budget_mins))]  # time in mins
  
  thetas = Get_CandidatePeriods(rvData , all_P, threshold_power = best_threshold)
  print(paste(">> Will test", length(thetas), "candidate periods."))
  if(null_samples < 100) {
    warning("Number of null samples should be larger than 100. e.g., 1,000")
  }
  # 1b. Observed periodogram
  ls0 = lombe_scragle_fast(rvData, all_P, v=F)
  cnt = 0
  time_start = proc.time()[3]
  for(theta0 in thetas) {
    cnt = cnt + 1
    # H0: theta*=theta0
    id0 = which.min(abs(all_P - theta0))
    if(length(id0)==0) {
      stop(paste(theta0, " is not in the Theta (universe set of periods)"))
    }
    # 2. For each theta0, sample from the null.
    # print(paste(">> Sampling", null_samples, "periodograms. Please wait.."))
    out = Test_Period_Theta0_np(rvData, all_P, theta0, null_samples = null_samples)
    S_null = out$tvals
    sobs = out$obs
    
    abline(v=sobs, col="red", lwd=2)
    pvals_m = rbind(pvals_m, c(theta0, mean(S_null <= sobs), 
                               mean(S_null >= sobs),
                               mean(S_null == sobs), 
                               out$hasWarning))
    print(pvals_m)
    print(paste("> [", cnt,"/", length(thetas)," ", round(proc.time()[3]-time_start, 2), "secs.] Last test, H0: theta*=", theta0))
    
  }
  # Return p-values.
  return(pvals_m)
}
