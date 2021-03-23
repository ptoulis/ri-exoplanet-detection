#' Code from Section 4.2 in Toulis, P. and Bean, J. (2021). "Randomization Inference of Periodicity in Unequally Spaced Time Series with Application to Exoplanet Detection."
#' https://www.ptoulis.com/s/astro_main.pdf
#' 
#' Panos Toulis, ptoulis@chicagobooth.edu, 03/2021
#' 
#' Simulation based on Hall and Li (2000)
#' 
rm(list=ls())
source("rv_lib.R")

h = function(x) 1 - cos(2*pi*x) 
theta_true = sqrt(2) # true period

#' Observation designs. Correspond to the three choices in Section 4.2.
sample_days = function(n, design_no=2) {
  # Design 1.
    # Design 2 (periodic)
    theta_f=1
    all_x = seq(0, n, length.out=10000) 
    if(design_no==1) {
      f1 = function(x) (x <= 1 & x >=0.6)  # it is night?
    } else if(design_no==2) {
      # more probability near midnight
      f1 = function(x) max(c(0, -pi*sin(2*pi*(x))* (x <= 1 & x >=0)))
    } else if(design_no==3) {
      # +- .5hr around midnight.
      f1 = function(x) (abs(x-0.75) < 0.0416) * (x <= 1 & x >=0)
    } else {
      stop("Design not supported.")
    }
    
    #f1 = function(x) abs((x - 3/4) < 0.001)
    f = sapply(all_x, function(x)  sum(f1(x - seq(0, n-1)*theta_f)))
    f = f/sum(f)  # make into density
    
    return(sort(sample(all_x, size=n, replace=F, prob = f)))
}

#' Sample radial velocity data.
sample_y = function(days, sigma) {
  n = length(days)
  y = h(days/theta_true) + rnorm(n, sd=sigma)
  return(y)
}

Figure2 = function(nsamples, budget) {
  
  set.seed(43)  # problematic for Sim 2
  # Simulation 1 or 2.
  design_no=2
  n = 200
  sigma = 1.5
  days = sample_days(n, design_no)
  y = sample_y(days, sigma)
  all_P = 10^seq(-1, 2, length.out=10000)
  # lines(days, h(days/theta0), col="red", type="b")
  rvData = data.frame(days=days, val=y, se=sigma)
  ls0 = lombe_scragle_fast(rvData, all_P, v=T)
  
  Test_Period_Theta0(rvData, all_P, theta0 = theta_true, null_samples=100)
  # out1 = Test_Period_Theta0_np(rvData, all_P, theta0=ls0$Phat, null_samples = 100, verbose = T)
  # out = Test_Period_Theta0_np(rvData, all_P, theta0=theta_true, null_samples = 100, verbose = T)
  # 9.6471368626
  # 12.7177546117
  result = Test_Period_Theta0(rvData, all_P, 9.6471368626, null_samples = 100)
  result =  Test_Period_Theta0(rvData, all_P, 12.7177546117, null_samples = 1000)

  #confset_np = Build_ConfidenceSet_theta0_np(rvData, all_P, null_samples = 200, time_budget_mins = 10)
  out = Build_ConfidenceSet_theta0(rvData, all_P, null_samples = nsamples,
                                   time_budget_mins = budget)
  
  save(out, file="results_Confidence_Sets.rda")
}


Figure2b = function(nreps=100, nsamples=100) {
  
  # Simulation 1 or 2.
  design_no=2
  n = 200
  sigma = 1.5
  pvals = c()
  all_P = 10^seq(-1, 2, length.out=10000)
  mode = c()
  same_alt = c()
    
  for(j in 1:nreps) {
    days = sample_days(n, design_no)
    y = sample_y(days, sigma)
    # lines(days, h(days/theta0), col="red", type="b")
    rvData = data.frame(days=days, val=y, se=sigma)
    ls0 = lombe_scragle_fast(rvData, all_P, v=T)
    
    cand = Get_CandidatePeriods(rvData, all_P)
    
    id0 = which.min(abs(cand- theta_true))
    
    theta0 = cand[id0]
    id1 = which.min(abs(all_P - theta0))
    isMode = (ls0$Pf[id1] == max(ls0$Pf))
    # Test: H0:theta*=theta0
    out = Test_Period_Theta0(rvData, all_P, theta0, null_samples=nsamples)
    pvals = c(pvals, out$pval)
    mode = c(mode, isMode)
    
    cover95 = ifelse(abs(theta0 - theta_true) > 0.2, 0, mean(pvals >= 0.05))
    cover99 = ifelse(abs(theta0 - theta_true) > 0.2, 0, mean(pvals >= 0.01))
    print(paste(">> ", j,"/", nreps,",", "dt=", round(out$dt, 1), "secs.",
                "95% cover=",round(100*cover95, 1),"%  -  99% cover=", round(100*cover99, 1),"%",
                " - id0=Pf mode in", round(100*mean(mode), 1),"% of cases."))
  }
  return(pvals)
}

## Sampling distribution of periodogram peak.
## Simulation based on Hall and Li (2000)
Figure1 = function(nsamples) {
  
  # theta_f = 7 # observations per week
  all_n = c(50, 100, 200, 500)
  all_designs = c(1,2, 3)
  all_P = 10^seq(-1, 2, length.out=10000)

  RESULTS = matrix(0, nrow=0, ncol=3)
  colnames(RESULTS) = c("n", "design", "theta_hat")
  noise_level = 1.5
  
  for(n in all_n) {
    for(j in all_designs) {
      set.seed(43)
      days = sample_days(n, j) # sample T^n
      print(paste("theta^ distribution for n=", n,"design=",j))
      peaks = replicate(nsamples, {
        y = sample_y(days, sigma=noise_level)
        # lines(days, h(days/theta0), col="red", type="b")
        rvData = data.frame(days=days, val=y, se=noise_level)
        ls0 = lombe_scragle_fast(rvData, all_P, v=F)
        ls0$Phat
      })
      hist(peaks, breaks=50, main=sprintf("peaks n=%d design=%d", n, j))
      RESULTS = rbind(RESULTS, cbind(rep(n, nsamples), rep(j, nsamples), peaks))
      # print(RESULTS)
    }
  }
  
  out = as.data.frame(RESULTS)
  save(out, file="sampling.rda")
  A = subset(out, n==50 & design==2)
  hist(A$theta_hat, breaks=100)
  A = subset(out, theta_hat <=5)
  library(ggplot2)
  g = ggplot(data=A)
  g = g + geom_histogram(aes(x=theta_hat, y=..density..), binwidth=0.2)
  design.labs = c("(i)", "(ii)", "(iii)")
  names(design.labs) = c("1", "2","3")
  g = g +   facet_grid(design ~ n, 
                       labeller=labeller(design=design.labs)) 
  g = g+theme(text = element_text(size=20),
              axis.text.x = element_text( hjust=1))
  g = g+ labs(x=expression(paste(theta,"^")))
  plot(g)
  
  return(as.data.frame(RESULTS))

}

