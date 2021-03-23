#' Load all RV datasets in one object.
#'
rm(list=ls())
source("exoplanet_functions.R")
options(digits=10)

create_GenFunction= function(arg_rv, X, all_P) {
  #' Generates the synthetic sampling function.
  #' rv = RV data
  #' X = activity data.
  #' all_P = period search space.
  #'    
  CHECK_rv(arg_rv)
  # 1. activity component
  rv = remove_signal(arg_rv, X)
  A_t = arg_rv$val - rv$val  # activity component.
  
  # 2. planet component
  pg = new_LS(rv, all_P, verbose=FALSE)
  phat = pg$Phat;  t1 = rv$days
  X_per = cbind(1, cos(2*pi*t1/phat), sin(2*pi*t1/phat))  # planet signal
  fit = lm(rv$val ~ X_per + 0)
  bhat = coef(fit)  # planet  coefficients.
  sigma_hat = summary(fit)$sigma
  
  rm(rv)
  
  fn = function(P0) {
    
    Xnew = cbind(1, cos(2*pi*t1/P0), sin(2*pi*t1/P0))
    sigma = sqrt(sigma_hat^2 + arg_rv$se^2)
    y_synth = A_t + (Xnew %*% bhat) + rnorm(length(t1), sd=sigma)
    
    rv_ret = arg_rv
    rv_ret$val = y_synth  # same t, se but different "val".
    return(rv_ret)
  }
  
  return(fn)
}

get_rv_blankSep = function(out) {
  rv = t(apply(out,  1, function(row) {
    s = as.numeric(strsplit(row, split=" ")[[1]])
    s = s[!is.na(s)]
  }))
  colnames(rv) = c("days", "val", "se")
  return(as.data.frame(rv))
}

load_basic = function(fname) {
  if(!file.exists(fname)) {
    stop(sprintf("Filename for %s does not exist.", fname))
  }
  out = read.csv(fname, header=F, stringsAsFactors = F)
  rv = get_rv_blankSep(out)
  rv$days = rv$days - min(rv$days)  ## offset  days to ease LS.
  ## Main Export.
  Xrv = cbind(1, rv$days)  # only use time for 51Pegb or GJ436b 
  ## Add sampling function
  list(rv=rv, X=Xrv)
}

load_HD85512b = function() {
  out = read.csv("data/HD85512b/HD85512b.txt", header=F, stringsAsFactors = FALSE)
  out = t(apply(out,  1, function(row) {
    s = as.numeric(strsplit(row, split=" ")[[1]])
    s[!is.na(s)]
  }))
  
  rv = data.frame(days=out[,1], val=out[,2], se=out[,3])
  rv$days = rv$days -  min(rv$days)
  Xrv = out[, c(4, 6, 8)]
  
  list(rv=rv, X=Xrv)
}

load_damasso2020 = function() {
  print("> Loading Damasso 2020 data...")
  fname = "data/Damasso2020/GL551_UVES_binned.txt"
  ## 
  out = read.csv(fname, header=T, stringsAsFactors = F)
  rvUVES = get_rv_blankSep(out)
  fname = "data/Damasso2020/Proxima_pre2017_HARPS.txt"
  out = read.csv(fname, header=T, stringsAsFactors = F)
  rvHARPSpre = get_rv_blankSep(out)
  rvHARPSpre$days = rvHARPSpre$days + 2400000	
  fname = "data/Damasso2020/Proxima_2017_HARPS.txt"
  out = read.csv(fname, header=T, stringsAsFactors = F)
  rvHARPSpost = get_rv_blankSep(out)
  rvHARPSpost$days = rvHARPSpost$days + 2400000	
  rvHARPS = rbind(rvHARPSpre, rvHARPSpost)
  ##
  ## same start
  rvHARPS$days =  rvHARPS$days - 2.45*1e6
  rvUVES$days  = rvUVES$days - 2.45*1e6
 ## Combined dataset.
  rv = rbind(rvUVES, rvHARPS)

  ## Ha: spectroscopy
  out =  read.csv("data/Damasso2020/halpha_harps_unbinned.txt")  # HARPS, unbinned
  Ha = get_rv_blankSep(out)
  ## Need to bin observations
  new_Ha = matrix(0, nrow=1, ncol=ncol(Ha))
  colnames(new_Ha) = c("days", "val", "se")
  new_Ha = as.data.frame(new_Ha)
  new_Ha[1,] = Ha[1,]
  for(i in 2:nrow(Ha)) {
    # print(i)
    last = nrow(new_Ha)
    if(Ha$days[i] - new_Ha$days[last] < 0.5) {
      new_Ha[last,] = .5*(new_Ha[last,] + Ha[i,]) # merge
    } else {
      new_Ha[last+1,] = Ha[i,]
    }
  }
  stopifnot(all(diff(new_Ha$days) > .5))
  Ha = new_Ha
  
  # UVES Ha index.
  Ha_uves =  read.csv("data/Damasso2020/UVES_HALPHA.txt", sep = ",", header=F)
  colnames(Ha_uves) = c("days", "val", "se")
  Ha_uves$days =  Ha_uves$days - 2.45*1e6
  Ha_uves$val = Ha_uves$val - 1.26 #  as described in the paper.
  # All Ha
  Ha = rbind(Ha_uves, Ha)
  
  ## Merge Ha data with RV.
  abs_d = sapply(1:nrow(rv), function(i) { 
    j = which.min(abs(Ha$days - rv$days[i])); 
    a= Ha[j, ]; 
    abs(rv$days[i]-a$days) 
  })
  rmv_id = which(abs_d > 0.5)
  # Update rv dataset by removing days for which there is no Ha data.
  rv = rv[-rmv_id, ]
  
  new_Ha = matrix(0, nrow=1, ncol=ncol(Ha))
  colnames(new_Ha) = c("days", "val", "se")
  new_Ha = as.data.frame(new_Ha)
  for(i in 1:nrow(rv)) {
    j = which.min(abs(Ha$days - rv$days[i]))
    new_Ha[i,] = Ha[j,]
  }
  Ha = new_Ha
  # z = sapply(1:nrow(rv), function(i) { a= find_close(i); abs(rv$days[i]-a$days) })
  ## Now rv, Ha  have the same length and close time index.
  # plot(Ha$days, rv$days, type="l")
  # abline(0, 1, col="red")
  print(sprintf(">> Mean abs diff in time %.5f", mean(abs(Ha$days-rv$days))))
  
  ## Create X matrix.
  X = cbind(1, rv$days, rv$days^2, Ha$val)
  
  return(list(rv=rv, X=X))
}

load_dumusque2012 = function() { 
  print("Dumusque data.")
  fname = "datasets/Dumusque2012.txt"
  astro = read.csv(fname, sep = "\t", stringsAsFactors = FALSE)
  astro = astro[-c(1),]
  for(cc in colnames(astro)) {
    astro[[cc]] = as.numeric(astro[[cc]])
  }
  
  rv = data.frame(val=1000*astro$RV..km.s., se=1000*astro$RV.error..km.s.,
                  days=astro$jdb..days.)
  rv$days = rv$days - min(rv$days)  ## offset  days to ease LS.
  
  # Determine year of  observation
  block  = which(diff(rv$days) > 200)
  start = 1
  years = rep(NA, length(rv$days))
  for(b  in 1:length(block)) {
    bb = block[b]
    years[seq(start, bb)] = 2008 + b - 1
    start = bb + 1
  }
  years[seq(start,  length(years))] = 2011
  
  ## 1. Binary signature.
  BinX = cbind(1, rv$days, rv$days^2)
  
  # 2. Magnetic cycle
  MagX = matrix(astro$Low.pass.filter.log.R.hk., ncol=1)
  
  ## 3. Rotational Periods
  
  
  RotPer = c(39.76, 37.83, 36.75)
  what_to_fit = list()
  what_to_fit[[1]] = c(1, 2)
  what_to_fit[[2]] = c(1, 3, 4)
  what_to_fit[[3]] = c(1, 2, 3)
  
  rotlist = list()
  cols = c()
  Nt = length(years)
  RotX = matrix(0, nrow=length(years), ncol=0)
  # For all years 2009, 2010, 2011
  for(i in 1:3) {
    yr = 2008 + i
    Ind = (years==yr)
    P = RotPer[i]
    T1 = rv$days
    for(j in what_to_fit[[i]]) {
      Pj = P/j
      x1 = cos(2*pi*T1 / Pj) * Ind
      x2 = sin(2*pi*T1 / Pj) * Ind
      RotX = cbind(RotX, x1)
      RotX = cbind(RotX, x2)
      cols = c(cols, c(paste(yr, "cosP/",j, sep=""), 
                       paste(yr, "sinP/", j, sep="")))
      
    }
  }
  colnames(RotX) = cols
  
  # 
  Xrv = cbind(BinX, MagX, RotX)
  return(list(rv=rv, X=Xrv))
}

load_escude2016 = function() {
  spec_files = list.files("data/Escude2016/Spectroscopy", full.names = T)
  spec_names = c("Ha", "HE", "m2", "m3", "s-index")
  Xlist = list()
  stopifnot(length(spec_files) > 0)
  warning(">> Removing UVES file for Escude 2016.")
  spec_files = setdiff(spec_files, c("data/Escude2016/Spectroscopy/UVES_HALPHA.csv"))
  for(i in 1:length(spec_files)) {
    file = spec_files[i]
    out = read.csv(file, header=F)
    colnames(out) = c("days", "val", "se", "pre2016")
    Xlist[[i]] = out
  }
  
  # 2. Load RV data.
  rvHARPS = read.csv("data/Escude2016/RV/RVTERRA.csv", header=F)
  colnames(rvHARPS) = c("days", "val", "se", "pre2016")
  rvHARPS_post = subset(rvHARPS, pre2016==2)
  rvHARPS_pre = subset(rvHARPS, pre2016==1)
  rvHARPS_pre$pre2016 = NULL
  rvHARPS_post$pre2016 = NULL
  
  
  # 3. Create the "X"
  Xspec_pre = matrix(0, nrow=nrow(rvHARPS_pre), ncol=0)
  Xspec_post = matrix(0, nrow=nrow(rvHARPS_post), ncol=0)
  for(j in 1:length(Xlist)) {
    Xspec_pre = cbind(Xspec_pre, subset(Xlist[[j]], pre2016==1)$val)
    Xspec_post = cbind(Xspec_post, subset(Xlist[[j]], pre2016==2)$val)
  }
  colnames(Xspec_pre) = spec_names
  colnames(Xspec_post) = spec_names
  
  
  
  # FINALLY: Difference out time.
  t0 = min(rvHARPS_post$days)
  rvHARPS_post$days = rvHARPS_post$days - t0
  
  ## All covariates pre-2016
  Xall_post = as.matrix(cbind(days=rvHARPS_post$days, Xspec_post))
  
  ##
  rv =  rvHARPS_post
  Xrv = Xall_post
  return(list(rv=rv, X=Xrv))
}

#' Add RV data entrty.
#'
gen_rvObj = function(name, all_P=NULL) {
  if(is.null(all_P)) {
    all_P = 10^seq(-1, 3, length.out=4000)  # classical.
  }
  # all_P = 10^seq(-1, 3, length.out=nPeriods)  # classical.
  
  dat = NA ## needs to be list(rv, X)
  dname = sprintf("./data/%s/", name)
  
  if(!file.exists(dname)) {
    stop(sprintf("Data %s does not exist...!", dname))
  } else {
    
    print(sprintf("Loading data for %s", dname))
    if(name == "51Pegb" || name=="GJ436b" || name=="Damasso2017") {
      ## 
      fname = sprintf("data/%s/%s.txt", name, name)
      dat = load_basic(fname)
      #
    } else if(name=="Dumusque2012") {
      # Dumusque 2012 data.
      dat = load_dumusque2012()
      #
    } else if(name=="HD85512b") {
      ##
     dat = load_HD85512b()
      ##
    } else if(name=="Escude2016")  {
     dat = load_escude2016() 
    } else if(name=="Damasso2020") {
      dat = load_damasso2020()
    }
  }
  
  # >>>  At this point should produce (rv, Xrv, all_P)
  rv = dat$rv; Xrv = dat$X
  print(sprintf("> Period search between [%.3f, %.3f] days over %d total periods.", min(all_P), max(all_P), length(all_P)))
  
  print(sprintf("> Finding planet for Data=%s..", name))
  Phat_0 = find_planet(rv, Xrv, all_P, verbose = T) # original RV analysis
  print("> Adding gen_synth function..")
  
  # fn that samples synthetic  RV.
  fn = create_GenFunction(rv, Xrv, all_P)
  
  find_planet_synth = function(P0) {
    synth_rv = fn(P0)  # sample synthetic rv.
    phat = find_planet(synth_rv, Xrv, all_P, verbose=FALSE) #  synthetic RV analysis
    phat  # return fitted period.
  }
  
  list(rv=rv, X=Xrv, 
       find_planet_synth=find_planet_synth, 
       Phat=Phat_0, 
       all_P=all_P)
  
}

load_all_Data = function() {
  rvList = list()
  data_files = list.dirs("./data/", full.names = FALSE, recursive = F)
 cnt = 1
  for(Df in data_files) {
    t0 = proc.time()[3]
    if(Df=="Damasso2020") {
      all_P = 10^seq(log(2, base=10), log(3000, base=10), length.out=6000)
    } else {
      all_P = 10^seq(-1, 3, length.out=4000)
    }
    rvList[[cnt]] = gen_rvObj(Df, all_P)
    cnt = cnt + 1
    print(sprintf(">>> Loading time for %s =  %.3f secs..", Df, proc.time()[3] - t0)) 
    print("")
  }
  
  saveRDS(rvList, file="rvList.rdata")
}


require(nortest)
require(ks)
require(e1071)
D = rnorm(5000)

tstat = function(e_arg)  { 
  # c(mean(u^3), mean(u^4))
 # c(shapiro.test(e_arg)$statistic, lillie.test(e_arg)$statistic)
  # c(shapiro.test(e_arg)$statistic, pearson.test(e_arg)$statistic)
  # as.numeric(c(shapiro.test(e_arg)$statistic, ks.test(e_arg/sd(e_arg), D)$statistic))
  # c(shapiro.test(e_arg)$statistic, max(e_arg) / sd(e_arg))
  
  # c(lillie.test(e_arg)$statistic, ad.test(e_arg)$statistic)
  # c(skewness(e_arg), kurtosis(e_arg), shapiro.test(e_arg)$statistic)
  c(shapiro.test(e_arg)$statistic, ks.test(e_arg/sd(e_arg), D)$statistic)
}

gen_H0 = function(n) {
  vals = t(replicate(10000, {
    tstat(rnorm(n))
  }))
  
  m = colMeans(vals)
  S_1 = solve(cov(vals))
  
  t_comb = function(v) { 
    u = matrix(as.numeric((v - m)), nrow=1)
    as.numeric(u %*% S_1 %*% t(u))
  }
  
  if(ncol(vals)==2) {
    all_t  = expand.grid(t0=seq(min(vals[,1]), max(vals[,1]), length.out=500), 
                         t1=seq(min(vals[,2]), max(vals[,2]), length.out=500))
    z_vals = apply(vals, 1, t_comb)
    q1 = quantile(z_vals, 0.95)
    
    z_all = apply(all_t, 1, t_comb)
    i = which(abs(z_all -  q1) < 1e-1)
    plot(vals[,1], vals[,2], pch=20, cex=0.1)
    points(all_t[i,1], all_t[i,2], pch=20, cex=0.03, col="red")
    points(m[1], m[2], col="purple", pch="o", cex=2)
  } else {
    pairs(vals, pch=20, cex=0.05)
  }
  
  H0 = apply(vals, 1, t_comb)

  return(list(H0=H0, t_comb=t_comb))
}

test_e = function(e, null, verbose=FALSE) {
  # test for residuals
  tobs = tstat(e)
  pval = mean(null$H0 >= null$t_comb(tstat(e)))
  2*min(pval, 1-pval)
}

CHECK_P = function(P0, null_obj) {
  P1 = 4.2308977400789418155
  
  y = rv$val
  ts = rv$days
  X = cbind(t=ts, 
            xt1=cos(2*pi*ts/P0), 
            xt2=sin(2*pi*ts/P0))
  
  fit = lm(y ~ X)
  e0 = fit$residuals
  print(test_e(e0, null_obj, T) )
  print(shapiro.test(e0)$p.value)
  
  X = cbind(t=ts, 
            xt1=cos(2*pi*ts/P1), 
            xt2=sin(2*pi*ts/P1))
  
  fit = lm(y ~ X)
  e1 = fit$residuals # for P1
  
  
  print(sum(dnorm(e0, sd=sd(e0), log=T)))
  print(sum(dnorm(e1, sd=sd(e0), log=T)))
  
  par(mfrow=c(2, 1))
  hist(e0/sd(e0), col="blue", breaks=40, main=sprintf("P=%.4f", P0))
  hist(e1/sd(e1), col="blue", breaks=40, main=sprintf("P=%.4f", P1))
  
}


test_P = function(P0, use_basic=F) {
  # P0 = 4.2308977400789418155
  # P0=0.30821153249670213414
  # P0=0.80703187659959596534
  # P0=5.13063379048395784565
  # P0=0.12811818013837131258
  y = rv$val
  ts = rv$days
  X = cbind(out$X,
            xt1=cos(2*pi*ts/P0), 
            xt2=sin(2*pi*ts/P0))
  fit = NA
  if(all(X[,1] == 1)) {
    fit = lm(y ~ X + 0)
  } else {
    fit = lm(y ~ X)
  }
  e = fit$residuals
  
  # test_e(e, null_obj, T) 
  # shapiro.test(e)$p.value
 # hist(e/sd(e), breaks=40, col="blue")
  
  # rf = randomForest(y ~ ., data=df)
  
  # Test for residduals
  if(use_basic) {
    return(shapiro.test(e)$p.value) 
  } else {
    return(test_e(e, null_obj))
  }
  # lillie.test(e)
}

# hist(null$H0, col="blue")
t0 = proc.time()[3]
# out = load_basic("data/51Pegb/51Pegb.txt")
out = load_dumusque2012()
rv = out$rv
all_P = 10^seq(-1, 3, length.out=5000)
# Pegb = gen_rvObj("51Pegb", all_P)
Dumu = gen_rvObj("Dumusque2012", all_P)
Phat = Dumu$Phat
n  = length(rv$val)
# null_obj = gen_H0(n)

out = gen_rvObj("51Pegb", all_P)
out$find_planet_synth(1.3)

fit_peg = function(P0) {
  P0 = 0.80703187659959596534
   P0 = 4.2308977400789418155
  ts = rv$days
  y = rv$val
  Xp = cbind(b0=1, t=ts, cos=cos(2*pi*ts/P0), sin=sin(2*pi*ts/P0))
  fit = lm(y ~ Xp + 0)
  
  ei = fit$residuals
  b = tail(coef(fit), 2) # (b_cos, b_sin)
  plot(y, fit$fitted.values, pch=20)
  abline(0, 1, col="red")
  # sum e_i * t_i * (cos(2pi*ti/P0) b_sin - sin(2*pi*ti/P0) * b_cos)
  wi = ei * ts * (Xp[, 3] * b[2] - Xp[, 4] * b[1])
  sum(wi)
}


run_test = function() {
  pv = c()
  shapiro_pv = c()
  par(mfrow=c(1, 1))
  for(j in 1:length(all_P)) {
    P0 = all_P[j]
    # pv = c(pv, test_P(P0))
    pv = c(pv, 1)
    shapiro_pv = c(shapiro_pv, test_P(P0, use_basic = T))
    
    t1 = proc.time()[3]
    alpha_0 = 0.01
    
    if(t1- t0 > 5 | j==length(all_P)) {
      h = head(all_P, j)
      plot(h, pv, log="x", type="l", ylim=c(0, 0.1), col="blue",
           main="blue=combined test, red=shapiro test")
      lines(h, shapiro_pv, col="red")
      abline(h=0.01, col="blue", lty=3)
      print("Combined Test: Theta^ = ")

      print(sort(h[which(pv > 0.01)]))
      print("Shapiro Test: Theta^ = ")
      print(sort(h[which(shapiro_pv > 0.01)]))
      
      print(paste(" j  = ", j, " / ", length(all_P)))
      t0 = t1
    }
    
  }
  
  PlanetObj = Dumu
  pCand = all_P[which(shapiro_pv >= 0.05)]
  for(j in 1:length(pCand)) {
    print(paste("Testing ", pCand[j]))
    t0 = proc.time()[3]
    h1 = replicate(100, { PlanetObj$find_planet_synth(pCand[j])})
    hist(h1, breaks=30, col="blue", main=sprintf("P0=%.3f", pCand[j]))
    print(h1)
    f_obs = sum(abs(h1 - Phat) <= 1e-3) / 100
    
    t1 = proc.time()[3]
    print(paste("Took ", t1-t0, " secs., f_obs = ", f_obs))
    
    abline(v=Phat, col="red", lwd=2)
  }
}


rm(list=ls())
food = read.csv("Food_Inspections.csv")
head(food)
colnames(food)
summary(food)
