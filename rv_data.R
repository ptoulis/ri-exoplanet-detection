#' Load all RV datasets in one object.
#'
print("> LOADING rv_data.R")
options(digits=10)
source("rv_lib.R")

DATASETS = list(list(name="51Pegb", nper=25000), 
                list(name="GJ436b", nper=30000), 
                list(name="HD85512b", nper=25000), 
                list(name="Dumusque2012", nper=10000),
                list(name="ESPRESSO", nper=10000))

load_Dataset_Index = function(IDX) {
  stopifnot(IDX %in% seq(1, length(DATASETS)))
  all_P = 10^seq(-1, 3, length.out=DATASETS[[IDX]]$nper) # periods
  rvData = load_Dataset(DATASETS[[IDX]]$name, all_P)  # RV dataset.
  return(list(RV=rvData, Periods=all_P))
}
#' RV dataset is 
#' (days, val, se, X1, X2, .....)

get_rv_blankSep = function(out) {
  # removes Blank separator from RV data file.
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
  ## Add sampling function
  return(rv)
}

load_HD85512b = function() {
  out = read.csv("datasets/HD85512b.txt", header=F, stringsAsFactors = FALSE)
  out = t(apply(out,  1, function(row) {
    s = as.numeric(strsplit(row, split=" ")[[1]])
    s[!is.na(s)]
  }))
  
  rv = data.frame(days=out[,1], val=out[,2], se=out[,3])
  rv$days = rv$days -  min(rv$days)
  Xrv = out[, c(4, 6, 8)]
  
  rv = cbind(rv, FWH=Xrv[,1], Spa=Xrv[,2], logRhK=Xrv[,3])
  return(rv)
}


load_dumusque2012 = function() { 
  # print("Dumusque data.")
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
  BinX = cbind(ts2=rv$days^2)
  
  # 2. Magnetic cycle
  MagX = matrix(astro$Low.pass.filter.log.R.hk., ncol=1)
  colnames(MagX) = "logRhk"
  
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
  rv = cbind(rv, Xrv)
  return(rv)
}

load_ESPRESSO = function() {
  rvData = read.csv("datasets/Proxima_RV.csv", header=T)
  rvData = subset(rvData, days > 8524) # ESPRESSO?
  nrow(rvData)
  
  # rvData$days = rvData$days - min(rvData$days)
  fwData = read.csv("datasets/Proxima_FWHM.csv", header=T)
  fwData = subset(fwData, days > 8524) # ESPRESSO?
  nrow(fwData)
  fwData$Spec=NULL
  # rvData$days = rvData$days - min(rvData$days)
  # fwData$days = fwData$days - min(fwData$days)
  stopifnot(fwData$days == rvData$days)
  
  # rvData = cbind(rvData, fw=fwData$val, fw_se=fwData$se)
  return(list(rv=rvData, FWHM=fwData))
}

load_Dataset = function(name, all_P=c()) {
  print(paste(">> Loading dataset: ", name))
  if(name=="51Pegb") {
    return(load_basic("datasets/51Pegb.txt"))
  } else if(name=="GJ436b") {
    return(load_basic("datasets/GJ436b.txt"))
  } else if(name=="HD85512b") {
    out = load_HD85512b()
    X = as.matrix(out[,seq(4, ncol(out))])
    ts = out$days
    fit1 = lm(out$val ~ ts + X) # fits linear trend.
    rvData = data.frame(val=fit1$residuals, days=ts, se=out$se)
    return(rvData)
    
  } else if(name=="Dumusque2012") {
    
    out = load_dumusque2012()
    X = as.matrix(out[,seq(4, ncol(out))])
    ts = out$days
    fit1 = lm(out$val ~ ts + X)
    rvData = data.frame(val=fit1$residuals, days=ts, se=out$se)
    return(rvData)
    #
  } else if(name=="ESPRESSO") {
    ll = load_ESPRESSO()
    stopifnot(length(all_P) > 0) # requires some pre-processing
    print("> Loading FWHM data to remove from signal..")
    rvData = ll$rv
    fwData = ll$FWHM
    
    fw0 = lombe_scragle_fast(fwData, all_P, v=F)
    ts = rvData$days
    
    #X2 = cbind(cos(2*pi*ts/fw0$Phat), sin(2*pi*ts/fw0$Phat),
    #         cos(4*pi*ts/fw0$Phat), sin(4*pi*ts/fw0$Phat))
    X = cbind(cos(2*pi*ts/fw0$Phat), sin(2*pi*ts/fw0$Phat))
    
    y = rvData$val
    # Remove other stellar signal.
    fit1 = lm(y ~ ts + X + fwData$val)
    y = fit1$residuals
    rvData$val = y
    
    return(rvData) # FW data have been removed.
  }
  stop(sprintf("Dataset %s does not exist.", name))
}
