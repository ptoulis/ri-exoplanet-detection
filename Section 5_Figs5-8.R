#' Analysis of 4 exoplanets in Section 5 in Toulis, P. and Bean, J. (2021). "Randomization Inference of Periodicity in Unequally Spaced Time Series with Application to Exoplanet Detection."
#' https://www.ptoulis.com/s/astro_main.pdf
#' 
#' Panos Toulis, ptoulis@chicagobooth.edu, 03/2021
#' 
#' Here, we load results pre-computed on a cluster.
#' If you want to replicate a quick (and approximate) version of the figures 
#' look at the end of the script file (Line 66).
rm(list=ls())
source("rv_lib.R")
source("rv_data.R")
require(plyr)
require(dplyr)

#' Dataset to analyze
#' 1 = 51 PEG b  --- Figure 5
#' 2 = GJ 436 b  --- Figure 6
#' 3 = HD 85512 b
#' 4 = α Centauri B  --- Figure 7
#' 5 = Proxima Centauri (ESPRESSO data)  --- Figure 8
I = 1

load("ALL_RESULTS_R=100000.rdata")  # Pre-computed data. Has 6 columns
ALL_RESULTS = as.data.frame(ALL_RESULTS)
A = ddply(ALL_RESULTS, .(dataset), summarize, theta=theta0[which(!reject)])

dname = DATASETS[[I]]$name

# Load data (rvData, periods)
D = load_Dataset_Index(I)
rvData = D$RV
all_P = D$Periods

# LS periodogram
ls0 = lombe_scragle_fast(rvData, all_P, v=T)
ls0$Phat
print(paste(">> Dataset", dname))
print(">> Significant theta0 at the 0.1% level:")
thetas = subset(A, dataset==I)$theta

# Grab results for particular dataset.
C = subset(ALL_RESULTS, dataset==I)
# 
Pf = sapply(C$theta0, function(theta) {
  id = which.min(abs(all_P-theta))
  ls0$Pf[id]
})

idx = head(rev(order(Pf)), 12)
C2 = C[idx,]
rownames(C2) = seq(1, nrow(C2))
idx2 = order(C2$theta0)
C3 = C2[idx2, ]

pval = C3[,5]
C4 = data.frame(theta0=round(C3$theta0, 4), pval=round(pval, 4),
                incl95=ifelse(pval >= 0.05, "yes", "no"),
                incl99=ifelse(pval >= 0.01, "yes", "no")
)
rownames(C4) = NULL
C4
require(xtable)
print(xtable(C4, digits = 4), include.rownames=FALSE)


#' Quick (and approximate) replication.
#' 

#' Simple version?
Theta_Hat = Build_ConfidenceSet(rvData, all_P, time_budget_mins = 5)

