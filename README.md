# Randomization Inference for Exoplanet Detection

This code implements the methods in: 
Toulis, P. and Bean, J. (2021). "Randomization Inference of Periodicity in Unequally Spaced Time Series with Application to Exoplanet Detection."

The main idea in this paper is to estimate the hidden periodicity, say of a candidate exoplanet, not by relying on standard statistical asymptotics, but by simulating the sampling distribution of the periodogram peak and then inverting a series of null hypotheses on the underlying period. 
With this approach, our method is robust to irregular observation times, especially in small samples.

See here https://www.ptoulis.com/s/astro_main.pdf for the latest version.

## Example 1: Create radial velocity data frame.

The appropriate format for radial velocity is as a `data.frame(days, val, se)` where `days` are the JB day measurements, `val` are the radial velocity measurements and `se` are the standard errors.

    set.seed(1)
    n = 100
    var_days = 0.05
    Tn = seq(1, n) + var_days*runif(n, min=-1, max=1)
    theta_true = sqrt(2) 

    Yn = 1.5 * cos(2*pi*Tn/theta_true) + rnorm(n)    # synthetic RV measurements (true period is ~=1.41 days)
    rv = data.frame(val=Yn, days=Tn, se=rep(1, n))   # RV data frame.

## Example 2: Generate periodogram for the RV data.

We use the generalized version of the Lombe-Scragle periodogram. See [(Zechmeister and M. Kürster, 2009)](https://www.aanda.org/articles/aa/pdf/2009/11/aa11296-08.pdf)

    source("rv_lib.R")
    all_P = 10^seq(-1, log(50, base=10), length.out=10000)
    ls0 = lombe_scragle_fast(rv, all_P, v=T)
    
    ls0$Phat                                        # Periodogram peak should be at 3.14 days.

## Example 3. Hypothesis testing

Let's test whether the true period could be equal to sqrt(2), in notation H0: θ* = sqrt(2). Since the truth is equal to sqrt(2) we should **not** reject the null.

    test = Test_Period_Theta0(rv, all_P, theta0=theta_true, null_samples=100)
    test$pval                                       # pval = 0.38. We don't reject.
    
## Example 4. Build 99% confidence set.

Here, we will construct a 99% confidence set for the unknown periodicity. The function `Build_ConfidenceSet` implements Procedure 1 in the [paper]( https://www.ptoulis.com/s/astro_main.pdf). 

    ci = Build_ConfidenceSet(rv, all_P, null_samples=100, time_budget_mins = 2)

This confidence set will contain the values `{0.586, 0.774, 1.413, 3.417}`. We see that the true value is included (as expected from the above test).

The construction below is designed to be approximate and fast. In a full application we need to increase `null_samples` and the `time_budget_mins` so that inference relies on more samples. If `m` is the number of periods to be tested for inclusion in the confidence set, the code calculates `m` by solving the following equation:

    m * null_samples * C = time_budget_mins

where `C` is wall-clock time / periodogram calculation --- this is estimated through a quick simulation so that we can adjust to the current computing environment.

## Example 5. Real exoplanet detection: 51 Pegasi b

    source("rv_data.R")
    all = load_Dataset_Index(1)                                 # 51 Peg b. For indexing details see `rv_data.R`.
    rvPeg = all$RV                                              # load RV dataset.
    all_P = all$Periods                                         # universe of periods (denoted Θ in paper)
    ls0 = lombe_scragle_fast(rvPeg, all_P, v=T)                 # LS periodogram of 51Pegasi B
    ls0$Phat                                                    # Peak should be at 4.23-days.
    
    ci = Build_ConfidenceSet(rv, all_P, time_budget_mins = 5)   # Build (fast) 99% confidence set for unknown period.

The confidence set will be a singleton `{4.23}` indicating that the underlying periodicity can be sharply identified.

## Example 6. Real exoplanet detection: candidate exoplanet around α Centauri B

See [(Dumusque et al, 2012)](https://www.nature.com/articles/nature11572) for details.

    all = load_Dataset_Index(1)                                 # α Centauri B data.
    rvD = all$RV                                              
    all_P = all$Periods                                        
    ls0 = lombe_scragle_fast(rvD, all_P, v=T)                 
    ls0$Phat                                                    # Peak should be at 3.24 days.
    
    ci = Build_ConfidenceSet(rvD, all_P, time_budget_mins = 5)  # Build (fast) 99% confidence set for unknown period.

The confidence set will be `{0.762,  3.236,  8.118, 61.133}` indicating that the unnderlying period cannot be identified with this dataset.

## Example 7. Real exoplanet detection: candidate exoplanet around Proxima Centauri

See [(Anglada-Escude, 2016)](https://www.nature.com/articles/nature19106) for details.
For the analysis here we use the more recent and precise ESPRESSO measurements found here (https://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/A+A/639/A77)
See also [(Suarez Mascareno et al, 2020)](https://arxiv.org/abs/2005.12114) for details on the ESPRESSO data.

    all = load_Dataset_Index(5)                                 # Proxima Centauri data.
    rvE = all$RV 
    all_P = all$Periods               
    ls0 = lombe_scragle_fast(rvE, all_P, v=T)                 
    ls0$Phat                                                    # Peak should be at 11.17 days.
    
    ci = Build_ConfidenceSet(rvE, all_P, time_budget_mins = 10)  # Requires more samples than before.

The confidence set will be `{0.916, 11.17}` indicating that the detection appears to be robust. However, there is a nuisance signal at 0.916-days that cannot be rejected at the 1% level. With few more observations (additional 10-15) this nuisance signal could be eliminated. See Section 6 in the paper.

