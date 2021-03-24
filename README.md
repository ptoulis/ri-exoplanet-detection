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

Here, we will construct a 99% confidence set for the unknown periodicity. The function `Build_ConfidenceSet` implements Procedure 1 in the [paper]( https://www.ptoulis.com/s/astro_main.pdf). The construction below is designed to be approximate and fast. In a full application we need to increase `num_samples` and the `time_budget_mins` so that inference relies on more samples. See paper for details.

    ci = Build_ConfidenceSet(rv, all_P, time_budget_mins = 2)

This confidence set will contain the values `{0.586, 0.774, 1.413, 3.417}`. We see that the true value is included (as expected from the above test).

## Example 5. Real exoplanet detection problem: 51 Pegasi b

    source("rv_data.R")
    rvPeg = load_Dataset("51Pegb") # load RV dataset.
    all_P = 10^seq(-1, 2.5, length.out=25000)                   # universe of periods (denoted Θ in paper)
    ls0 = lombe_scragle_fast(rvPeg, all_P, v=T)                 # LS periodogram of 51Pegasi B
    ls0$Phat                                                    # Peak should be at 4.23-days.
    
    ci = Build_ConfidenceSet(rv, all_P, time_budget_mins = 5)   # Build (fast) 99% confidence set for unknown period.

The confidence set will be a singleton `{4.23}` indicating that the underlying periodicity can be sharply identified.

