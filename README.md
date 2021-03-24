# Randomization Inference for Exoplanet Detection

This code implements the methods in: 
Toulis, P. and Bean, J. (2021). "Randomization Inference of Periodicity in Unequally Spaced Time Series with Application to Exoplanet Detection."

The main idea in this paper is to estimate the hidden periodicity, say of a candidate exoplanet, not by relying on standard statistical asymptotics, but by simulating the sampling distribution of the periodogram peak and then inverting a series of null hypotheses on the underlying period. 
With this approach, our method is robust to irregular observation times, especially in small samples.

See here https://www.ptoulis.com/s/astro_main.pdf for the latest version.
