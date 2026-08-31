# HRM 1.3.0

[![CRANstatus](https://www.r-pkg.org/badges/version/HRM)](https://cran.r-project.org/package=HRM)
[![](https://cranlogs.r-pkg.org/badges/HRM)](https://cran.r-project.org/package=HRM)
[![R-CMD-check](https://github.com/happma/HRM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/happma/HRM/actions/workflows/R-CMD-check.yaml)


R package for analysing high-dimensional repeated measures for factorial designs. A description of this package can be found in [1], theoretical derivations of the test statistics are in [2] and [3].



To install the current development version:

``` r
## install devtools package
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
# install package
devtools::install_github("happma/HRM", ref = "test", dependencies = TRUE)
library(HRM)
```

With this package it is possible to test for main and interaction effects of up to three whole- or subplot-factors. In total, a maximum of four factors can be used. There are two different S3 methods available. The first method requires a list of matrices in the wide table format. The second method requires a data.frame in the long table format.

``` r
## hrm_test with a list of matrices

# number patients per group
n = c(10,10)
# number of groups
a=2
# number of variables
d=40

# defining the list consisting of the samples from each group
mu_1 = mu_2 = rep(0,d)
# autoregressive covariance matrix
sigma_1 = diag(d)
for(k in 1:d) for(l in 1:d) sigma_1[k,l] = 1/(1-0.5^2)*0.5^(abs(k-l))
sigma_2 = 1.5*sigma_1
X = list(MASS::mvrnorm(n[1],mu_1, sigma_1), MASS::mvrnorm(n[2],mu_2, sigma_2))
X=lapply(X, as.matrix)

hrm_test(data=X, alpha=0.05)


## hrm.test with a data.frame using a 'formula' object

# using the EEG dataset
hrm_test(value ~ group*region*variable, subject = "subject", data = EEG)
```

## References

[1] Happ, M., Harrar, S. W., and Bathke, A. C. (2018). HRM: An R Package for Analysing High-dimensional Multi-factor Repeated Measures. The R Journal 10(1), 534--548. <a href="https://journal.r-project.org/articles/RJ-2018-032/index.html">https://journal.r-project.org/articles/RJ-2018-032/index.html</a>

[2] Happ, M., Harrar S. W. and Bathke, A. C. (2017). High-dimensional Repeated
  Measures. Journal of Statistical Theory and Practice. 11(3), 468-477. URL:
  <a href="https://doi.org/10.1080/15598608.2017.1307792">doi:10.1080/15598608.2017.1307792</a>.
  
[3] Happ, M., Harrar, S. W., & Bathke, A. C. (2016). Inference for low‐and high‐dimensional multigroup repeated measures designs with unequal covariance matrices. Biometrical Journal, 58(4), 810-830. <a href = "https://doi.org/10.1002/bimj.201500064">doi:10.1002/bimj.201500064</a>
