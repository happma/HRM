# HRM

[![CRANstatus](https://www.r-pkg.org/badges/version/HRM)](https://cran.r-project.org/package=HRM)
<a href="https://www.rpackages.io/package/HRM"><img src="https://www.rpackages.io/badge/HRM.svg" /></a>
[![](https://cranlogs.r-pkg.org/badges/HRM)](https://cran.r-project.org/package=HRM)

R package for analysing high-dimensional repeated measures for factorial designs. 

To install the current development version:

``` r
## install devtools package
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
# install package
devtools::install_github("happma/HRM")
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
X = list(mvrnorm(n[1],mu_1, sigma_1), mvrnorm(n[2],mu_2, sigma_2))
X=lapply(X, as.matrix)

hrm_test(data=X, alpha=0.05)


## hrm.test with a data.frame using a 'formula' object

# using the EEG dataset
hrm_test(value ~ group*region*variable, subject = "subject", data = EEG)
```



Additionally, the package can be used with a GUI.
``` r
hrm_GUI()
```

For more information, see

Happ, M., Harrar S. W. and Bathke, A. C. (2018). HRM: An R Package for Analysing High-dimensional Multi-factor Repeated Measures. The R Journal. URL:
  <a href="https://journal.r-project.org/archive/2018/RJ-2018-032/index.html">https://journal.r-project.org/archive/2018/RJ-2018-032/index.html</a>.

Happ, M., Harrar S. W. and Bathke, A. C. (2017). High-dimensional Repeated
  Measures. Journal of Statistical Theory and Practice. 11(3), 468-477. URL:
  <a href="http://www.tandfonline.com/doi/full/10.1080/15598608.2017.1307792">http://www.tandfonline.com/doi/full/10.1080/15598608.2017.1307792</a>.
