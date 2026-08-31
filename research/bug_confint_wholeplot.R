## ---------------------------------------------------------------------------
## Reproducer: confint() fails for any design with a whole-plot factor
##
## Status : open
## Found  : 2026-08-31, during the CRAN resubmission audit (HRM 1.3.0)
## Affects: parametric AND nonparametric; every design with >= 1 whole-plot factor
## ---------------------------------------------------------------------------

library(HRM)

## --- 1. What works: no whole-plot factor -----------------------------------
z_ok <- hrm_test(value ~ region * variable, subject = "subject", data = EEG)
head(confint(z_ok), 3)          # fine


## --- 2. What fails: add a whole-plot factor --------------------------------
z_bad <- hrm_test(value ~ group * dimension, subject = "subject", data = EEG)
confint(z_bad)
#> Error: missing value where TRUE/FALSE needed

## Same failure without the nonparametric branch, so it is NOT rank-specific:
confint(hrm_test(value ~ group * dimension, subject = "subject",
                 data = EEG, nonparametric = TRUE))
#> Error: missing value where TRUE/FALSE needed

## This is the example printed in README.md, verbatim -- it is broken too:
z_readme <- hrm_test(value ~ group * region * variable,
                     subject = "subject", data = EEG)
confint(z_readme, level = 0.99)
#> Error: missing value where TRUE/FALSE needed


## --- 3. Upstream bug, plus an observation about the input ------------------
## The failure is an UPSTREAM BUG in mvtnorm::qmvnorm(), not a legitimate
## rejection of the input: this call worked in the past and is unchanged
## since 2018-07-13. Reproduced under mvtnorm 1.3.3 and 1.3.6.
##
## Separately, R happens to be singular for these designs. That is recorded
## below because it is useful context, but it is NOT the established cause.

corr_of <- function(z) {
  m <- dim(z$var)[1]; R <- diag(m); c <- rep(0, m)
  for (i in 1:m) for (j in 1:m) {
    ci <- c * 0; cj <- c * 0; ci[i] <- 1; cj[j] <- 1
    R[i, j] <- (t(ci) %*% z$var %*% cj) *
      1 / sqrt(t(ci) %*% z$var %*% ci) * 1 / sqrt(t(cj) %*% z$var %*% cj)
  }
  R
}

for (nm in c("z_ok", "z_bad")) {
  R  <- corr_of(get(nm))
  ev <- eigen(R, only.values = TRUE)$values
  cat(sprintf("%-6s dim=%3d rank=%3d min.eigenvalue=%.3e\n",
              nm, nrow(R), qr(R)$rank, min(ev)))
}
#> z_ok   dim= 40 rank= 40 min.eigenvalue=1.467e-03   <- positive definite
#> z_bad  dim=160 rank=136 min.eigenvalue=-1.690e-15  <- singular

## Confirm qmvnorm is the failing call, not anything upstream:
R_bad <- corr_of(z_bad)
stopifnot(!anyNA(R_bad), isTRUE(all.equal(R_bad, t(R_bad))))  # R itself is clean
mvtnorm::qmvnorm(0.95, corr = R_bad, tail = "both")
#> Error: missing value where TRUE/FALSE needed

## And it is the singularity, not the size -- a 160x160 identity is fine:
mvtnorm::qmvnorm(0.95, corr = diag(160), tail = "both")$quantile
#> 3.5979


## --- 4. Why the test suite never caught it ---------------------------------
## confint() IS asserted 6 times, but only in test_0w_1s.R, test_0w_2s.R and
## test_0w_3s.R -- all "0w" designs, i.e. zero whole-plot factors. The whole-
## plot code path in conf_int.R (the `object$factors[[1]] != "none"` branch)
## has no coverage at all.
##
## inst/examples/confint.R has the call commented out inside "## Not run:",
## so R CMD check never executes confint() either.
