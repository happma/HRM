context("confint")

dat <- EEG
dat$factor1 <- as.factor(rep(1:2, 3200))
dat$factor2 <- as.factor(rep(c(1,1,2,2), 1600))
dat$factor3 <- as.factor(rep(kronecker(1:10, rep(1,4)), 160))


# ------------------------------------------------------------------------
# No whole-plot factor: qmvnorm() succeeds on the full correlation matrix,
# so these take the global-FWER path. Values frozen from HRM 1.2.1.
# ------------------------------------------------------------------------

test_that("confint, 0 wholeplot factors, parametric", {
  expect_equivalent(as.numeric(confint(HRM::hrm_test(
    value ~ dimension, subject = "subject", data = EEG))[1, ]),
    c(1, 3.9380174, 3.802278, 4.073756), tol = 1e-4)
  expect_equivalent(as.numeric(confint(HRM::hrm_test(
    value ~ region*variable, subject = "subject", data = EEG))[1, ]),
    c(1, 1, 3.9380174, 3.802331, 4.073704), tol = 1e-4)
  expect_equivalent(as.numeric(confint(HRM::hrm_test(
    value ~ factor3*factor2*factor1, subject = "subject", data = dat))[1, ]),
    c(1, 1, 1, 3.9380174, 3.802278, 4.073756), tol = 1e-4)
})

test_that("confint, 0 wholeplot factors, nonparametric", {
  expect_equivalent(as.numeric(confint(HRM::hrm_test(
    value ~ dimension, subject = "subject", data = EEG,
    nonparametric = TRUE, np.correction = FALSE))[1, ]),
    c(1, 0.91993262, 0.9196157, 0.9202495), tol = 1e-4)
  expect_equivalent(as.numeric(confint(HRM::hrm_test(
    value ~ region*variable, subject = "subject", data = EEG,
    nonparametric = TRUE))[1, ]),
    c(1, 1, 0.91993262, 0.9196158, 0.9202494), tol = 1e-4)
  expect_equivalent(as.numeric(confint(HRM::hrm_test(
    value ~ factor3*factor2*factor1, subject = "subject", data = dat,
    nonparametric = TRUE, np.correction = FALSE))[1, ]),
    c(1, 1, 1, 0.91993262, 0.9196157, 0.9202495), tol = 1e-4)
})

test_that("confint reports global FWER when qmvnorm succeeds", {
  ci <- confint(HRM::hrm_test(value ~ region*variable, subject = "subject", data = EEG))
  expect_match(attr(ci, "status"), "over all factor level combinations")
  expect_false(any(is.na(ci$CI_lower)))
  expect_false(any(is.na(ci$CI_upper)))
})


# ------------------------------------------------------------------------
# With a whole-plot factor the correlation matrix is large enough that
# qmvnorm() fails upstream. confint() must not error: it falls back to
# FWER within each group and says so.
# ------------------------------------------------------------------------

test_that("confint does not error on wholeplot designs", {
  expect_error(confint(HRM::hrm_test(
    value ~ group*dimension, subject = "subject", data = EEG)), NA)
  expect_error(confint(HRM::hrm_test(
    value ~ group*region*variable, subject = "subject", data = EEG)), NA)
  expect_error(confint(HRM::hrm_test(
    value ~ group*dimension, subject = "subject", data = EEG,
    nonparametric = TRUE)), NA)
})

test_that("confint falls back to per-group FWER on wholeplot designs", {
  z  <- HRM::hrm_test(value ~ group*region*variable, subject = "subject", data = EEG)
  ci <- confint(z)

  expect_match(attr(ci, "status"), "within each group")
  # one row per factor level combination, and every interval is finite
  expect_equal(nrow(ci), dim(z$var)[1])
  expect_false(any(is.na(ci$CI_lower)))
  expect_false(any(is.na(ci$CI_upper)))
  # intervals must bracket the estimate
  est <- ci[, ncol(ci) - 2]
  expect_true(all(ci$CI_lower <= est))
  expect_true(all(est <= ci$CI_upper))
})

test_that("confint per-group quantiles differ between groups", {
  ci <- confint(HRM::hrm_test(value ~ group*region*variable,
                              subject = "subject", data = EEG))
  est   <- ci[, ncol(ci) - 2]
  halfw <- ci$CI_upper - est
  # 4 groups x 40 combinations; each group gets its own quantile, so the
  # half-widths must not all come from a single global value
  grp <- rep(1:4, each = nrow(ci) / 4)
  expect_gt(length(unique(round(tapply(halfw, grp, max), 8))), 1)
})

test_that("confint level argument widens the intervals", {
  z   <- HRM::hrm_test(value ~ group*dimension, subject = "subject", data = EEG)
  c95 <- confint(z, level = 0.95)
  c99 <- confint(z, level = 0.99)
  est <- c95[, ncol(c95) - 2]
  expect_true(all(c99$CI_lower <= c95$CI_lower))
  expect_true(all(c95$CI_upper <= c99$CI_upper))
})
