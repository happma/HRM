## ---------------------------------------------------------------------------
## confint() regression assertions, removed from tests/testthat/ in HRM 1.3.0
## when confint.HRM was withdrawn (see bug_confint_wholeplot.R).
##
## All six passed against 0-whole-plot designs. Restore them alongside
## conf_int.R if/when the singular-correlation-matrix issue is resolved --
## and add whole-plot cases, which were never covered.
## ---------------------------------------------------------------------------

context("confint (withdrawn)")

## --- from test_0w_1s.R -----------------------------------------------------
result_CI <- c(1, 3.9380174, 3.802278, 4.073756)
expect_equivalent(as.numeric(confint(HRM::hrm_test(
  value ~ dimension, subject = "subject", data = EEG))[1, ]), result_CI, tol = 1e-4)

result_CI <- c(1, 0.91993262, 0.9196157, 0.9202495)
expect_equivalent(as.numeric(confint(HRM::hrm_test(
  value ~ dimension, subject = "subject", data = EEG,
  nonparametric = TRUE, np.correction = FALSE))[1, ]), result_CI, tol = 1e-4)

## --- from test_0w_2s.R -----------------------------------------------------
result_CI <- c(1, 1, 3.9380174, 3.802331, 4.073704)
expect_equivalent(as.numeric(confint(HRM::hrm_test(
  value ~ region * variable, subject = "subject", data = EEG))[1, ]), result_CI, tol = 1e-4)

result_CI <- c(1, 1, 0.91993262, 0.9196158, 0.9202494)
expect_equivalent(as.numeric(confint(HRM::hrm_test(
  value ~ region * variable, subject = "subject", data = EEG,
  nonparametric = TRUE))[1, ]), result_CI, tol = 1e-4)

## --- from test_0w_3s.R (needs the factor1/2/3 columns built there) ---------
result_CI <- c(1, 1, 1, 3.9380174, 3.802278, 4.073756)
expect_equivalent(as.numeric(confint(HRM::hrm_test(
  value ~ factor3 * factor2 * factor1, subject = "subject", data = dat))[1, ]), result_CI, tol = 1e-4)

result_CI <- c(1, 1, 1, 0.91993262, 0.9196157, 0.9202495)
expect_equivalent(as.numeric(confint(HRM::hrm_test(
  value ~ factor3 * factor2 * factor1, subject = "subject", data = dat,
  nonparametric = TRUE, np.correction = FALSE))[1, ]), result_CI, tol = 1e-4)
