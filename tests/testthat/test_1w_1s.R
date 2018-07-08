context("1 Wholeplot, 1 Subplot")
true_result <- c(2.993584, 130.107638, 2.676095, 1.881486)
result <- as.numeric(HRM::hrm_test(value ~ group*dimension, subject = "subject", data = EEG)$result[1, 2:5])

test_that("function hrm_test", {
  expect_equal(result, true_result, tol = 1e-4)
})

true_result <- c(3.077845, 137.450499, 2.648275, 5.318450)
dat <- EEG
EEG$value2 <- exp(EEG$value)
result <- as.numeric(HRM::hrm_test(value ~ group*dimension, subject = "subject", data = EEG, nonparametric = TRUE)$result[1, 2:5])
result2 <- as.numeric(HRM::hrm_test(value2 ~ group*dimension, subject = "subject", data = dat, nonparametric = TRUE)$result[1, 2:5])

test_that("function hrm_test, nonparametric", {
  expect_equal(result, true_result, tol = 1e-4)
  expect_equal(result2, true_result, tol = 1e-4)
})