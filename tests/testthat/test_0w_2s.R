context("0 Wholeplot, 2 Subplot")
dat <- EEG
dat$value2 <- exp(dat$value)

true_result <- c(5.542208, 881.211121, 2.157672, 225.554716)
result <- as.numeric(HRM::hrm_test(value ~ region*variable, subject = "subject", data = dat)$result[1, 2:5])

test_that("function hrm_test", {
  expect_equal(result, true_result, tol = 1e-4)
})

true_result <- c(4.571516, 726.871040, 2.287492, 46.875637)
result <- as.numeric(HRM::hrm_test(value ~ region*variable, subject = "subject", data = dat, nonparametric = TRUE)$result[1, 2:5])
result2 <- as.numeric(HRM::hrm_test(value2 ~ region*variable, subject = "subject", data = dat, nonparametric = TRUE)$result[1, 2:5])

test_that("function hrm_test, nonparametric", {
  expect_equal(result, true_result, tol = 1e-4)
  expect_equal(result2, true_result, tol = 1e-4)
})