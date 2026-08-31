context("Plot Function")

object_hrm <- HRM::hrm_test(value ~ group*dimension, subject = "subject", data = EEG)
object_hrm2 <- HRM::hrm_test(value ~ dimension, subject = "subject", data = EEG)


test_that("function plot", {
  expect_equivalent(class(object_hrm), "HRM")
  expect_equivalent(class(object_hrm2), "HRM")
  expect_true(inherits(plot(object_hrm), "ggplot"))
  expect_true(inherits(plot(object_hrm2), "ggplot"))
  expect_true(inherits(plot(object_hrm, xlab = "time", ylab = "mean", legend = FALSE,
                            legend.title = "", ggplot2::theme_bw() + ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position="none")),
                       "ggplot"))
})