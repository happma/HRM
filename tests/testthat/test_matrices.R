context("hrm_test with matrices")

# number patients per group
n <- c(10,10)
# number of groups
a<-2
# number of variables
d<-40

# defining the list consisting of the samples from each group
mu_1 <- rep(0,d)
mu_2 <- rep(0,d)
# autoregressive covariance matrix
sigma_1 <- diag(d)
for(k in 1:d) for(l in 1:d) sigma_1[k,l] <- 1/(1-0.5^2)*0.5^(abs(k-l))
sigma_2 <- 1.5*sigma_1
set.seed(12345)
X <- list(mvrnorm(n[1],mu_1, sigma_1), mvrnorm(n[2],mu_2, sigma_2))
X<-lapply(X, as.matrix)

true_result <- c(1, 19.8267291, 4.3561230, 0.3923543)
t_matrices <- hrm_test(data=X, alpha=0.05)
result <- as.numeric(t_matrices$result[1,2:5])

X[[1]] <- data.frame(value = vec(t(X[[1]])), group = 1, time = rep(1:40,n[1]))
X[[2]] <- data.frame(value = vec(t(X[[2]])), group = 2, time = rep(1:40,n[2]))
X <- rbind(X[[1]], X[[2]])
X$group <- as.factor(X$group)
X$subject <- gl(sum(n), d)

result2 <- as.numeric(hrm_test(value ~ group*time, data = X, subject = "subject")$result[1,2:5])

test_that("function hrm_test with matrices", {
  expect_equal(result, true_result, tol = 1e-4)
  expect_equal(result2, true_result, tol = 1e-4)
  expect_output(summary(t_matrices))
  expect_output(print(t_matrices))
  expect_equal(class(plot(t_matrices)), c("gg", "ggplot"))
})
