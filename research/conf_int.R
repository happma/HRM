## ---------------------------------------------------------------------------
## WITHDRAWN from the package in HRM 1.3.0. Kept here for future restoration.
##
## Two open issues, both measured on 2026-08-31 (R 4.4.2, mvtnorm 1.3.3):
##
## 1) BROKEN for designs with a whole-plot factor.
##    qmvnorm() below fails with "missing value where TRUE/FALSE needed"
##    because the correlation matrix R is SINGULAR whenever a group has
##    fewer subjects than repeated measures (rank <= n_i - 1 < d).
##    For value ~ group*region*variable: 160x160, rank 136.
##    See research/bug_confint_wholeplot.R for the full reproducer.
##    NOTE: this reportedly worked in the past, so it is likely an mvtnorm
##    regression -- the qmvnorm call here is unchanged since 2018-07-13.
##    Also verified against mvtnorm 1.3.6 (newer than the 1.3.3 in use):
##    STILL ERRORS. So this is not fixed by simply upgrading mvtnorm; the
##    singular case needs handling on our side.
##
## 2) SLOW -- this is why the example was commented out for CRAN.
##    Two independent costs, only one of which is fixable here:
##
##      a) The R[i, j] double loop below is O(m^4): m^2 cells, each doing
##         three full m x m matrix products. REPLACE IT WITH cov2cor():
##
##             R <- cov2cor(object$var)
##
##         Verified identical (max abs diff 3.3e-16) and ~29000x faster.
##             m=40 : loop 0.025 s -> cov2cor 0.000 s
##             m=160: loop 2.92  s -> cov2cor 0.000 s
##
##      b) qmvnorm() itself costs ~2.1 s at m=40 and is NOT fixable here --
##         it is Genz-Bretz Monte Carlo integration over m dimensions.
##         This dominated the old test suite (36 s -> under 1 s once the six
##         confint assertions were removed).
##
##    So cov2cor() alone will not make the example CRAN-fast. On restoration,
##    wrap the example in \donttest{} rather than commenting it out -- that is
##    the sanctioned mechanism for a correct-but-slow example.
##
## Removed test assertions are preserved in research/test_confint.R.
## ---------------------------------------------------------------------------

#' Function to calculate confidence intervals
#'
#' @description Function to calculate simultaneous, asymptotic (1-alpha) confidence intervals for an object of class 'HRM'.
#' @rdname confint.HRM
#' @param object an object from class 'HRM' returned from the function hrm_test
#' @param parm currently ignored; all possible confidence intervals are calculated
#' @param level confidence level (FWER) used for calculating the inverals
#' @param ... Further arguments passed to 'hrm_test' will be ignored
#' @return Returns a data.frame with mean and 1-alpha confidence interval for each factor combintation
#' @example inst/examples/confint.R
#' @keywords export
confint.HRM <- function(object, parm, level = 0.95, ...) {
  stopifnot(!is.null(object$formula))

  output <- summaryBy(object$formula, data = object$data, FUN = mean)
  ss <- summaryBy(object$formula, data = object$data, FUN = length)
  ss <- ss[[dim(ss)[2]]]

  ncol <- dim(output)[2]
  m <- dim(object$var)[1]
  CI_lower <- rep(0, m)
  CI_upper <- rep(0, m)
  c <- rep(0, m)
  alpha <- 1 - level

  R <- diag(m)
  # calculate correlation matrix
  #
  # FIXME (see header): this whole double loop is equivalent to
  #     R <- cov2cor(object$var)
  # which is ~29000x faster and identical to 3.3e-16. Kept as-is only so the
  # withdrawn code matches what shipped in 1.2.1.
  for(i in 1:m) {
    for(j in 1:m) {
      ci <- c*0
      cj <- c*0
      ci[i] <- 1
      cj[j] <- 1
      R[i, j] <- (t(ci)%*%object$var%*%cj)*1/sqrt(t(ci)%*%object$var%*%ci)*
        1/sqrt(t(cj)%*%object$var%*%cj)
    }
  }
  # FIXME (see header): fails here on singular R (whole-plot designs), and
  # costs ~2.1 s even at m=40. Not fixable without changing the approach.
  quantiles <- qmvnorm(level, corr = R, tail = "both")$quantile

  if(object$nonparametric) {

    if(object$factors[[1]] != "none") {
      grp <- object$data[, object$factors[[1]][1]]
      if(length(object$factors[[1]])  > 1) {
        for(i in 2:length(object$factors[[1]])) {
          grp <- paste(grp, object$data[, object$factors[[1]][i]], sep="")
        }
      }

      object$data$grouping <- as.factor(grp)

      setDT(object$data)
      object$data[,"prank" := 1/(dim(object$data)[1])*(pseudorank(object$data[[as.character(object$formula[[2]])]], object$data[, grouping]) - 1/2)]
    } else {
      setDT(object$data)
      object$data[,"prank" := 1/(dim(object$data)[1])*(rank(object$data[[as.character(object$formula[[2]])]], ties.method="average") - 1/2)]
    }
    new_formula <- as.formula(paste("prank ~", split(as.character(object$formula), "~")[[1]][3]))
    output <- as.data.frame(summaryBy(new_formula, data = object$data, FUN = mean))

    for(i in 1:m) {
      c <- c*0
      c[i] <- 1
      sdi <- sqrt(t(c)%*%object$var%*%c)
      CI_lower[i] <- output[i, ncol] - quantiles*sdi
      CI_upper[i] <- output[i, ncol] + quantiles*sdi
    }
  } else {
    output <- as.data.frame(output)
    for(i in 1:m) {
      c <- c*0
      c[i] <- 1
      sdi <- sqrt(t(c)%*%object$var%*%c)
      CI_lower[i] <- output[i, ncol] - quantiles*sdi
      CI_upper[i] <- output[i, ncol] + quantiles*sdi
    }
  }

  output$CI_lower <- CI_lower
  output$CI_upper <- CI_upper

  return(output)
}
