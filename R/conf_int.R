####################################################################################################################################
### Filename:    conf_int.R
### Description: Simultaneous confidence intervals for an object of class 'HRM'.
###
###              mvtnorm::qmvnorm() fails with "missing value where TRUE/FALSE
###              needed" once the correlation matrix gets large. This is an
###              upstream bug, not a rejection of the input: it depends on the
###              dimension together with the correlation structure, not on rank
###              (a full-rank 60x60 fails, a singular 40x40 does not), and it is
###              deterministic. Rather than abort, we fall back through weaker
###              forms of multiplicity control and report in attr(, "status")
###              exactly which one was used. See research/bug_confint_wholeplot.R.
####################################################################################################################################


#' Function to calculate confidence intervals
#'
#' @description Function to calculate simultaneous, asymptotic (1-alpha) confidence intervals for an object of class 'HRM'.
#' @rdname confint.HRM
#' @param object an object from class 'HRM' returned from the function hrm_test
#' @param parm currently ignored; all possible confidence intervals are calculated
#' @param level confidence level (FWER) used for calculating the inverals
#' @param ... Further arguments passed to 'hrm_test' will be ignored
#' @return Returns a data.frame with mean and 1-alpha confidence interval for each factor level combination. The attribute \code{"status"} states which form of multiplicity control was used: family-wise error rate over all factor level combinations, family-wise error rate within each whole-plot group, a Sidak correction within each group, or none.
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

  # correlation matrix. Up to 1.2.1 this was a double loop doing three m x m
  # matrix products per cell, i.e. O(m^4); cov2cor() is identical to 3.3e-16
  # and ~29000x faster (2.9 s -> 0 s at m = 160).
  R <- stats::cov2cor(object$var)
  # qmvnorm() can fail with "missing value where TRUE/FALSE needed". This is
  # an upstream bug in mvtnorm, not a legitimate rejection of the input: it is
  # driven by the dimension of the correlation matrix together with its
  # correlation structure, not by rank (a full-rank 60x60 fails, a singular
  # 40x40 does not). The failure is deterministic.
  #
  # We therefore try progressively weaker forms of multiplicity control and
  # report in $status exactly which one was used:
  #   1. FWER over all factor level combinations   (qmvnorm on all of R)
  #   2. FWER within each whole-plot group          (qmvnorm per block of R)
  #   3. Sidak correction within each group         (closed form, no qmvnorm)
  #   4. nothing                                    (CI columns are NA)

  # object$var is block diagonal with one block per whole-plot group, the
  # blocks in the same order as the rows of 'output'.
  wholeplot <- object$factors[[1]]
  a <- if(identical(as.character(wholeplot)[1], "none")) {
    1L
  } else {
    as.integer(prod(vapply(wholeplot,
      function(f) nlevels(as.factor(object$data[[f]])), integer(1L))))
  }
  if(is.na(a) || a < 1L || m %% a != 0L) {
    a <- 1L
  }
  bs     <- m / a
  blocks <- split(seq_len(m), rep(seq_len(a), each = bs))

  status    <- NULL
  quantiles <- rep(NA_real_, m)

  # 1. FWER over everything
  q_global <- tryCatch(qmvnorm(level, corr = R, tail = "both")$quantile,
                       error = function(e) NULL)

  if(!is.null(q_global) && is.finite(q_global)) {
    quantiles[] <- q_global
    status <- "Simultaneous confidence intervals; FWER controlled over all factor level combinations."
  } else if(a > 1L) {
    # 2. FWER within each whole-plot group
    q_block <- lapply(blocks, function(idx)
      tryCatch(qmvnorm(level, corr = R[idx, idx, drop = FALSE], tail = "both")$quantile,
               error = function(e) NULL))
    ok <- vapply(q_block, function(q) !is.null(q) && is.finite(q), logical(1L))
    if(all(ok)) {
      for(i in seq_len(a)) {
        quantiles[blocks[[i]]] <- q_block[[i]]
      }
      status <- "Simultaneous confidence intervals; FWER controlled within each group. FWER over all factor level combinations was not computable."
    }
  }

  if(is.null(status)) {
    # 3. Sidak within each group
    q_sidak <- stats::qnorm(1 - (1 - (1 - alpha)^(1/bs))/2)
    if(is.finite(q_sidak)) {
      quantiles[] <- q_sidak
      status <- paste0("Sidak-corrected confidence intervals within each group (",
                       bs, " comparisons per group). FWER was not computable.")
    }
  }

  if(is.null(status)) {
    # 4. nothing worked
    quantiles[] <- NA_real_
    status <- "Calculation of confidence interval was not possible."
  }

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

  } else {
    output <- as.data.frame(output)
  }

  for(i in 1:m) {
    c <- c*0
    c[i] <- 1
    sdi <- sqrt(t(c)%*%object$var%*%c)
    CI_lower[i] <- output[i, ncol] - quantiles[i]*sdi
    CI_upper[i] <- output[i, ncol] + quantiles[i]*sdi
  }

  output$CI_lower <- CI_lower
  output$CI_upper <- CI_upper

  attr(output, "status") <- status

  return(output)
}
