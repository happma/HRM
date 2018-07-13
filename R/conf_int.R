#' Function to calculate confidence intervals
#' 
#' @description Function to calculate confidence intervals for an object of class 'HRM'.
#' @rdname confint.HRM
#' @param object an object from class 'HRM' returned from the function hrm_test
#' @param parm currently ignored; all possible confidence intervals are calculated
#' @param level confidence level used for calculating the inverals
#' @param ... Further arguments passed to 'hrm_test' will be ignored
#' @return Returns a data.frame with mean and 1-alpha confidence interval for each factor combintation
#' @example R/example_confint.txt
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
      object$data[,prank:= 1/(dim(object$data)[1])*(psrank(object$data[[as.character(object$formula[[2]])]], object$data[, grouping]) - 1/2)]
    } else {
      setDT(object$data)
      object$data[,prank:= 1/(dim(object$data)[1])*(rank(object$data[[as.character(object$formula[[2]])]], ties.method="average") - 1/2)]
    }
    new_formula <- as.formula(paste("prank ~", split(as.character(object$formula), "~")[[1]][3]))
    output <- as.data.frame(summaryBy(new_formula, data = object$data, FUN = mean))
    
    for(i in 1:m) {
      c <- c*0
      c[i] <- 1
      sdi <- sqrt(t(c)%*%object$var%*%c)
      CI_lower[i] <- output[i, ncol] - qt(1-alpha/2, ss[i])*sdi
      CI_upper[i] <- output[i, ncol] + qt(1-alpha/2, ss[i])*sdi
    }
  } else {
    output <- as.data.frame(output)
    for(i in 1:m) {
      c <- c*0
      c[i] <- 1
      sdi <- sqrt(t(c)%*%object$var%*%c)
      CI_lower[i] <- output[i, ncol] - qt(1-alpha/2, ss[i])*sdi
      CI_upper[i] <- output[i, ncol] + qt(1-alpha/2, ss[i])*sdi
    }
  }
  
  output$CI_lower <- CI_lower
  output$CI_upper <- CI_upper
  
  return(output)
}
