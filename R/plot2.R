#' Plotting Profile Curves
#' 
#' @description Plotting profile curves for up to one whole- or subplot-factor
#' @param x An object of class 'HRM' from the function 'hrm.test'
#' @param xlab label of the x-axis of the plot
#' @param ylab label of the y-axis of the plot
#' @param legend logical indicating if a legend should be plotted
#' @param legend.title title of the legend
#' @param ... Further arguments passed to the 'plot' function
#' @example R/example_plot.txt
#' @keywords export
plot.HRM <- function(x, xlab = "dimension", ylab = "mean", legend = TRUE, legend.title = NULL, ...) {
  
  stopifnot(is.logical(legend), is.character(xlab), is.character(ylab), class(x) == "HRM")
  
  if(length(x$factors[[1]]) == 1 & length(x$factors[[2]]) == 1) {
    group <- x$factors[[1]][1]
    factor1 <- x$factors[[2]][1]
    if(group == "none" & factor1 != "none") {
      x$data$group = as.factor(1)
      group <- "group"
    }
    if(group != "none" & factor1 == "none"){
      x$data$factor1 == as.factor(1)
      factor1 <- "factor1"      
      p <- NULL
      stop("At least one subplot-factor is needed!")
    }
    
    pl <- hrm.plot(data = x$data, group = group , factor1 = factor1, subject = x$subject, response = as.character(x$formula[[2]]), xlab=xlab, ylab=ylab, legend = legend, legend.title = legend.title )
    return(pl)
    
  } else if(is.null(x$formula)){
    
    m <- lapply(x$data, colMeans)
    m <- unlist(m)
    a <- length(x$data)
    d <- dim(x$data[[1]])[2]
    
    if(d <= 1){
      stop("At least two measurements per subject are needed!")
    }
    
    means <- data.frame(value = m, group = gl(a, d), dimension = as.factor(rep(1:d, a)))
    
    pl <- ggplot() +
          geom_line(data=means, aes(x=means$dimension, y=means$value,group=means$group,colour=means$group)) +
          geom_point(data=means, aes(x=means$dimension, y=means$value,group=means$group,colour=means$group),size=1.5) +
          xlab(xlab) + 
          ylab(ylab)   
    
    if(!legend){
      pl <- pl + theme(legend.position = "none") 
    } else {
        if(!is.null(legend.title) & is.character(legend.title)){
          pl <- pl + scale_colour_hue(name=legend.title)
        } else {
            pl <- pl + theme(legend.title = element_blank())
        }
        pl <- pl + theme(legend.background = element_rect()) 
    }

    
    return(pl)
    
  } else {
    print("Plot function only supports one whole- and one subplotfactor.")
  }

}