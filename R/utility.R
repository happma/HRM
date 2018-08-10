####################################################################################################################################
### Filename:    utility.R
### Description: trace estimator functions, dual empirical covariance matrix function
###              
###
###
####################################################################################################################################


calcU <- function(X, n, grp, K){
  Q1c <- 0
  Q2c <- 0
  for(i in 1:n[grp]){
    j <- i + 1
    while(j <= n[grp]){
      #ii <- min(setdiff(1:n[grp],c(i,j)))
      #jj <- min(setdiff(1:n[grp],c(i,j,ii)))
      #Z1 <- K%*%(X[[grp]][i,]- X[[grp]][ii,] )
      #Z2 <- K%*%(X[[grp]][j,]- X[[grp]][jj,] )
      #Z1 <- K%*%(X[[grp]][i,]- 1/(n[grp]-1)*(colMeans(X[[grp]])*n[grp]-X[[grp]][j,]) )
      #Z2 <- K%*%(X[[grp]][j,]- 1/(n[grp]-1)*(colMeans(X[[grp]])*n[grp]-X[[grp]][i,]) )
      Z1 <- K%*%(X[[grp]][i,]- colMeans(X[[grp]]) )
      Z2 <- K%*%(X[[grp]][j,]- colMeans(X[[grp]]) )
      Q2c <- Q2c + (t(Z1)%*%Z2)^2
      Q1c <- Q1c + t(Z1)%*%Z1*t(Z2)%*%Z2
      j <- j + 1
    }
  }
  #corr <- (1-1/n[grp])^2+(n[grp]-2)*1/(n[grp]-1)^2
  #corr <- 2
  corr <- (1-1/n[grp])
  Q2c <- Q2c*2/(n[grp]*(n[grp]-1))*1/corr^2 #(1-1/n[grp])^(-2)#1/4
  Q1c <- Q1c*2/(n[grp]*(n[grp]-1))*1/corr^2 #(1-1/n[grp])^(-2)#1/4
  return(c(Q1c,Q2c))
}


calcU_onegroup <- function(X, n, K){
  Q1c <- 0
  Q2c <- 0
  for(i in 1:n[1]){
    j <- i + 1
    while(j <= n[1]){
      #ii <- min(setdiff(1:n[1],c(i,j)))
      #jj <- min(setdiff(1:n[1],c(i,j,ii)))
      #Z1 <- K%*%(X[i,]- X[ii,] )
      #Z2 <- K%*%(X[j,]- X[jj,] )
      Z1 <- K%*%(X[i,]-colMeans(X) )
      Z2 <- K%*%(X[j,]-colMeans(X) )

      Q2c <- Q2c + (t(Z1)%*%Z2)^2
      Q1c <- Q1c + t(Z1)%*%Z1*t(Z2)%*%Z2
      j <- j + 1
    }
  }
  #corr <- (n[1]^2-2*n[1]+2)/n[1]^2
  Q2c <- Q2c*2/(n[1]*(n[1]-1))*(1-1/n[1])^(-2)#1/4
  Q1c <- Q1c*2/(n[1]*(n[1]-1))*(1-1/n[1])^(-2)#1/4
  return(c(Q1c,Q2c))
}


#' Unbiased estimator
#'
#' @param i group index
#' @param M a matrix
#' @param n vector of sample size
#' @keywords internal
.E1 <- function(n,i, M, nonparametric, Q) {
  
  trace_estimator <- ifelse(nonparametric, Q[i,1], (n[i]*(n[i]-1))/((n[i]-2)*(n[i]+1))*(matrix.trace(M)^2-2/n[i]*matrix.trace(M%*%M)))
  return (trace_estimator)
}
#' Unbiased estimator
#'
#' @param i group index
#' @param M a matrix
#' @param n vector of sample size
#' @keywords internal
.E2 <- function(n,i, M, nonparametric, Q) {
  trace_estimator <- ifelse(nonparametric, Q[i,2], (n[i]-1)^2/((n[i]-2)*(n[i]+1))*(matrix.trace(M%*%M)-1/(n[i]-1)*matrix.trace(M)^2))
  return (trace_estimator)
}
#' Unbiased estimator
#' 
#' @param i group index
#' @param M a matrix
#' @keywords internal
.E3 <- function(M_i, M_j) {
  return (matrix.trace(M_i)*matrix.trace(M_j))
}
#' Unbiased estimator
#' 
#' @param i group index
#' @param M a matrix
#' @keywords internal
.E4 <- function(M_i,M_j) {
  return (matrix.trace(M_i%*%M_j))
}


#' Function for the output: significant p-values have on or more stars
#' 
#' @param value p-value
#' @keywords internal
.hrm.sigcode <- function(value) {
  
  char <- ""
  if(value <= 0.1 & value > 0.05) { char = "."}
  if(value <= 0.05 & value > 0.01) {char = '*'}
  if(value <= 0.01 & value > 0.001) {char = "**"}
  if(value <= 0.001 & value >= 0) {char = "***"}
  return (char)
}


#' Function for the indentity matrix
#' 
#' @param size dimension of the matrix
#' @keywords internal
I <- function(size){
  return (diag(rep(1,size)))
}
#' Function for a matrix with entries 1
#' 
#' @param size dimension of the matrix
#' @keywords internal

J <- function(size){
  return (rep(1,size)%*%t(rep(1,size)))
}
#' Function for the centering matrix
#' 
#' @param size dimension of the matrix
#' @keywords internal
P <- function(size){
  return (I(size)-J(size)*1/size)
}

#' Function for the dual empirical matrix
#' 
#' @param Data data.frame
#' @param B not used
#' @keywords internal
DualEmpirical <- function(Data, B){
  n <- dim(Data)[1]
  B <- B(dim(Data)[2]) # B = e.g. t(J_d)%*%J_d
  return(1/(n-1)*P(n)%*%Data%*%B%*%t(Data)%*%P(n))
}
#' Function for the dual empirical matrix
#' 
#' @param Data data.frame
#' @param B part of the hypothesis matrix
#' @keywords internal
DualEmpirical2 <- function(Data, B){
  n <- dim(Data)[1]
  return(1/(n-1)*P(n)%*%Data%*%B%*%t(Data)%*%P(n))
}