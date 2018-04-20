####################################################################################################################################
### Filename:    f2.R
### Description: Function for calculating the test statistic for one whole- and one subplot factor
###              
###
###
####################################################################################################################################




#' Test for interaction of factor A and B
#' 
#' @param X dataframe containing the data in the long table format
#' @param alpha alpha level used for the test
#' @param group column name of the data frame X specifying the groups
#' @param factor1 column name of the data frame X of the first factor variable
#' @param subject column name of the data frame X identifying the subjects
#' @param data column name of the response variable
#' @param H string specifying the hypothesis
#' @param text a string, which will be printed in the output
#' @return Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
#' @keywords internal
hrm.1w.1f <- function(X, alpha, group , factor1, subject, data, H, text, nonparametric, ranked ){
  
  stopifnot(is.data.frame(X),is.character(subject), is.character(group),is.character(factor1), alpha<=1, alpha>=0, is.logical(nonparametric))
  f <- 0
  f0 <- 0
  crit <- 0
  test <- 0  
  
  group <- as.character(group)
  factor1 <- as.character(factor1)
  subject <- as.character(subject)
  
  X <- as.data.table(X)
  setnames(X, c(data, group, factor1, subject), c("data", "group", "factor1", "subject"))

  a <- nlevels(X[,group])
  d <- nlevels(X[,factor1])
  c <- 1
  n <- table(X[,group])/d
  KGV <- Reduce(Lcm, n)
  lambda <- KGV/n
  
  if(max(lambda) <= 100 & max(n) <= 30 & nonparametric & is.null(ranked)){
    len <- dim(X)[1]
    prData <- list(X,0)
    z <- levels(X[,group])
    
    # amplify data to artificially create balanced groups
    for(i in 1:a){
      prData[[i+1]] <- X[group==z[i]][rep(1:(n[i]*d), each = (lambda[i]-1)), ]
    }
    X <- rbindlist(prData)
    X[,data]<- (rank(X[,data], ties.method = "average")-1/2)*1/(KGV*a*d)

    # select original observations from amplified data
    X <- X[1:len,]
  }
  
  X <- split(X, X[,group], drop=TRUE)
  for(i in 1:a){
    X[[i]] <- X[[i]][ order(subject, factor1), ]
    X[[i]] <- X[[i]][,data]
    X[[i]] <- matrix(X[[i]],ncol=d*c,byrow=TRUE)
    n[i] <- dim(X[[i]])[1]
  }

  if((max(lambda) > 100 | max(n) > 30) & nonparametric & is.null(ranked)){
    X <- pseudorank(X)
    for(i in 1:a){
      X[[i]] <- 1/(sum(n)*d)*(X[[i]] - 1/2)
    }
  }

  if(is.null(ranked)){
    eval.parent(substitute(ranked<-X))
  } else {
    X <- ranked
  }
  
  # group <- as.character(group)
  # factor1 <- as.character(factor1)
  # subject <- as.character(subject)
  # X <- split(X, X[,group], drop=TRUE)
  # a <- length(X)
  # d <- nlevels(X[[1]][,factor1])
  # c <- 1
  # n <- rep(0,a) 
  # 
  # for(i in 1:a){
  #   X[[i]] <- X[[i]][ order(X[[i]][,subject], X[[i]][,factor1]), ]
  #   X[[i]]<-X[[i]][,data]
  #   X[[i]] <- matrix(X[[i]],ncol=d*c,byrow=TRUE)
  #   n[i] <- dim(X[[i]])[1]
  # }
  
  # dataOriginal <<- X
  # 
  # startR <- Sys.time()
  # X <- pseudorank(X)
  # endeR <<- Sys.time() - startR
  # 
  # dataRanks <<-  X
  
  # print("orig")
  # print( sum(var(X[[1]])[1,] ))
  # 
  # Y <- X[[1]]
  # 
  # for(i in 1:(n[1])){
  #   Y[i,] <- Y[i,]-(sum(n)*d+1)/2
  # }
  # 
  # print("centered")
  # print(sum(var(Y)[1,] ))
  


  
  # fehl <- rep(0,a)
  # for(j in 1:1){
  #   g1 <- rep(0, n[j])
  #   for(i in 1:n[j]){
  #     g1[i] <- sum((X[[j]][i,]-1/2)^2)^2
  #   }
  #   fehl[j] <- mean(g1)-matrix.trace(var(X[[j]]-1/2 ))^2 - 2*matrix.trace(var(X[[j]]-1/2)%*%var(X[[j]]-1/2) )
  # }
  # print("fehler")
  # print(fehl[1])

  
  # creating X_bar (list with a entries)
  X_bar <<- as.matrix(vec(sapply(X, colMeans, na.rm=TRUE))) #- (sum(n)*d+1)*1/2

  if(H=="A"){
    K <- 1/d*J(d)
    S <- diag(n)-1/sum(n)*n%*%t(n)
  } else if(H=="Au"){
    K <- 1/d*J(d)
    S <- P(a)
  } else if(H=="B"){
    K <- P(d)
    S <- J(a)
  } else if(H=="AB"){
    K <- P(d)
    S <- P(a)
  }

  # creating dual empirical covariance matrices
  K_AB <- kronecker(S, K)
  V <- lapply(X, DualEmpirical2, B=K)
  
  zaehler <<- t(X_bar)%*%K_AB%*%X_bar
  
  sp1 <<- matrix.trace(V[[1]])
  cov1 <<- t(X[[1]][1,])%*%(X[[1]][1,])
  cov2 <<- t(X[[1]][2,])%*%(X[[1]][2,])
  
  ##########################
  ### U statistics
  #########################
  
  Q = data.frame(Q1 = rep(0,a), Q2 = rep(0,a))
  if(nonparametric){
    for(i in 1:a){
      Q[i,] <- calcU(X,n,i,K)
    }    
  }


  # for(i in 1:n[1]){
  #   j <- i + 1
  #   while(j <= n[1]){
  #     ii <- sample(setdiff(1:n[1],c(i,j)),1)
  #     jj <- sample(setdiff(1:n[1],c(i,j,ii)),1)
  #     Z1 <- K%*%(X[[1]][i,]- X[[1]][ii,] )
  #     Z2 <- K%*%(X[[1]][j,]- X[[1]][jj,] )
  #     # Z1 <- K%*%(X[[1]][i,]-colMeans(X[[1]]) )
  #     # Z2 <- K%*%(X[[1]][j,]-colMeans(X[[1]]) )
  #     Q2[1] <- Q2[1] + (t(Z1)%*%Z2)^2
  #     Q1[1] <- Q1[1] + t(Z1)%*%Z1*t(Z2)%*%Z2
  #     j <- j + 1
  #   }
  # }
  #corr <- (n[1]^2-2*n[1]+2)/n[1]^2
  # Q2[1] <- Q2[1]*2/(n[1]*(n[1]-1))*1/corr[1]^2
  # Q1[1] <- Q1[1]*2/(n[1]*(n[1]-1))*1/corr[1]^2
  # 
  # for(i in 1:n[2]){
  #   j <- i + 1
  #   while(j <= n[2]){
  #     ii <- sample(setdiff(1:n[2],c(i,j)),1)
  #     jj <- sample(setdiff(1:n[2],c(i,j,ii)),1)
  #     Z1 <- K%*%(X[[2]][i,]- X[[2]][ii,] )
  #     Z2 <- K%*%(X[[2]][j,]- X[[2]][jj,] )
  #     # Z1 <- K%*%(X[[2]][i,]-colMeans(X[[2]]) )
  #     # Z2 <- K%*%(X[[2]][j,]-colMeans(X[[2]]) )
  #     Q2[2] <- Q2[2] + (t(Z1)%*%Z2)^2
  #     Q1[2] <- Q1[2] + t(Z1)%*%Z1*t(Z2)%*%Z2
  #     j <- j + 1
  #   }
  # }
  # #corr <- (n[2]^2-2*n[2]+2)/n[2]^2
  # Q2[2] <- Q2[2]*2/(n[2]*(n[2]-1))*1/corr[2]^2
  # Q1[2] <- Q1[2]*2/(n[2]*(n[2]-1))*1/corr[2]^2
  # 
  # for(i in 1:n[3]){
  #   j <- i + 1
  #   while(j <= n[1]){
  #     ii <- sample(setdiff(1:n[3],c(i,j)),1)
  #     jj <- sample(setdiff(1:n[3],c(i,j,ii)),1)
  #     Z1 <- K%*%(X[[3]][i,]- X[[3]][ii,] )
  #     Z2 <- K%*%(X[[3]][j,]- X[[3]][jj,] )
  #     # Z1 <- K%*%(X[[3]][i,]-colMeans(X[[3]]) )
  #     # Z2 <- K%*%(X[[3]][j,]-colMeans(X[[3]]) )
  #     Q2[3] <- Q2[3] + (t(Z1)%*%Z2)^2
  #     Q1[3] <- Q1[3] + t(Z1)%*%Z1*t(Z2)%*%Z2
  #     j <- j + 1
  #   }
  # }
  # #corr <- (n[3]^2-2*n[3]+2)/n[3]^2
  # Q2[3] <- Q2[3]*2/(n[3]*(n[3]-1))*1/corr[3]^2
  # Q1[3] <- Q1[3]*2/(n[3]*(n[3]-1))*1/corr[3]^2


  # .E1 <- function(n,i, M) {
  #   return (Q[i,1])
  # }
  # .E2 <- function(n,i, M) {
  #   return (Q[i,2])
  # }
  
  # .E1 = function(n,i, M) {
  #   return ((matrix.trace(M)^2))
  # }
  # .E2 = function(n,i, M) {
  #   return ((matrix.trace(M%*%M)))
  # }
  
  # .E11 <- function(n,i, M) {
  #   return ((n[i]*(n[i]-1))/((n[i]-2)*(n[i]+1))*(matrix.trace(M)^2-2/n[i]*matrix.trace(M%*%M)))
  # }
  # .E21 <- function(n,i, M) {
  #   return ((n[i]-1)^2/((n[i]-2)*(n[i]+1))*(matrix.trace(M%*%M)-1/(n[i]-1)*matrix.trace(M)^2))
  # }
   
  #################################################################################################

  
  # f
  f_1 <- 0
  f_2 <- 0
  
  for(i in 1:a){
    f_1 <- f_1 + (S[i,i]*1/n[i])^2*.E1(n,i,V[[i]],nonparametric,Q)
    j <- i+1
    while(j<=a){
      f_1 <- f_1 + 2*(S[i,i]*1/n[i])*(S[j,j]*1/n[j])*.E3(V[[i]],V[[j]])
      j <- j+1
    }
  }
  
  for(i in 1:a){
    f_2 <- f_2 + (S[i,i]*1/n[i])^2*.E2(n,i,V[[i]],nonparametric,Q)
    j <- i+1
    while(j<=a){
      f_2 <- f_2 + 2*S[i,j]*S[j,i]*1/(n[i]*n[j])*.E4(1/(n[i]-1)*P(n[i])%*%X[[i]],1/(n[j]-1)*K%*%t(X[[j]])%*%P(n[j])%*%X[[j]]%*%K%*%t(X[[i]])%*%P(n[i]))
      j <- j+1
    }
  }
  
  f <- f_1/f_2
  

  ##################################################################################################
  
  
  
  #################################################################################################
  # f0
  f0_1 <- f_1
  f0_2 <- 0
  
  
  for(i in 1:a){
    f0_2 <- f0_2 + (S[i,i]*1/n[i])^2*1/(n[i]-1)*.E2(n,i,V[[i]],nonparametric,Q)
  }

  #debuging <<-abs(fehl[1])*d/n[1]^4*1/f0_2
  f0 <- f0_1/f0_2
  
  ##################################################################################################

  # critical value
  crit <- qf(1-alpha,f,f0)
  
  # Test
  direct <- direct.sum(1/n[1]*var(X[[1]]),1/n[2]*var(X[[2]]))
  if(a>2){
    for(i in 3:a) {
      direct <- direct.sum(direct, 1/n[i]*var(X[[i]]))
    }
  }

  test <- (t(X_bar)%*%K_AB%*%X_bar)/(t(rep(1,dim(K_AB)[1]))%*%(K_AB*direct)%*%(rep(1,dim(K_AB)[1])))
  p.value <- 1-pf(test,f,f0)
  output <- data.frame(hypothesis=text,df1=f,df2=f0, crit=crit, test=test, p.value=p.value, sign.code=.hrm.sigcode(p.value))
  
  
  return (output)
}

# End ------------------------------------------------------------
