library(matrixcalc)
library(MASS)
library(tidyr)
library(HRM)


pseudorank <- function(liste){
  # defining variables und preparing data
  a <- length(liste)
  g <- list(0)
  n <- vapply(liste, FUN=nrow, FUN.VALUE = 0L  )
  g <- lapply(liste, FUN=vec)
  gg <- unlist(g)
  t <- ncol(liste[[1]])

  # defining normalized ecdf's
  F <- function(x, h, n, i){
    return( 1/2*1/n[i]*(sum(h[[i]]<x) + sum(h[[i]]<=x)  ) )
  }

  # unweighted average ecdf
  G <- function(x, h, n, ngroups) {
    return(1/ngroups*sum(vapply(1:ngroups, FUN=F, FUN.VALUE = 1, h=g, n=n, x=x)))
  }

  # calculate pseudoranks
  pseudoranks <- sum(n)*vapply(gg, FUN = G, h = g, n = n, ngroups = a, FUN.VALUE = vector(mode="double", length=1L)  ) + 1/2

  # put all observations from one subject in one row
  result <- NULL
  for(i in 1:a){
    j <- sum(n[0:(i-1)])*t+1
    jj <- sum(n[0:i])*t
    result <- rbind(result, matrix(pseudoranks[j:jj], nrow = n[i], ncol=t, byrow = FALSE) )
  }
  return(result)
}

I = function(size){
  return (diag(rep(1,size)))
}

J = function(size){
  return (rep(1,size)%*%t(rep(1,size)))
}

P = function(size){
  return (I(size)-J(size)*1/size)
}


.E1 = function(n,i, M) {
  return ((n[i]*(n[i]-1))/((n[i]-2)*(n[i]+1))*(matrix.trace(M)^2-2/n[i]*matrix.trace(M%*%M)))
}
.E2 = function(n,i, M) {
  return ((n[i]-1)^2/((n[i]-2)*(n[i]+1))*(matrix.trace(M%*%M)-1/(n[i]-1)*matrix.trace(M)^2))
}

# .E1 = function(n,i, M) {
#   return ((matrix.trace(M)^2))
# }
# .E2 = function(n,i, M) {
#   return ((matrix.trace(M%*%M)))
# }


# Setting:
k = 3 # groups
p = 3 # variables
n = c(5,7,5) # sample sizes
N = sum(n)
alpha=0.05

t <- 30
tt = c(2,10,50,100,300,400)



sim = 10*10^3
g11 = rep(0,sim)

alternative = seq(0, 0, 0.2)
reject = rep(0, length(tt))
reject1 = rep(0, length(tt))
reject11 = rep(0, length(tt))
reject2 = rep(0, length(tt))
reject3 = rep(0, length(tt))
reject31 = rep(0, length(tt))
reject32 = rep(0, length(tt))

start = Sys.time()
#set.seed(12345)
for(l in 1:length(tt)){
  print(tt[l])
  y <- list(reject2, reject3, tt)
  set.seed(tt[l])
  #save(file="result_normal_ind_10_200.RData", y)
  for(j in 1:sim){
    if(j%%100==0){
      print(paste(j,"/",sim))
      print(reject[l]/j)
      print(reject2[l]/j)
      print(reject3[l]/j)
    }
    t = 50 #tt[l]
    n <- c(20,20,20)
    #C = J(k)/k
    C = P(k)/(k-1)
    #D = cbind(I(t-1), rep(-1,t-1))
    D = t(as.matrix(rep(1,t)))*1/t

    #gg <- kronecker(rep(1,n[1]), t(seq(from=1,to=0,length.out = t)))


    a =  matrix(kronecker( rnorm(n[[1]], sd=1)  , t(rep(1,t))  ),ncol=t,nrow=n[1] )
    #a <- a*gg
    group1p1 =  matrix( rnorm(n[[1]]*t,mean = 1, sd=2) , ncol=t,nrow=n[1])       + a
    group1p2 =  matrix(( rnorm(n[[1]]*t, mean = 0, sd = 1) ), ncol=t,nrow=n[1] ) + a
    group1p3 =  matrix(( rnorm(n[[1]]*t, mean = 10, sd = 1) ), ncol=t,nrow=n[1] ) + a

    #gg <- kronecker(rep(1,n[2]), t(seq(from=1,to=0,length.out = t)))
    b =  matrix(kronecker(rnorm(n[[2]],sd=1 )  , t(rep(1,t))  ), ncol=t,nrow=n[2])#*gg
    group2p1 =  matrix(( rnorm(n[[2]]*t, mean = 1, sd=1) ), ncol=t,nrow=n[2])   + b
    group2p2 =  matrix(( rnorm(n[[2]]*t, mean = 0, sd = 2) ), ncol=t,nrow=n[2]) + b
    group2p3 =  matrix(( rnorm(n[[2]]*t, mean = 10,sd = 1) ), ncol=t, nrow=n[2]) + b

    #gg <- kronecker(rep(1,n[3]), t(seq(from=1,to=0,length.out = t)))
    c =  matrix(kronecker(rnorm(n[[3]],sd=1)  , t(rep(1,t))  ), ncol=t,nrow=n[3])#*gg
    group3p1 =  matrix(( rnorm(n[[3]]*t, mean = 1, sd=1) ), ncol=t,nrow=n[3])   + c
    group3p2 =  matrix(( rnorm(n[[3]]*t, mean = 0, sd = 1) ), ncol=t,nrow=n[3]) + c
    group3p3 =  matrix(( rnorm(n[[3]]*t,mean = 10, sd = 1.4) ), ncol=t, nrow=n[3]) + c

    p1=(rbind(group1p1[,],group2p1[,],group3p1[,]))
    p2=(rbind(group1p2[,],group2p2[,],group3p2[,]))
    p3=(rbind(group1p3[,],group2p3[,],group3p3[,]))



    R = data.frame(1:(p*t))
    for(i in 1:sum(n)){
      R[,i] = vec(rbind(p1[i,],p2[i,],p3[i,]))
    }
    R  = as.matrix(R) # in columns are the Y_ik
    Z = kronecker(D, I(p))%*%R
    #Z = R
    Z = as.matrix(Z)



    df <- data.frame(response = vec(as.matrix(R)), subject = gl(sum(n),t*p),
                     group = as.factor(c(rep(1,n[1]*t*p),rep(2,n[2]*t*p),rep(3,n[3]*p*t))),
                     time = as.factor(rep(1:t,p*sum(n))),
                     variable = gl(p,t*sum(n)))

    # hrm.mv.internal(df, "group" , "time", "subject", "response", "variable", C, D )
    # hrm.mv.1w.1f(df, "group" , "time", "subject", "response", "variable" )
    hrm_test(response~group*time, subject=subject, variable=variable, data = df)

        N <- sum(n)
    S1 = var(t(Z[,1:n[[1]] ]))
    S2 = var(t(Z[,(n[[1]]+1):(n[[1]]+n[[2]]) ]))
    S3 = var(t(Z[,((n[[1]]+n[[2]]+1):N) ]))
    S  = 1/k*(1/n[1]*S1 + 1/n[2]*S2 + 1/n[3]*S3)

    Zbar1 = as.matrix(1/n[1]*rowSums(Z[,1:n[1]]))
    Zbar2 = as.matrix(1/n[2]*rowSums(Z[,(n[1]+1):(n[1]+n[2])]))
    Zbar3 = as.matrix(1/n[3]*rowSums(Z[,(n[1]+n[2]+1):N]))
    Zbar = as.matrix(cbind(Zbar1, Zbar2, Zbar3))

    Q <- data.frame(Q1 = rep(0,k), Q2 = rep(0,k))
    Q[1,] <- calcUCpp(t(Z[,1:(n[1])]),n[1],I(dim(S1)[1]),colMeans(t(Z[,1:(n[1])])))
    Q[2,] <- calcUCpp(t(Z[,(n[1]+1):(n[1]+n[2])]),n[2],I(dim(S1)[1]),colMeans(t(Z[,(n[1]+1):(n[1]+n[2])])))
    Q[3,] <- calcUCpp(t(Z[,(n[1]+n[2]+1):(n[1]+n[2]+n[3])]),n[3],I(dim(S1)[1]),colMeans(t(Z[,(n[1]+n[2]+1):(n[1]+n[2]+n[3])])))

    # .E1 = function(n,i, M) {
    #   return ((n[i]*(n[i]-1))/((n[i]-2)*(n[i]+1))*(matrix.trace(M)^2-2/n[i]*matrix.trace(M%*%M)))
    # }
    # .E2 = function(n,i, M) {
    #   return ((n[i]-1)^2/((n[i]-2)*(n[i]+1))*(matrix.trace(M%*%M)-1/(n[i]-1)*matrix.trace(M)^2))
    # }
    .E1 = function(n,i, M) {
      return (Q[i,1])
    }
    .E2 = function(n,i, M) {
      return (Q[i,2])
    }

    #numerator = 1/n[1]^2*(.E2(n,1,S1) + .E1(n,1,S1))+1/n[2]^2*(.E2(n,2,S2) + .E1(n,2,S2))+2*1/n[1]*1/n[2]*(matrix.trace(S1%*%S2) + matrix.trace(S1)*matrix.trace(S2))
    #C = diag(c(1,1,1))
    fHdenominator = 1/n[1]^2*C[1,1]^2*(.E2(n,1,S1)+.E1(n,1,S1))+
      1/n[2]^2*C[2,2]^2*(.E2(n,2,S2)+.E1(n,2,S2))+
      1/n[3]^2*C[3,3]^2*(.E2(n,3,S3)+.E1(n,3,S3))+
      2*1/n[1]*1/n[2]*C[1,2]*C[2,1]*(matrix.trace(S1%*%S2)+matrix.trace(S1)*matrix.trace(S2)) +
      2*1/n[1]*1/n[3]*C[1,3]*C[3,1]*(matrix.trace(S1%*%S3)+matrix.trace(S1)*matrix.trace(S3)) +
      2*1/n[2]*1/n[3]*C[2,3]*C[3,2]*(matrix.trace(S2%*%S3)+matrix.trace(S2)*matrix.trace(S3))

    numeratorH = 1/n[1]^2*C[1,1]^2*(.E2(n,1,S1) + .E1(n,1,S1))+
      1/n[2]^2*C[2,2]^2*(.E2(n,2,S2) + .E1(n,2,S2))+
      1/n[3]^2*C[3,3]^2*(.E2(n,3,S3) + .E1(n,3,S3))  +
      2*1/n[1]*1/n[2]*C[1,1]*C[2,2]*(matrix.trace(S1%*%S2) + matrix.trace(S1)*matrix.trace(S2)) +
      2*1/n[1]*1/n[3]*C[1,1]*C[3,3]*(matrix.trace(S1%*%S3) + matrix.trace(S1)*matrix.trace(S3)) +
      2*1/n[2]*1/n[3]*C[2,2]*C[3,3]*(matrix.trace(S2%*%S3) + matrix.trace(S2)*matrix.trace(S3))

    fGdenominator = 1/k^2*(1/n[1]^2*1/(n[1]-1)^2*(n[1]*(1-1/n[1])^2 + n[1]*(n[1]-1)*1/n[1]^2)*(.E2(n,1,S1)+.E1(n,1,S1))+
                             1/n[2]^2*1/(n[2]-1)^2*(n[2]*(1-1/n[2])^2 + n[2]*(n[2]-1)*1/n[2]^2)*(.E2(n,2,S2)+.E1(n,2,S2)) +
                             1/n[3]^2*1/(n[3]-1)^2*(n[3]*(1-1/n[3])^2 + n[3]*(n[3]-1)*1/n[3]^2)*(.E2(n,3,S3)+.E1(n,3,S3)) )

    # fGdenominator = 1/k^2*(1/n[1]^2*1/(n[1]-1)*(.E2(n,1,S1)+.E1(n,1,S1))+
    #                          1/n[2]^2*1/(n[2]-1)*(.E2(n,2,S2)+.E1(n,2,S2)) +
    #                          1/n[3]^2*1/(n[3]-1)*(.E2(n,3,S3)+.E1(n,3,S3)) )

    numerator = 1/n[1]^2*(.E2(n,1,S1)+.E1(n,1,S1))+
      1/n[2]^2*(.E2(n,2,S2)+.E1(n,2,S2))+
      1/n[3]^2*(.E2(n,3,S3)+.E1(n,3,S3))+
      2*1/n[1]*1/n[2]*(matrix.trace(S1)*matrix.trace(S2)+matrix.trace(S1%*%S2))+
      2*1/n[1]*1/n[3]*(matrix.trace(S1)*matrix.trace(S3)+matrix.trace(S1%*%S3))+
      2*1/n[2]*1/n[3]*(matrix.trace(S2)*matrix.trace(S3)+matrix.trace(S2%*%S3))
    numerator <- numerator*1/k^2




    # numerator = 1/n[1]^2*(.E2(n,1,S1) + .E1(n,1,S1))+1/n[2]^2*(.E2(n,2,S2) + .E1(n,2,S2))+1/n[3]^2*(.E2(n,3,S3) + .E1(n,3,S3))  +2*1/n[1]*1/n[2]*(matrix.trace(S1%*%S2) + matrix.trace(S1)*matrix.trace(S2)) + 2*1/n[1]*1/n[3]*(matrix.trace(S1%*%S3) + matrix.trace(S1)*matrix.trace(S3)) + 2*1/n[2]*1/n[3]*(matrix.trace(S2%*%S3) + matrix.trace(S2)*matrix.trace(S3))

    fH = numeratorH/fHdenominator
    fG <- numerator/fGdenominator

    #C = P(k)*1/(k-1)
    H = fH*Zbar%*%C%*%t(Zbar)
    G = fG*S
    d =  p*matrix.rank(D%*%t(D))
    TD = sqrt(d)*fG*matrix.trace(H)/matrix.trace(G) - sqrt(d)*fH
    #TDsd1 = sqrt(2*(fH*fG+fH^2)*matrix.trace(S%*%S)/d)/(matrix.trace(S)/d*sqrt(fG))
    TDsd1 = (2*(fH)*matrix.trace(S%*%S)/d-matrix.trace(S)^2/(fG*d))/(matrix.trace(S)^2/d^2)*(1+fH/fG)*fG^2/(fG^2+fG-2)

    TDsdalt = (2*(fH)*matrix.trace(S%*%S)/d-matrix.trace(S)^2/(fG*d))/(matrix.trace(S)^2/d^2)*fG^2/(fG^2+fG-2)

    TDsrv = TD/sqrt(TDsd1)
    if(abs(TDsrv)>=qnorm(1-0.05/2)){reject3[l] = reject3[l] + 1}

    TDsrvalt = TD/sqrt(TDsdalt)
    if(abs(TDsrvalt)>=qnorm(1-0.05/2)){reject[l] = reject[l] + 1}

    TDsdN <- 1/n[1]^2*.E2(n,1,S1)+
      1/n[2]^2*.E2(n,2,S2)+
      1/n[3]^2*.E2(n,3,S3)+
      2*1/n[1]*1/n[2]*(matrix.trace(S1%*%S2)) +
      2*1/n[1]*1/n[3]*(matrix.trace(S1%*%S3)) +
      2*1/n[2]*1/n[3]*(matrix.trace(S2%*%S3))
    TDsdN <- sqrt(2*fH*TDsdN*d)

    TDsdD <- 1/n[1]^2*.E1(n,1,S1)+
      1/n[2]^2*.E1(n,2,S2)+
      1/n[3]^2*.E1(n,3,S3)+
      2*1/n[1]*1/n[2]*matrix.trace(S1)*matrix.trace(S2)+
      2*1/n[1]*1/n[3]*matrix.trace(S1)*matrix.trace(S3)+
      2*1/n[2]*1/n[3]*matrix.trace(S2)*matrix.trace(S3)
    TDsdD <- sqrt(TDsdD)

    TDsd <- TDsdN/(matrix.trace(S)/d*k)*sqrt(1+fH/fG)
    TDsd2 <- TDsdN/TDsdD*sqrt(1+fH/fG)

    # TDu = TD/TDsd
    # if(abs(TDu)>=qnorm(1-0.05/2)){reject[l] = reject[l] + 1}
    TDu2 = TD/TDsd2
    if(abs(TDu2)>=qnorm(1-0.05/2)){reject2[l] = reject2[l] + 1}




    # M <- fG - 0.5*(d+fH+1)
    # gamma <- 1/48*d*fH*(d^2+fH^2-5)
    # f <- d*fH
    # prob <- function(x) {
    #   pchisq(x, df = f) + gamma/M^2*(pchisq(x,df=f+4) - pchisq(x,df=f))-0.95
    # }
    # crit <- uniroot(prob, interval = c(0,200))$root
    # if(-M*log(Lambda)>=crit){
    #   reject1 <- reject1 + 1
    # }

    # Muirhead verwirft immer dann, wenn auch Rao F verwirft
    # if(-M*log(Lambda)>=crit | test >= qf(1-alpha, v1, v2)){
    #   reject <- reject + 1
    # }
    #
    # if(-M*log(Lambda)>=crit & test < qf(1-alpha, v1, v2)){
    #   print("...")
    # }

    # large fG, small p , small a
    # Bathke Harrar (2008), p.149
    # Anderson (2003), p.332
    # nu = d*fH
    # qWilks = qchisq(1-alpha,nu) + 1/fG*( (d-fH+1)/2 - 0 )*qchisq(1-alpha,nu) + k*(d+fH+1)/(2*(nu+2))*qchisq(1-alpha,nu)^2*1/fG
    # if( matrix.trace(1/(2-1)*H%*%solve(G,tol=1e-30)) >= qWilks/(fG+0)){
    #   reject[l] = reject[l] + 1
    # }
    # qWilks = qchisq(1-alpha,nu) + 1/fG*( (d-fH+1)/2 + (d-fH+1)/2 )*qchisq(1-alpha,nu) + k*(d+fH+1)/(2*(nu+2))*qchisq(1-alpha,nu)^2*1/fG
    # if( -log(det(G%*%solve(G+H,tol=1e-30))) >= qWilks/(fG-(d-fH+1)/2*1/fG)){
    #   reject1[l] = reject1[l] + 1
    # }
    # qWilks = qchisq(1-alpha,nu) + 1/fG*( (d-fH+1)/2 - (fH-1) )*qchisq(1-alpha,nu) + k*(d+fH+1)/(2*(nu+2))*qchisq(1-alpha,nu)^2*1/fG
    # if( (matrix.trace(H%*%solve(G+H,tol=1e-30))) >= qWilks/(fG +(fH-1)*1/fG)){
    #   reject11[l] = reject11[l] + 1
    # }
    #

    # alternative degrees of freedom
    # fHdenominator = 1/n[1]^2*C[1,1]^2*(.E21(n,1,S1)+.E11(n,1,S1))+1/n[2]^2*C[2,2]^2*(.E21(n,2,S2)+.E11(n,2,S2))+1/n[3]^2*C[3,3]^2*(.E21(n,3,S3)+.E11(n,3,S3))+2*1/n[1]*1/n[2]*C[1,2]*C[2,1]*(matrix.trace(S1%*%S2)+matrix.trace(S1)*matrix.trace(S2)) + 2*1/n[1]*1/n[3]*C[1,3]*C[3,1]*(matrix.trace(S1%*%S3)+matrix.trace(S1)*matrix.trace(S3)) + 2*1/n[2]*1/n[3]*C[2,3]*C[3,2]*(matrix.trace(S2%*%S3)+matrix.trace(S2)*matrix.trace(S3))
    # fGdenominator = (1/n[1]^2*1/(n[1]-1)*(.E21(n,1,S1)+.E11(n,1,S1))+1/n[2]^2*1/(n[2]-1)*(.E21(n,2,S2)+.E11(n,2,S2)) + 1/n[3]^2*1/(n[3]-1)*(.E21(n,3,S3)+.E11(n,3,S3)) )
    #
    # numeratorH = 1/n[1]^2*C[1,1]^2*(.E21(n,1,S1) + .E11(n,1,S1))+1/n[2]^2*C[2,2]^2*(.E21(n,2,S2) + .E11(n,2,S2))+1/n[3]^2*C[3,3]^2*(.E21(n,3,S3) + .E11(n,3,S3))  +2*1/n[1]*1/n[2]*C[1,1]*C[2,2]*(matrix.trace(S1%*%S2) + matrix.trace(S1)*matrix.trace(S2)) + 2*1/n[1]*1/n[3]*C[1,1]*C[3,3]*(matrix.trace(S1%*%S3) + matrix.trace(S1)*matrix.trace(S3)) + 2*1/n[2]*1/n[3]*C[2,2]*C[3,3]*(matrix.trace(S2%*%S3) + matrix.trace(S2)*matrix.trace(S3))
    # numerator = 1/n[1]^2*(.E21(n,1,S1) + .E11(n,1,S1))+1/n[2]^2*(.E21(n,2,S2) + .E11(n,2,S2))+1/n[3]^2*(.E21(n,3,S3) + .E11(n,3,S3))  +2*1/n[1]*1/n[2]*(matrix.trace(S1%*%S2) + matrix.trace(S1)*matrix.trace(S2)) + 2*1/n[1]*1/n[3]*(matrix.trace(S1%*%S3) + matrix.trace(S1)*matrix.trace(S3)) + 2*1/n[2]*1/n[3]*(matrix.trace(S2%*%S3) + matrix.trace(S2)*matrix.trace(S3))
    #
    # fH = numeratorH/fHdenominator
    # fG =numerator/fGdenominator
    #
    # H = fH*Zbar%*%C%*%t(Zbar)
    # G = fG*S
    #
    # Lambda = 1/det(I(dim(H)[1])+H%*%solve(G,tol=1e-30))
    #
    # m1 = fH
    # m2 = fG
    # d =  p*matrix.rank(D%*%t(D)) # dim(Zbar)[1]
    # if(d>fG){
    #   d=fG
    #   print("v2<0")
    #   print(j)
    # }
    # a = sqrt((d^2+m1^2-5)/((m1*d)^2-4))
    # if((d^2+m1^2-5) < 0){
    #   stop("<0")
    # }
    # v1 = m1*d
    # v2 = a^(-1)*(m2-1/2*(d-m1+1))-1/2*(m1*d-2)
    # # finite sample approximation with an F(v1, v2) distribution
    # test = (Lambda^(-a)-1)*v2/v1
    # if(v2 < 0){
    #   print("v2<0")
    # }
    # if( test >= qf(1-alpha, v1, v2)){
    #   reject31[l] = reject31[l] + 1
    # }





    # alternative degrees of freedom 2
    # fHdenominator = 1/n[1]^2*C[1,1]^2*(.E22(n,1,S1)+.E12(n,1,S1))+1/n[2]^2*C[2,2]^2*(.E22(n,2,S2)+.E12(n,2,S2))+1/n[3]^2*C[3,3]^2*(.E22(n,3,S3)+.E12(n,3,S3))+2*1/n[1]*1/n[2]*C[1,2]*C[2,1]*(matrix.trace(S1%*%S2)+matrix.trace(S1)*matrix.trace(S2)) + 2*1/n[1]*1/n[3]*C[1,3]*C[3,1]*(matrix.trace(S1%*%S3)+matrix.trace(S1)*matrix.trace(S3)) + 2*1/n[2]*1/n[3]*C[2,3]*C[3,2]*(matrix.trace(S2%*%S3)+matrix.trace(S2)*matrix.trace(S3))
    # fGdenominator = (1/n[1]^2*1/(n[1]-1)*(.E22(n,1,S1)+.E12(n,1,S1))+1/n[2]^2*1/(n[2]-1)*(.E22(n,2,S2)+.E12(n,2,S2)) + 1/n[3]^2*1/(n[3]-1)*(.E22(n,3,S3)+.E12(n,3,S3)) )
    #
    # numeratorH = 1/n[1]^2*C[1,1]^2*(.E22(n,1,S1) + .E12(n,1,S1))+1/n[2]^2*C[2,2]^2*(.E22(n,2,S2) + .E12(n,2,S2))+1/n[3]^2*C[3,3]^2*(.E22(n,3,S3) + .E12(n,3,S3))  +2*1/n[1]*1/n[2]*C[1,1]*C[2,2]*(matrix.trace(S1%*%S2) + matrix.trace(S1)*matrix.trace(S2)) + 2*1/n[1]*1/n[3]*C[1,1]*C[3,3]*(matrix.trace(S1%*%S3) + matrix.trace(S1)*matrix.trace(S3)) + 2*1/n[2]*1/n[3]*C[2,2]*C[3,3]*(matrix.trace(S2%*%S3) + matrix.trace(S2)*matrix.trace(S3))
    # numerator = 1/n[1]^2*(.E22(n,1,S1) + .E12(n,1,S1))+1/n[2]^2*(.E22(n,2,S2) + .E12(n,2,S2))+1/n[3]^2*(.E22(n,3,S3) + .E12(n,3,S3))  +2*1/n[1]*1/n[2]*(matrix.trace(S1%*%S2) + matrix.trace(S1)*matrix.trace(S2)) + 2*1/n[1]*1/n[3]*(matrix.trace(S1%*%S3) + matrix.trace(S1)*matrix.trace(S3)) + 2*1/n[2]*1/n[3]*(matrix.trace(S2%*%S3) + matrix.trace(S2)*matrix.trace(S3))
    #
    # fH = numeratorH/fHdenominator
    # fG =numerator/fGdenominator
    #
    # H = fH*Zbar%*%C%*%t(Zbar)
    # G = fG*S
    #
    # Lambda = 1/det(I(dim(H)[1])+H%*%solve(G,tol=1e-30))
    #
    # m1 = fH
    # m2 = fG
    # d =  p*matrix.rank(D%*%t(D)) # dim(Zbar)[1]
    # if(d>fG){
    #   d=fG
    #   print("v2<0")
    #   print(j)
    # }
    # a = sqrt((d^2+m1^2-5)/((m1*d)^2-4))
    # v1 = m1*d
    # v2 = a^(-1)*(m2-1/2*(d-m1+1))-1/2*(m1*d-2)
    # # finite sample approximation with an F(v1, v2) distribution
    # test = (Lambda^(-a)-1)*v2/v1
    # if(v2 < 0){
    #   print("v2<0")
    # }
    # if( test >= qf(1-alpha, v1, v2)){
    #   reject32[l] = reject32[l] + 1
    # }



  }
}

Sys.time()-start
#reject/sim # chi^2 Approximation LH
# reject1/sim # chi^2 Approximation LR
# reject11/sim # chi^2 Approximation BNP
reject2/sim # Wilk
reject3/sim # Fuji
#reject31/sim # Rao F, U Sch?tzer
#reject32/sim # Rao F, Plug-In Sch?tzer
n
p*tt

# y <- list(reject, reject2, reject3, tt)
# save(file="result_normal_ind_10_400_101510.RData", y)
