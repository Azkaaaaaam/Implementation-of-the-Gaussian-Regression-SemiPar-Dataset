  
################################# BASIC  #####################################################################


library(SemiPar)
library(caTools)
library(MASS)
library(ggplot2)
library(Rcpp)
library(SCtools)
library(cvTools)
library(robustbase)
set.seed(10)
??gaupro

data(age.income)
input <- age.income$age
output <- age.income$log.income
  summary(exp(output))
  summary(input)
  
  #####################Exploration
#### SIMPLE 
plot(input,output)

  
  
#### methohd 1
pf <- ggplot(age.income, aes(age, log.income)) + geom_point() + xlab("Age") + ylab("log.income")
pf+geom_smooth(method="gam", formula = y ~ s(x, bs = "cs"))





#### method 2

#############################Sim Parameters##########################################################################################

N <- 500 #discretization
nbsim <- 100 #nb of simulations
n <- length(input)
n_t <- length(input)
tt <- 10 ### Theta
sig <- sd(output) # std of noise
u <- seq(20,70,length=N)

  
  


###################################################################################################

####################################Function
k <- function(x,xp,theta){
  (1+(sqrt(3)*abs(x-xp)/theta))*exp(-sqrt(3)*abs(x-xp)/theta)
}


####################################Convariance Matrix & Function

####################################Convariance Matrix & Function
kx <- matrix(NA,nrow=N,ncol=n)
for(j in 1 : n ){
  for( i in 1 :N ){
    kx[i,j]=k(u[i],input[j], tt)
  }
}

kxfct <- function(x,theta){
  V <- matrix(NA,length(x),length(input))
  for(i in 1 : length(x)){
    for( j in 1 :length(input) ){
      V[i,j]=k(x[i],input[j], theta)
    }
  }
  return(V)
}

####################################Kxx Matrix & Function

Kxx <- matrix(NA,n,n)
for(j in 1 : n ){
  for( i in 1 :n ){
    Kxx[i,j]=k(input[i],input[j], tt)
  }
}

Kxxfct <- function(x,theta){
  W <- matrix(NA,length(input),length(input))
  for(i in 1 : length(input)){
    for( j in 1 :length(input) ){
      W[i,j]=k(input[i],input[j], theta)
    }
  }
  return(W)
}


#################################### kxx Function

kxx <- matrix(NA,N,N)
for(j in 1 : N ){
  for( i in 1 :N ){
    kxx[i,j]=k(u[i],u[j], tt)
  }
}

################################### Conditional Cov and Mean + Simulations
C <- kxx - kx%*%chol2inv(chol(Kxx+ sig^2*diag(n)))%*%t(kx) #conditional covariance matrix
zeta <- kx%*%chol2inv(chol(Kxx+ sig^2*diag(n)))%*%output #conditional mean vector  --> ANYTHING IN THE 
Y <- t(mvrnorm(nbsim,zeta,C)) ### sims
Y 



################################### Kriging : Transforming the mean vector into the mean function 
krigMean <- function(x,theta,sigN){
  kxfct(x,theta)%*%chol2inv(chol(Kxxfct(x,theta)+ sigN^2*diag(n)))%*%output
}

################################### Output estimate 
Y_t <- krigMean(input, tt,sig)

q1 <- zeta + qnorm(0.05, 0, sqrt(diag(C)))
q2 <- zeta + qnorm(0.95, 0, sqrt(diag(C)))

matplot(u, Y, type='l', lwd=1,col='gray',lty=1,ylab='', xlab='')
lines(u,zeta,col="black",lwd=2)
lines(input,Y_t, type='l', lwd=1,col='red',lty=1,ylab='', xlab='')
title(xlab="x",ylab="Y(x)",line=2)

points(input,output,pch=19, col='black')
lines(u, q1, lwd=2, lty=2, col=2)
lines(u, q2, lwd=2, lty=2, col=2)

  MSPE <- function(t){
    t[1] <- theta
    t[2] <- sigN
    mean((krigMean(input,t[1],t[2])-output)^2)
  }
  
