################################# Matern 5/2 covariance function  #####################################################################

library(SemiPar)
library(caTools)
library(MASS)
library(ggplot2)
library(Rcpp)
library(SCtools)
library(cvTools)
library(robustbase)
set.seed(10)
##par(mfrow = c(1,2))


#############################Training VS Test Set##########################################################################################
data(age.income)
sample <- sample.split(age.income$log.income, SplitRatio = 0.7)
train  <- subset(age.income, sample == TRUE)
test   <- subset(age.income, sample == FALSE)
?sample
train.input <- train$age
train.output <- train$log.income
test.input <-  test$age
test.output <- test$log.income


#############################Sim Parameters##########################################################################################

N <- 500 #discretization
nbsim <- 100 #nb of simulations
n <- length(train.input)
n_t <- length(test.input)
tt <- 30 ### Theta
sig <- sd(train.output) # std of noise
u <- seq(20,70,length=N)

###############################################TRAINING####################################################

####################################Function
k <- function(x,xp,theta){
  exp(-abs(x-xp)/theta)
  }


####################################Convariance Matrix & Function
kx <- matrix(NA,nrow=N,ncol=n)
for(j in 1 : n ){
  for( i in 1 :N ){
    kx[i,j]=k(u[i],train.input[j], tt)
  }
}

kxfct <- function(x,theta){
  V <- matrix(NA,length(x),length(train.input))
  for(i in 1 : length(x)){
    for( j in 1 :length(train.input) ){
      V[i,j]=k(x[i],train.input[j], theta)
    }
  }
  return(V)
}

####################################Kxx Matrix & Function

Kxx <- matrix(NA,n,n)
for(j in 1 : n ){
  for( i in 1 :n ){
    Kxx[i,j]=k(train.input[i],train.input[j], tt)
  }
}

Kxxfct <- function(x,theta){
  W <- matrix(NA,length(train.input),length(train.input))
  for(i in 1 : length(train.input)){
    for( j in 1 :length(train.input) ){
      W[i,j]=k(train.input[i],train.input[j], theta)
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
zeta <- kx%*%chol2inv(chol(Kxx+ sig^2*diag(n)))%*%train.output #conditional mean vector  --> ANYTHING IN THE 
Y <- t(mvrnorm(nbsim,zeta,C)) ### sims
Y 


################################### Kriging : Transforming the mean vector into the mean function 
krigMean <- function(x,theta,sigN){
  kxfct(x,theta)%*%chol2inv(chol(Kxxfct(x,theta)+ sigN^2*diag(n)))%*%train.output
}

################################### test set applied to the mean fct 
Y_t <- krigMean(test.input, tt,sig)


matplot(u, Y, type='l', lwd=1,col='gray',lty=1,ylab='', xlab='')
lines(u,zeta,col="black",lwd=2)
lines(test.input,Y_t, type='l', lwd=1,col='red',lty=1,ylab='', xlab='')
title(xlab="x",ylab="Y(x)",line=2)

points(train.input,train.output,pch=19, col='black')
points(test.input,test.output,pch=19, col='red')





###############################################Evaluation####################################################
### User Defined
MSPE <- mean((Y_t -test.output)^2)
MSPE
### Pre-Defined
mspe(Y_t,test.output)


