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
par(mfrow = c(1,2))


#############################Training VS Test Set##########################################################################################
data(age.income)
sample <- sample.split(age.income$log.income, SplitRatio = 0.7)
train  <- subset(age.income, sample == TRUE)
test   <- subset(age.income, sample == FALSE)
train.input <- train$age
train.output <- train$log.income
test.input <-  test$age 
test.output <- test$log.income


#############################Sim Parameters##########################################################################################

N <- 500 #discretization
nbsim <- 100 #nb of simulations
n <- length(train.input)
n_t <- length(test.input)
tt <- 60 ### Theta
sig <- sd(train.output) # std of noise
u <- seq(20,70,length=N)
t <- c(tt,sig)

###############################################TRAINING####################################################

####################################Function
k <- function(x,xp,theta){
  (1+(sqrt(5)*abs(x-xp)/theta)+(5*(x-xp)^2)/(3*theta^2))*exp(-sqrt(5)*(abs(x-xp))/theta)
}

#mean((krigMean(test.input,t[1],t[2])-test.output)^2) 
####################################Convariance Matrix & Function


kx <- matrix(NA,nrow=N,ncol=n)
for(j in 1 : n ){
  for( i in 1 :N ){
    kx[i,j]=k(u[i],train.input[j], t[1])
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
    Kxx[i,j]=k(train.input[i],train.input[j], t[1])
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
    kxx[i,j]=k(u[i],u[j], t[1])
  }
}
t[2]

################################### Conditional Cov and Mean + Simulations
C <- kxx - kx%*%chol2inv(chol(Kxx+ t[2]^2*diag(n)))%*%t(kx) #conditional covariance matrix
zeta <- kx%*%chol2inv(chol(Kxx+ t[2]^2*diag(n)))%*%train.output #conditional mean vector  --> ANYTHING IN THE 
Y <- t(mvrnorm(nbsim,zeta,C)) ### sims
Y 

################################### Kriging : Transforming the mean vector into the mean function 
krigMean <- function(x,theta,sigN){
  kxfct(x,theta)%*%chol2inv(chol(Kxxfct(x,theta)+ sigN^2*diag(n)))%*%train.output
}

#mean((krigMean(test.input,t[1],t[2])-test.output)^2)
################################### test set applied to the mean fct 
Y_t <- krigMean(test.input, t[1],t[2])


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


#####################################Method 1- Optimization##########################################

MSPE1 <- function(t){
  
  mean((krigMean(test.input,t[1],t[2])-test.output)^2)
  
}

res_BFGS <- optim(t, MSPE1 ,method = "BFGS" ) #,control = list(), hessian = FALSE)
res_BFGS
res_SANN <- optim(t, MSPE1 ,method = "SANN", control = list(maxit=500)) #,control = list(), hessian = FALSE)

res_SANN

res_CG <- optim(t, MSPE1 ,method = "CG" ) 
res_CG

#####################################Method 2- Visualization##########################################

################## Optimizing Theta

MSPE2 <- function(t){
  res <- matrix(NA,length(t),2)
  for (i in 1: length(t)){ 
   
    res[i,1] <- t[i]
    res[i,2] <-  mspe(krigMean(test.input, t[i],sig),test.output)
  }
  return(res)
}
TT <- c(30:100)

plot(MSPE2(TT)[,1],MSPE2(TT)[,2], xlab='Theta' , ylab='MSPE')


################## Optimizing Sigma 

MSPE3 <- function(t){
  res <- matrix(NA,length(t),2)
  for (i in 1: length(t)){ 
    
    res[i,1] <- t[i]
    res[i,2] <-  mspe(krigMean(test.input, 60,t[i]),test.output)
  }
  return(res)
}
SS <- seq(0.05,0.5,by=0.01)
plot(SS,MSPE3(SS)[,2], xlab='Sigma' , ylab='MSPE')



