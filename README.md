# Implementation-of-the-Gaussian-Regression-SemiPar-Dataset

## Repository Tour 
- [Basic Algorithm.R](https://github.com/Azkaaaaaam/Implementation-of-the-Gaussian-Regression-SemiPar-Dataset/blob/7a6987f4382f4ae789c9eebdf23caea98dfb3e08/Basic.R) : As the main goal of the project is to work on the optimization of the gaussian regression, we started by developing a first basic code, that will later get enriched.
→ The “basic algorithm” took random parameters, and the covariance function was also randomly chosen ( Matern 3/2).

- [Maten32-80.R](https://github.com/Azkaaaaaam/Implementation-of-the-Gaussian-Regression-SemiPar-Dataset/blob/7a6987f4382f4ae789c9eebdf23caea98dfb3e08/Maten32-80.R) : Given the first results obtained from the basic code, our aim is to work on its optimization. Since we are dealing with a supervised learning technique, one of the first optimizations we may do is dividing the dataset into test and training sets. 

### Implemented Algorithms:
Now that we have divided our dataset, we are going to focus on the optimization of the covariance functions since they are the main determinants of the regression output. We opted for 4 main functions: Matern 3/2 - Matern 5/2 - Exponential and Squared Exponential

- [Maten32.R](https://github.com/Azkaaaaaam/Implementation-of-the-Gaussian-Regression-SemiPar-Dataset/blob/7a6987f4382f4ae789c9eebdf23caea98dfb3e08/Maten32.R)

- [Matern 5-2.R](https://github.com/Azkaaaaaam/Implementation-of-the-Gaussian-Regression-SemiPar-Dataset/blob/7a6987f4382f4ae789c9eebdf23caea98dfb3e08/Matern%205-2.R)

- [Squared_Expo.R](https://github.com/Azkaaaaaam/Implementation-of-the-Gaussian-Regression-SemiPar-Dataset/blob/7a6987f4382f4ae789c9eebdf23caea98dfb3e08/Squared_Expo.R)

### Parameters Optimization & Final Optimized Code:
It should be known that our optimum function depends enormously on the initialized value. This was deducted from running the function with different initialized values manually when all the optimized parameters were close to the initial values . 

- [Optimization](https://github.com/Azkaaaaaam/Implementation-of-the-Gaussian-Regression-SemiPar-Dataset/blob/7a6987f4382f4ae789c9eebdf23caea98dfb3e08/Optimization.R)

- [Most Optimized Algo.R](https://github.com/Azkaaaaaam/Implementation-of-the-Gaussian-Regression-SemiPar-Dataset/blob/7a6987f4382f4ae789c9eebdf23caea98dfb3e08/Most%20Optimized%20Algo.R)

# General Overview
This project dealt with the implementation of a gaussian process to the data set of canadian log income in regards to age. In order to get the best results, after inspecting the data we split our data into two sets; the training and the testing to test our model. And then we started with the optimization of the model by changing the covariance function. The mean of squared predicted error was used to compare how the model is performing each time. At a first instance we got the optimal MSPE from the covariance function of squared exponential. After that we needed to do further optimization on the hyperparameters using the optimization function in R with several methods such as SANN , CG and BFGS. This function  was used on all cases of the covariance to compare if after the optimization will get different results. Even Though we got closer values of MSPE, the squared exponential gave us the best results and optimal model with the use of CG as an optimization method. 












______________________________________________________________________________________________________________________________________________________________________________

In collaboration with Mouna Hmani
