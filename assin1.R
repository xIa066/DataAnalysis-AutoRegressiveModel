library(ISLR)
library(MASS)
library(corrplot)
library(tidyverse)
library(dplyr)
library(glmnet)
library(stats)
library(ggplot2)
require(reshape2)
require(directlabels)

# Q1.1
?Boston
data("Boston")
plot(Boston)
summary(Boston)
#Q1.2
y.fit <- lm(medv~rm, data= Boston) 
plot(Boston$rm, Boston$medv)
abline(y.fit, col="red")

#Q1.3
summary(y.fit)

cor(Boston) %>%
  corrplot()

lm(medv~lstat, data= Boston) %>%
  summary()

#Q1.4
y.mfit <- lm(medv~., data=Boston)
summary(y.mfit)

#Q1.5
plot(Boston)

coef(y.mfit)
#Q2

#Q2.1
data(Hitters)
# this function takes a df and a column name
# return a df that filtered out NA value of that column
f1 <- function(data, desiredCols) {
  print(desiredCols)
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec,])
}

cleanedHitter <-f1(Hitters,"Salary")
dim(cleanedHitter)

#Q2.2
x <- model.matrix(~.,data = subset(cleanedHitter, select = -c(Salary)))
dim(x)
# trim first column for fitting
x<- x[,-1]
y <- cleanedHitter %>% select(Salary) %>% unlist() %>% as.numeric()
lambdaSeq <- sort(10^seq(10,-2,length=100), decreasing=TRUE)
modQ20<- glmnet(x,y,family = "gaussian",lambda = lambdaSeq, alpha = 0)


modQ20$beta[,1]
modQ20$beta[,100]

#Q2.3

betaL2Norm <- modQ20$beta[1:19,] %>%
  apply(. , MARGIN=2, FUN=function(x) sqrt(sum(x^2)))

df <- data.frame(betaL2Norm=betaL2Norm, lambda=log(lambdaSeq))
ggplot(df, aes(lambda, betaL2Norm)) + geom_line() + geom_point()

# no, because l2 norm of the parameters does not reflect the value the residuals
# of the model fit. Namely, the ridge regression minimize the sum of residuals 
# and the l2 norm of estimator times lambda. 

#Q2.4
set.seed(10)
sample <- sample.int(n = nrow(cleanedHitter), size = 131, replace = F)
Train <- cleanedHitter[sample,]
Test <- cleanedHitter[-sample,]

xTrain <- model.matrix(~.,data = subset(Train, select = -c(Salary)))[,-1]
yTrain <- Train$Salary
xTest <- model.matrix(~.,data = subset(Test, select = -c(Salary)))[,-1]
yTest <- Test$Salary

cv_fitRidge <- cv.glmnet(xTrain,
                    yTrain,
                    alpha=0,
                    lambda = lambdaSeq,
                    type.measure = "mse") 

plot(cv_fitRidge)
opt_lambda <- cv_fitRidge$lambda.min
opt_mse <- min(cv_fitRidge$cvm)

fit <- cv_fitRidge$glmnet.fit
y_predictedRidge <- predict(fit, s = opt_lambda, newx = xTest)
mean((yTest - y_predictedRidge)^2)

modQ2_4 <- glmnet(x,y,family = "gaussian", alpha = 0)

predict(modQ2_4, type = "coefficients", s = opt_lambda) # Display coefficients using lambda chosen by CV

#Q2.4

cv_fitLasso <- cv.glmnet(xTrain,
                    yTrain,
                    alpha=1,
                    lambda = lambdaSeq,
                    type.measure = "mse") 

plot(cv_fitLasso)
opt_lambda1 <- cv_fitLasso$lambda.min
opt_mse1 <- min(cv_fitLasso$cvm)

y_predictedLasso <- 
  cv_fitLasso$glmnet.fit %>%
  predict(s = opt_lambda1, newx = xTest) 

mean((yTest - y_predictedLasso)^2)
modQ2_4 <- glmnet(x,y,family = "gaussian", alpha = 1)
lasso_coef <- predict(modQ2_4, type = "coefficients", s = opt_lambda1) # Display coefficients using lambda chosen by C
lasso_coef[lasso_coef!=0]

# Q3.4

# Generate eps, a vector of 1100 draws from the Gaussian distribution;
# Initialize the first 3 elements in the vector x with zeros: x[1] = x[2] = x[3] = 0.
# Run the loop: for i=4 to i = 1100 do x[i] = a1 * x[i-1] + a2 * x[i-2] + a3 * x[i-3]
# + eps[i] and discard the first 100 observations in x to reduce the effect of 
# initialization in step 2.


#  rep(0, 5)


M1 <- function(y_t, i){
  return(0.434*y_t[i-1]+0.217*y_t[i-2]+0.145*y_t[i-3]
         +0.108*y_t[i-4]+0.087*y_t[i-5])
}

M2 <- function(y_t, i){
  return(0.682*y_t[i-1]+0.346*y_t[i-2])
}


simQ3_4 <- function(y_t, M, p){
  y_t <-  c( rep(0,p))
  eta <- rnorm(100+p, mean = 0, sd = 1)
  for(i in (p+1):(100+p)) {
    y_t <- c(y_t, M(y_t,i) + eta[i])
  }
  return(y_t[(p+1):(100+p)])
}

y_M1 <- simQ3_4(y_M1, M1, 5)
y_M2 <- simQ3_4(y_M2, M2, 2)

plot(y_M1)
plot(y_M2)

# Q3.5
library(ggplot2)
require(reshape2)
require(directlabels)

YDesign <- function(y_t,p){
  n <- length(y_t)
  y <- y_t[p:(n-1)]
  if(p ==1){
    return(as.matrix(y, nrow = length(y)))
  } else {
    for(i in 2:p){
      indexi <- p-i+1
      y <- cbind(y,y_t[indexi:(n-i)])
    } 
    return(as.matrix(y, nrow = n-p, ncol = p))    
  }
}

YResponse <- function(y_t, p){
  n <- length(y_t)
  return(as.matrix(y_t[(p+1):n], 
                   nrow= length(y_t[(p+1):n])))
}

IC1 <- function(fit){
  T <- length(fit$residuals)
  p <- length(fit$coefficients)
  RSS <- sum(resid(fit)^2)
  sigma <- RSS/T
  ic_1 <- log(sigma) + ((2*(p+1))/T)
  return(ic_1)
}

IC2 <- function(fit){
  T <- length(fit$residuals)
  p <- length(fit$coefficients)
  RSS <- sum(resid(fit)^2)
  sigma <- RSS/T
  ic_2 <- log(sigma) + (T+p)/(T-p-2)
  return(ic_2)
}

IC3 <- function(fit){
  T <- length(fit$residuals)
  p <- length(fit$coefficients)
  RSS <- sum(resid(fit)^2)
  sigma <- RSS/T
  ic_3 <- log(sigma) + ((p*log(T))/T)
  return(ic_3)
}

IC1_M1q35 <- numeric()
IC2_M1q35 <- numeric()
IC3_M1q35 <- numeric()
for(i in 1:10){
  yDesign <- YDesign(y_M1,i)
  yResponse <- YResponse(y_M1,i)
  fit <- lm(yResponse~yDesign+0)
  IC1_M1q35 <- append(IC1_M1q35,IC1(fit))
  IC2_M1q35 <- append(IC2_M1q35,IC2(fit))
  IC3_M1q35 <- append(IC3_M1q35,IC3(fit))
}


df <- data.frame(IC1=IC1_M1q35, IC2=IC2_M1q35, IC3=IC3_M1q35, p=1:10 )
df.m <- melt(df,id.vars="p")
p <- ggplot(df.m, aes(x=p, y=value, color=variable)) + geom_line() +geom_point()
direct.label(p)

# Q3.6
M1 <- function(y_t, i){
  return(0.434*y_t[i-1]+0.217*y_t[i-2]+0.145*y_t[i-3]+0.108*y_t[i-4]+0.087*y_t[i-5])
}

YDesign <- function(y_t,p){
  n <- length(y_t)
  y <- y_t[p:(n-1)]
  if(p ==1){
    return(as.matrix(y, nrow = length(y)))
  } else {
    for(i in 2:p){
      indexi <- p-i+1
      y <- cbind(y,y_t[indexi:(n-i)])
    } 
    return(as.matrix(y, nrow = n-p, ncol = p))    
  }
}

YResponse <- function(y_t, p){
  n <- length(y_t)
  return(as.matrix(y_t[(p+1):n], 
                   nrow= length(y_t[(p+1):n])))
}


IC1 <- function(fit){
  T <- length(fit$residuals)
  p <- length(fit$coefficients)
  RSS <- sum(resid(fit)^2)
  sigma <- RSS/T
  ic_1 <- log(sigma) + ((2*(p+1))/T)
  return(ic_1)
}

IC2 <- function(fit){
  T <- length(fit$residuals)
  p <- length(fit$coefficients)
  RSS <- sum(resid(fit)^2)
  sigma <- RSS/T
  ic_2 <- log(sigma) + (T+p)/(T-p-2)
  return(ic_2)
}

IC3 <- function(fit){
  T <- length(fit$residuals)
  p <- length(fit$coefficients)
  RSS <- sum(resid(fit)^2)
  sigma <- RSS/T
  ic_3 <- log(sigma) + ((p*log(T))/T)
  return(ic_3)
}


simQ3_6 <- function(M, p, n){
  y_t <- numeric()
  y_t <-  c( rep(0,p))
  eta <- rnorm(n+p, mean = 0, sd = 1)
  for(i in (p+1):(n+p)) {
    y_t <- c(y_t, M(y_t,i) + eta[i])
  }
  return(y_t[(p+1):(n+p)])
}

result3.6 <- c()
# fixed IC1, 100 sample, underlying process M1
for(ic in c(IC1, IC2, IC3)){
  dfQ3_6 <- data.frame(p1 = 0, p2 = 0, p3 = 0, p4 = 0, p5 = 0, p6 =0, p7 = 0,
                       p8 = 0, p9 = 0, p10 =0)
  for(i in 1:1000){
    sampleSize <- 100
    # simulate 100 M1 sample
    yQ3_6 <- simQ3_6(M1, 5,sampleSize)
    IC <- numeric()
    for(p in 1:10){
      yDesgin <-YDesign(yQ3_6,p)
      yResponse <-YResponse(yQ3_6,p)
      fit <- lm(yResponse~yDesgin+0)
      IC <- append(IC, ic(fit))
    }
    dfQ3_6[1,which(IC == min(IC))] <- dfQ3_6[1,which(IC == min(IC))] + 1
  }
  result3.6<-c(result3.6, dfQ3_6)
}

# IC1 counts
result3.6[1:10]
# IC2 counts
result3.6[11:20]
# IC3 counts
result3.6[21:30]

print(df3.6)
# Q3.7

result3.7 <- c()
# fixed IC1, 100 sample, underlying process M1
for(ic in c(IC1, IC2, IC3)){
  df <- data.frame(p1 = 0, p2 = 0, p3 = 0, p4 = 0, p5 = 0, p6 =0, p7 = 0,
                       p8 = 0, p9 = 0, p10 =0)
  for(i in 1:1000){
    sampleSize <- 15
    # simulate 100 M1 sample
    yQ3_6 <- simQ3_6(M1, 5,sampleSize)
    IC <- numeric()
    for(p in 1:10){
      yDesgin <-YDesign(yQ3_6,p)
      yResponse <-YResponse(yQ3_6,p)
      fit <- lm(yResponse~yDesgin+0)
      IC <- append(IC, ic(fit))
    }
    print(min(IC))
    df[1,which(IC == min(IC))] <- df[1,which(IC == min(IC))] + 1
  }
  result3.7<-c(result3.7, df)
}

?Boston
plot(Boston$black, Boston$medv)
summary(Boston)


# n = 25, 100

calculateProb <- function(IC, n, p){
  result <- numeric()
  for(L in 1:8){
    if(IC == "IC1"){
      criticalValue <- ((n-p)/L) * exp((2*p+2+2*L)/(n-p-L)- (2*p+2)/(n-p)) - ((n-p-L)/(L))
    }else if(IC == "IC2"){
      criticalValue <- ((n-p)/L) * exp((n)/(n-2*p-2*L-2)- (n)/(n-2*p-2)) - ((n-p-L)/(L))
    }else if(IC == "IC3"){
      criticalValue <- ((n-p)/L) * exp((p+L)*log(n-p-L)/(n-p-L)- (p*log(n-p))/(n-p)) - 
        ((n-p-L)/(L))
    }else{
      print("invalid input IC")
    }
    result <- c(result, pf(criticalValue, L, n-p-L))
  }
}
