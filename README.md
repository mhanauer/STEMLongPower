---
title: "Power Analysis in R for Multilevel Models"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
I took the easy way out and reassigned correlations to the correct places.

The multivariate normal distribution packages that I have seen in R do not take standard deviations they take a covariance matrix.  I think the covariance matrix is serving as the sd.   I do think a covariance ranging from .5 to .7 is probably ok, but I am not sure.  Also I assumed the variances were all the same with an sd of .7 (since it is the variance I took the sqrt(.7)).
```{r}

n = 20000
set.seed(123)

r0s = rnorm(n,0,.2)
r1s = rnorm(n,0,.2)
B00 = 1 + r0s
B10 = .5 + r1s
time = rep(1:5, 4000)
eis = rnorm(n,0,.2)
# Not sure what to do with the time part
Yis = B00 + B10 + .5*(time) +eis
Yis




error = Yobs - 
origEig <- eigen(corrData)
origEig
library(psych)
corrData = cor.smooth(corrData)
origEig <- eigen(corrData)
origEig

library(MASS)
set.seed(123)
ITP <- mvrnorm(n, mu = c(rep(3.5,5)), Sigma = corrData)
TIME = rep(1:5, 4000)
errorITP = ITP - 2 - .5*(TIME)
ITP = ITP + errorITP
colnames(ITP) = c("time1", "time2", "time3", "time4", "time5")
ITP  = as.data.frame(ITP)
ITP = apply(ITP, 2, function(x){ifelse(x >= 7,7,ifelse(x <= 1, 1, x))})
ITP = as.data.frame(ITP)
library(reshape)
dataSTEM = reshape(ITP, varying = list(c("time1", "time2", "time3", "time4", "time5")), times = c(1,2,3,4,5), direction = "long")
colnames(dataSTEM) = c("TIME", "ITP", "ID")
WHITE = c(rep(1,.7*20000), rep(0, .3*20000))
dataTest = cbind(dataSTEM, WHITE)
dataTest = as.data.frame(dataTest)
head(dataTest)
```
Group Level predictors for lme: http://www.rensenieuwenhuis.nl/r-sessions-16-multilevel-model-specification-lme4/

It seems like for R, level two and level one are placed in the "fixed" section if they not varying slopes.  For example, if I replaced 1 with time, then I would be evaluating how each time point varying for each person.

It seems like the data generating process is not recovering the parameters.  
```{r}
library(lme4)

model = lmer(ITP ~  WHITE + TIME + (1| ID), data = dataTest)
summary(model)
```
