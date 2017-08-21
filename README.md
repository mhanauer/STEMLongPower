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
set.seed(123)
corrData = matrix(runif(25, min = .5, max = .7), ncol = 5)
corrData[1,2] = corrData[2,1]; corrData[1,3] = corrData[3,1]; corrData[1,4] = corrData[4,1]; corrData[1,5] = corrData[1,5] 
corrData[2,3] = corrData[3,2]; corrData[2,4] = corrData[4,2]; corrData[2,5] = corrData[5,2]
corrData[3,4] = corrData[4,3]; corrData[3,5] = corrData[5,3]; corrData[4,5] = corrData[5,4]
diag(corrData) = sqrt(.7)
corrData
n = 20000

origEig <- eigen(corrData)
origEig
library(psych)
#corrData = cor.smooth(corrData)
#origEig <- eigen(corrData)
#origEig

library(MASS)
set.seed(123)
ITP <- mvrnorm(n, mu = c(rep(3.5,5)), Sigma = corrData)
TIME = rep(1:5, 4000)
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


It seems like the data generating process is recovering the intecept which I set at 3.5  
```{r}
library(lme4)

model = lmer(ITP ~  WHITE + TIME + (TIME| ID), data = dataTest)
summary(model)
```
With the simr package, I can set the parameter estimates.  
```{r}
library(simr)
fixef(model)["TIME"] <- 0.6
fixef(model)["WHITE"] = .5

powerSim(model, nsim = 10)

dataModel1 = getData(model)
write.csv(dataModel1, "dataModel1.csv")
# Added 10 white people every 10 n's is 50 data points.
model2 = extend(model, along = "id", n = 140)


fixef(model2)["AI"] <- 0.6
powerSim(model2, nsim = 100)

```

