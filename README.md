---
title: "Power Analysis in R for Multilevel Models"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
I took the easy way out and reassigned correlations to the correct places.

Cannot figure out how to include the standard deviation.  I think the correlation matrix is serving as the 

```{r}
set.seed(123)
corrData = matrix(runif(25, min = .5, max = .7), ncol = 5); testCorr
corrData[1,2] = corrData[2,1]; corrData[1,3] = corrData[3,1]; corrData[1,4] = corrData[4,1]; corrData[1,5] = corrData[1,5] 
corrData[2,3] = corrData[3,2]; corrData[2,4] = corrData[4,2]; corrData[2,5] = corrData[5,2]
corrData[3,4] = corrData[4,3]; corrData[3,5] = corrData[5,3]; corrData[4,5] = corrData[5,4]

corrData
n = 20000

origEig <- eigen(corrData)
origEig
#library(psych)
#corrData = cor.smooth(corrData)
#origEig <- eigen(corrData)
#origEig

library(MASS)
set.seed(123)
ITP <- mvrnorm(n, mu = c(rep(3.5,5)), Sigma = corrData, empirical = TRUE)
cor(ITP)
cov(ITP)
cov(corrData)
corrData
TIME = rep(1:5, 4000)
errorITP = ITP - 2 - .5*(TIME)
ITP = ITP + errorITP
colnames(ITP) = c("time1", "time2", "time3", "time4", "time5")
ITP  = as.data.frame(ITP)
head(ITP)
ITP = apply(ITP, 2, function(x){ifelse(x >= 7,7,ifelse(x <= 1, 1, x))})
apply(ITP, 2, mean)
apply(ITP, 2, sd)
ITP = as.data.frame(ITP)
head(ITP)
library(reshape)
dataSTEM = reshape(ITP, varying = list(c("time1", "time2", "time3", "time4", "time5")), times = c(1,2,3,4,5), direction = "long")
colnames(dataSTEM) = c("TIME", "ITP", "ID")
dataSTEM
head(ITP)

WHITE = c(rep(1,.7*20000), rep(0, .3*20000))
length(WHITE)

dataTest = cbind(dataSTEM, WHITE)
dataTest = as.data.frame(dataTest)
head(dataTest)
```
Group Level predictors for lme: http://www.rensenieuwenhuis.nl/r-sessions-16-multilevel-model-specification-lme4/

It seems like for R, level two and level one are placed in the "fixed" section if they not varying slopes.  For example, if I replaced 1 with time, then I would be evaluating how each time point varying for each person.
```{r}
library(lme4)

model = lmer(ITP ~  WHITE + TIME + (1| ID), data = dataTest)
summary(model)
```
Here is how you set the parameter estimates to what you want them to be to estimate possible power. Here I have .6 which is a little more than half of one of 1 option so maybe that is reasonable?

This software only does one parameter at a time, but if it works for AI, which is the lowest number of people, then it should be good for the other ethniticies, because those have more people.

Question: Is there anyway to figure out what a .4 effect size would be on this scale?  I tried transforming, but I don't think I did it right.
```{r}
#fixef(model)["AA"] <- 0.4
#fixef(model)["HIS"] <-0.4
#fixef(model)["BL"] <- 0.4
fixef(model)["AI"] <- 0.6

powerSim(model, nsim = 10)

dataModel1 = getData(model)
write.csv(dataModel1, "dataModel1.csv")
# Added 10 white people every 10 n's is 50 data points.
model2 = extend(model, along = "id", n = 140)


fixef(model2)["AI"] <- 0.6
powerSim(model2, nsim = 100)




```
