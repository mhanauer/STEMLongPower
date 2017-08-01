---
title: "Power Analysis in R for Multilevel Models"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
Here I am creating the correlation matrix.  I am assuming that we need to include correlations for all five time points as well as ethnicities, because four ethnicities (Hispanic, African American, American Indian, and Asian) will be included in the model.  However, maybe including the binary variables isn't a good idea, because a Pearson's correlation matrix is not a good description of a correlation between a binary and continuous variable.

I am assuming that all correlation between all variables is somewhere between .2 and .4 (I think somebody told me to assume that), which I selected from a uniform distribution.  

I tried to create the correlation matrix in R, but I had a really hard time, so I just downloaded it and manually created it in excel (see corrData.csv).
```{r}
#Generate 8*8 runif and then insert them one by one and repeat them as necessary.
set.seed(123)
corrData = runif(64/2, max = .4, min = .2)
write.csv(corrData, "corrData.csv")
corrData = as.matrix(read.csv("corrData.csv", header = FALSE))
```
I have a hard time working with correlation matrices and the power analyses packages in R, so I am trying to create the dependent variable from the correlation matrix.

I created a random multivariate distribution based upon the correlation matrix that I created.  I assumed the means were zero for the sake of simplicity in this example, but may need to change them based upon research.  I then took the first five variables, because those are the variables for the five time points which I believe is what we want for the outcome variable.  I placed the five variables representing the outcome variable into one long form variable.  

Because I have had difficulty sorting in R, I sorted the data by id and time in excel then reuploaded the data (see y.csv).

I think this essentially returns data are similar to z-scores, because the mean is zero and the data are assumed to be multivariate normally distributed.
```{r}
require(MASS)
corrData = as.matrix(corrData)
dim(corrData)
set.seed(123)
dataMulti <- mvrnorm(720/5, mu = c(rep(0,9)), Sigma = corrData, empirical = TRUE)
dataMulti = dataMulti[,1:5]; dataMulti
colnames(dataMulti) = c("time1", "time2", "time3", "time4", "time5"); head(dataMulti)
library(reshape)
dataMulti = as.data.frame(dataMulti)
dataMulti = reshape(dataMulti, varying = list(c("time1", "time2", "time3", "time4", "time5")), times = c(1,2,3,4,5), direction = "long")
dim(dataMulti)

```
Then I added the other four ethnicity variables.  The data set has a total of 70+74 = 144 participants with the following breakdowns for ethnicity:  20 Hispanic, 20 African American, 10 Native American, 20 Asian American, and 74 White. 
```{r}
library(simr)
library(lme4)
set.seed(123)
#y = dataMulti
#write.csv(y, "y.csv")
y = read.csv("y.csv", header = TRUE)
eth = c(rep("HIS", 5*20), rep("BL", 5*20), rep("AI", 5*10), rep("AA", 5*20), rep("WHITE", 74*5))
eth = as.data.frame(eth)
HIS = ifelse(eth == "HIS", 1, 0)
BL = ifelse(eth == "BL", 1, 0)
AI = ifelse(eth == "AI", 1,0)
AA = ifelse(eth == "AA", 1, 0)
WHITE = ifelse(eth == "WHITE", 1, 0)

dataTest = cbind(y, HIS, BL , AI , AA , WHITE)
names(dataTest) = c("time", "y", "id", "HIS", "BL", "AI" , "AA" , "WHITE")
dataTest = as.data.frame(dataTest)
```
This is my understanding of the model.  It makes sense that we want time to be random, because we do not want to assume that every time point starts at the same level (i.e. intercept).  It also makes sense, because people are nested in time points not time point nested in people.        
```{r}
library(lme4)
model = lmer(y~ AI + HIS + BL + AA + (1 | time), data = dataTest)
summary(model)
```
Given that the dependent variable is essentially, a z-score, which is comparable to a Cohen's D, I am assuming a .4 effect size for each ethnicity.  We are a little short with the current assumptions of 144 people. 
```{r}
#fixef(model)["AA"] <- 0.4
#fixef(model)["HIS"] <-0.4
#fixef(model)["BL"] <- 0.4
fixef(model)["AI"] <- 0.4

powerSim(model, nsim = 50)
```
Try to extend the model
```{r}
set.seed(123)
model = extend()
powerCurve(model, nsim =  10)
```

