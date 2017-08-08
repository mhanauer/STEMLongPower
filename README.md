---
title: "Power Analysis in R for Multilevel Models"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
Here I am creating the correlation matrix.  I am creating two sets of correlation matrices.  One representing the measures relationship with themeselves.  So there are 5 measures so there are 4 correlations with each measure so there are 20 correlations for the measures over time.  For these correlations I am assuming a relationship between .5 and .7 using a uniform distribution.  

For the rest of the correlations, I just filled them with corrs ranging from -.2 to .4


Question: Should we create a correlation matrix that includes the correlations between the outcomes variables of interest if they are being modeled seperately?  For example, I think the model will only include one dependent variable at a time with only ethnicity and time as independent variables.  So do we need to include all variables at once?

Question: Should we include the ethnicities correlations with the dependent variables as well?

Question: Is there a more efficient to create a correlation matrix in R?  I created the correlations, but then placed into the matrix in excel by hand, because I could not figure out a more efficient way.

Some of the eigenvalues were not positive so I used a smoothing package to make them positive.


```{r}
set.seed(123)
corrData = matrix(runif(20, max = .7, min = .5), ncol =5)
write.csv(corrData, "corrData.csv")

corrData = matrix(runif(280, max = .4, min = -.2), ncol =25)
write.csv(corrData, "corrData.csv")

setwd("~/Google Drive/Kerrie/LongStem/Data")
corrData = as.matrix(read.csv("STEMLongCorr.csv", header = FALSE))
head(corrData)
origEig <- eigen(corrData)
origEig
library(psych)
corrData = cor.smooth(corrData)
origEig <- eigen(corrData)
origEig
```
Let's just do an example with only time correlations for one variable.  So we will have about 5*5 -5 correlations.
```{r}
set.seed(123)
corrData = matrix(runif(25, max = .7, min = .5), ncol =5)
diag(corrData) = 1
write.csv(corrData, "corrDataITP.csv")
setwd("~/Google Drive/Kerrie/LongStem/Data")
corrData = as.matrix(read.csv("corrDataITP.csv", header = TRUE))
origEig <- eigen(corrData)
origEig
#library(psych)
#corrData = cor.smooth(corrData)
#origEig <- eigen(corrData)
#origEig
```



Test generating ordinal values from a correlation matrix.  Remember that marginal values are the probabilities of selecting that value not condition on anything else so need to be a probability.

Get the number of categories by using the cummaltive probability as the marginals with k-1 categories and the number of c's for the different number of variables.

Have the marginals to be .2 for each category and .1 at the ends.

Correlation is not exactly the same not sure if that is a big deal or not (all close and within the .5 .7 range we specified).  Probably not since these are estimates anyways and they are close.
```{r}
library(GenOrd)
marginal = rep(list(c(.1,.3,.5, .7, .9)),5)
corrcheck(marginal)
Sigma = matrix(corrData, 5, 5)
n = 650/5
m = ordsample(n, marginal, Sigma)
m
cor(m) # compare it with Sigma
# empirical marginal distributions
cumsum(table(m[,1]))/n
cumsum(table(m[,2]))/n # compare them with the two marginal distributions
```


I have a hard time working with correlation matrices and the power analyses packages in R, so I am trying to create the dependent variable from the correlation matrix.

I created a random multivariate distribution based upon the correlation matrix that I created.  I assumed the means were zero for the sake of simplicity in this example, with a sigma according the correlation matrix.  I then took the first five variables, because those are the variables for the five time points which I believe is what we want for the outcome variable.  I placed the five variables representing the outcome variable into one long form variable.  

Because I have had difficulty sorting in R, I sorted the data by id and time in excel then reuploaded the data (see y.csv).

I then transformed the dependent variable into a z-score, because this will make setting the regression coefficents easier later.  
```{r}
require(MASS)
corrData = as.matrix(corrData)
dim(corrData)
set.seed(123)
dataMulti <- mvrnorm(650/5, mu = c(rep(0,9)), Sigma = corrData, empirical = TRUE)
dataMulti = scale(dataMulti, center = TRUE, scale = TRUE)
dataMulti = dataMulti[,1:5]; dataMulti
colnames(dataMulti) = c("time1", "time2", "time3", "time4", "time5"); head(dataMulti)
library(reshape)
dataMulti = as.data.frame(dataMulti)
dataMulti = reshape(dataMulti, varying = list(c("time1", "time2", "time3", "time4", "time5")), times = c(1,2,3,4,5), direction = "long")
dim(dataMulti)
```
Then I added the other four ethnicity variables.  The data set has a total of 70 non-white and 60 white for a total of 130 participants with the following breakdowns for ethnicity:  20 Hispanic, 20 African American, 10 Native American, 20 Asian American, and 60 White. 
```{r}
library(simr)
library(lme4)
set.seed(123)
y = dataMulti
write.csv(y, "y.csv")
setwd("~/Google Drive/Kerrie/LongStem/Data")
y = read.csv("y.csv", header = TRUE)
eth = c(rep("HIS", 5*20), rep("BL", 5*20), rep("AI", 5*10), rep("AA", 5*20), rep("WHITE", 60*5))
eth = as.data.frame(eth)
dim(eth)
HIS = ifelse(eth == "HIS", 1, 0)
BL = ifelse(eth == "BL", 1, 0)
AI = ifelse(eth == "AI", 1,0)
AA = ifelse(eth == "AA", 1, 0)
WHITE = ifelse(eth == "WHITE", 1, 0)

dataTest = cbind(y, HIS, BL , AI , AA , WHITE)
names(dataTest) = c("time", "y", "id", "HIS", "BL", "AI" , "AA" , "WHITE")
dataTest = as.data.frame(dataTest)
```
This is my understanding of the model.  It makes sense that we want time to be random, because we want to allow time to have its own intercept instead of assuming that each time point is the same.  It also makes sense, because people are nested in time points not time point nested in people.  So I think level one is people and then level two is time.
```{r}
library(lme4)
model = lmer(y~ AI + HIS + BL + AA + (1 | id), data = dataTest)
summary(model)
```
Given that the dependent variable is a z-score, which is comparable to a Cohen's D, which is also on a standard normal scale, I am assuming a .4 effect size for AI (American Indian).  The regression coefficients are comparable to Cohen's D, because it is essentially the effect of the treatment of being a particular ethnicity.  For example, we may expect that being American Indian will result in a .4 effect size, which means we can reset the American Indian parameter estimate to .4.  The power can only be tested for one independent variable at a time; however, since American Indian is the lowest number group, if we have enough power for them and are assuming similar effects sizes across ethnicities, then we should have enough power for other ethnicities.  We are a little short with only about .4 with the current assumptions and 130 people. 
```{r}
#fixef(model)["AA"] <- 0.4
#fixef(model)["HIS"] <-0.4
#fixef(model)["BL"] <- 0.4
fixef(model)["AI"] <- 0.4

powerSim(model, nsim = 100)
```
