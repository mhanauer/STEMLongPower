---
title: "Power Analysis in R for Multilevel Models"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
Here I am creating the correlation matrices.  I am assuming that we need to include correlations for time points as well as ethnicities.  Although, a correlation matrix isn't best the measurement tool for a Pearson's R correlation isn't the best representation of a correlation between a binary variable (ethnicity) and a countinous variable (intent to persist). 

First step is generate correlations for each of the four ethnicitices.  Using a range of .2 to .4, but will need to get the some literature to back this up. 
```{r}
#Generate 8*8 runif and then insert them one by one and repeat them as necessary.
set.seed(123)
corrData = runif(64/2, max = .4, min = .2)
write.csv(corrData, "corrData.csv")
corrData = as.matrix(read.csv("corrData.csv", header = FALSE))
```
Now we have the correlation matrix and either need to use that directly or create a data from that.

I created a multivariate distribution based upon the correlation matrix that I created.  I assumed the means were zero for the sake of simplicity in this example, but may need to change them based upon research.  I then took the first five variables, because those are the variables for the first time points which I believe is what we want for the dependent variable.  
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
This is just me creating an artificial data set to demonstrate how the model works.  It has a total of 70+74 = 148 participants with the following breakdowns for ethnicity:  20 Hispanic, 20 African American, 10 Native American, 20 Asian American, and 74 White.  Also, for this example, I am assuming only 4 time points, which is what we want for the quantitative portion.  If we can get enough power for four time points, then we will have enough power for 10 time points, which is what we have for the qualitative data. 

Here I am creating the data set (need to make it for five times not four)

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
Here are some examples of what I think we want.  I believe our goal is to get estimates for each of the four ethnicities at each of the four time points.  So I put together a model that has each of the four ethnicities in the fixed effects, which gives us the average estimate of each ethnicity on the dependent variable over time, while the random estimates provide the effect each ethnicity has over each of the four time points.  

For example, ranef(model) extracts the random effects.  For HIS (i.e. Hispanic), there is a random intercept (the column below HIS) and a random slope (the HIS column) for each time point.  Therefore, we can evaluate Hispanics levels of a dependent variable at each time point and if the trajectory (i.e. the slope) for hispanics is different for each time point relative to whites (white is the reference category that is left out). 

All of these estimates are relative to white individuals, because they are left out of the model.  However, if we wanted to understand how different ethnicities vary with each other we can develop different models with different reference categories.  For example we could include white and drop black, then all the parameter estimates would be in reference to black students. With black as the reference category, the parameter estimate for Hispanics would be the difference in the dependent variable relative to blacks at different times points. 

My understanding is that in this model, level one is the ethnicity and level two is time.  I don’t think modeling person within time makes sense for several reasons.  First, I do not believe we are interested in the estimates for individuals only for ethnic groups.  Second, modeling person people over time will make the model more complex likely increasing the number of observations that we need.  Third, when I tried running person within time, the model failed (wouldn’t run), so with this set up, it may not be possible to model person within time.
```{r}
library(lme4)
model = lmer(y~ AI + HIS + BL + AA + (1 | id), data = dataTest)
summary(model)
```

```{r}
fixef(model)["AA"] <- 0.6
fixef(model)["HIS"] <-0.6
fixef(model)["BL"] <- 0.6
fixef(model)["AI"] <- 0.6
fixef(model)["AI"] <- 0.6

powerSim(model, nsim = 10)
```
