---
title: "Power Analysis in R for Multilevel Models"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
Here I am creating a correlation for the ITP variable with correlations ranging from .5 to .7 with each other.

I cannot figure out an easy to take the random values and place them into a correlation matrix so I did it in excel and then reuploaded it.


Question: Should we create a correlation matrix that includes the correlations between all the variables that will be included in the model?  For example should we include the correlation between the different ethnicities since that could affect the value of the y?  Or can we set the relationship that we want with the dependent variable to have with the independent variables later when we set the regression coefficients?


Question: Is there a more efficient to create a correlation matrix in R?  I created the correlations, but then placed into the matrix in excel by hand, because I could not figure out a more efficient way.

The eigen values were all positive, but there is some code to smooth if need be.


```{r}
set.seed(123)
corrData = matrix(runif(20, max = .7, min = .5), ncol =5)
write.csv(corrData, "corrData.csv")

setwd("~/Google Drive/Kerrie/LongStem/Data")
corrData = as.matrix(read.csv("corrDataITP.csv", header = TRUE))
head(corrData)
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

Not sure how to add an error term when producing ordinal values.

Here I am using the GenOrd package to produce samples of ordinal values based upon a multiple variables.  I have to specifcy a correlation matrix and cummaltive probabilites (called marginal below) for each of the categories within each of the variables.  I made the the first and last category less probability (only 10% each).  Then I used the corr matrix from before to create the data. 

It does not exactly reproduce the corr and marginals, but it is pretty close.

Question: I think if I add an error term, then the values will no longer by ordinal?
```{r}
library(GenOrd)
marginal = rep(list(c(.1,.3,.5, .7, .9)),5)
corrcheck(marginal)
Sigma = matrix(corrData, 5, 5)
n = 650/5
set.seed(123)
dataSTEM = ordsample(n, marginal, Sigma); head(dataSTEM)
cor(dataSTEM) # compare it with Sigma
# empirical marginal distributions
cumsum(table(dataSTEM[,1]))/n
cumsum(table(dataSTEM[,2]))/n # compare them with the two marginal distributions
```
Here I am just creating a long form of the data.
```{r}
colnames(dataSTEM) = c("time1", "time2", "time3", "time4", "time5"); head(dataSTEM)
library(reshape)
dataSTEM = as.data.frame(dataSTEM); head(dataSTEM)
dataSTEM = reshape(dataSTEM, varying = list(c("time1", "time2", "time3", "time4", "time5")), times = c(1,2,3,4,5), direction = "long")
colnames(dataSTEM) = c("time", "ITP", "id")
head(dataSTEM)
dim(dataSTEM)
library(dplyr)
arrange(dataSTEM, time, id)
dataSTEM
write.csv(dataSTEM, "dataSTEM.csv")
```
Now I need to add the ethnicities and assign them to the scores, which I am assuming is no different for any of them we are assuming of them are the same, because we did not differeniate between them.

Then I added the other four ethnicity variables.  The data set has a total of 70 non-white and 60 white for a total of 130 participants with the following breakdowns for ethnicity:  20 Hispanic, 20 African American, 10 Native American, 20 Asian American, and 60 White. 
```{r}
library(simr)
library(lme4)
eth = c(rep("WHITE", 60), rep("BL", 20), rep("AI", 10), rep("AA", 20), rep("HIS", 20))
eth = matrix(rep(eth,5),ncol = 5)
head(eth)
eth = as.data.frame(eth); head(eth)
dim(eth)

eth = reshape(eth, varying = list(c("V1", "V2", "V3", "V4", "V5")), times = c(1,2,3,4,5), direction = "long")
colnames(dataSTEM) = c("time", "eth", "id")
eth = eth[c(2)]
head(eth)
HIS = ifelse(eth == "HIS", 1, 0)
BL = ifelse(eth == "BL", 1, 0)
AI = ifelse(eth == "AI", 1,0)
AA = ifelse(eth == "AA", 1, 0)
WHITE = ifelse(eth == "WHITE", 1, 0)

dataTest = cbind(dataSTEM, WHITE,HIS, BL , AI , AA )
dataTest
names(dataTest) = c("time", "ITP", "id","WHITE", "HIS", "BL", "AI" , "AA" )
dataTest = as.data.frame(dataTest); head(dataTest)
dataTest
write.csv(dataTest, "dataTest.csv")
```
This is my understanding of the model. 

If we need to model as ordinal, then could use Bayesian approach: https://kevinstadler.github.io/blog/bayesian-ordinal-regression-with-random-effects-using-brms/

Group Level predictors for lme: http://www.rensenieuwenhuis.nl/r-sessions-16-multilevel-model-specification-lme4/

It seems like for R, level two and level one are placed in the "fixed" section if they not varying slopes.  For example, if I replaced 1 with time, then I would be evaluating how each time point varying for each person.
```{r}
library(lme4)

model = lmer(ITP~  AI + time + HIS + BL + AA + (1| id), data = dataTest)
summary(model)
```
Here is how you set the parameter estimates to what you want them to be to estimate possible power. Here I have .5 which is half of one of 1 option so maybe that is reasonable?

This software only does one parameter at a time, but if it works for AI, which is the lowest number of people, then it should be good for the other ethniticies, because those have more people.

Question: Is there anyway to figure out what a .4 effect size would be on this scale?  I tried transforming, but I don't think I did it right.
```{r}
#fixef(model)["AA"] <- 0.4
#fixef(model)["HIS"] <-0.4
#fixef(model)["BL"] <- 0.4
fixef(model)["AI"] <- 0.5

powerSim(model, nsim = 10)

dataModel1 = getData(model)
write.csv(dataModel1, "dataModel1.csv")
# Added 50 white people every 10 n's is 50 people.
model2 = extend(model, along = "id", n = 140)
nrow(getData(model2))
dadtaModel2 =getData(model2)
write.csv(dataModel2, "dataModel2.csv")

powerSim(model2, nsim = 100)

pc1 = powerCurve(model2, along ="id", nsim = 10)
pc1
plot(pc1)



```
