---
title: "Power Analysis in R for Multilevel Models"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
This is just me creating an artificial data set to demonstrate how the model works.  It has a total of 180 participants with the following breakdowns for ethnicity:  20 Hispanic, 20 African American, 10 Native American, 20 Asian American, and 110 White.  Also, for this example, I am assuming only 4 time points, which is what we want for the quantitative portion.  If we can get enough power for four time points, then we will have enough power for 10 time points, which is what we have for the qualitative data.   
```{r}
library(simr)
library(lme4)
set.seed(123)
y = rnorm(720)
eth = c(rep("HIS", 4*20), rep("BL", 4*20), rep("AI", 4*10), rep("AA", 4*20), rep("WHITE", 4*40+200+80))
eth = as.data.frame(eth)
HIS = ifelse(eth == "HIS", 1, 0)
BL = ifelse(eth == "BL", 1, 0)
AI = ifelse(eth == "AI", 1,0)
AA = ifelse(eth == "AA", 1, 0)
WHITE = ifelse(eth == "WHITE", 1, 0)
person = rep(1:180, each = 4)
person = as.data.frame(person)
dim(person)
time = rep(1:4, 180)
time = as.data.frame(time)
dataTest = cbind(y, HIS, BL , AI , AA , WHITE, person, time)
names(dataTest) = c("y", "HIS", "BL", "AI" , "AA" , "WHITE", "person", "time")
dataTest = as.data.frame(dataTest)
```
Here are some examples of what I think we want.  I believe our goal is to get estimates for each of the four ethnicities at each of the four time points.  So I put together a model that has each of the four ethnicities in the fixed effects, which gives us the average estimate of each ethnicity on the dependent variable over time, while the random estimates provide the effect each ethnicity has over each of the four time points.  

For example, ranef(model) extracts the random effects.  For HIS (i.e. Hispanic), there is a random intercept (the column below HIS) and a random slope (the HIS column) for each time point.  Therefore, we can evaluate Hispanics levels of a dependent variable at each time point and if the trajectory (i.e. the slope) for hispanics is different for each time point relative to whites (white is the reference category that is left out). 

All of these estimates are relative to white individuals, because they are left out of the model.  However, if we wanted to understand how different ethnicities vary with each other we can develop different models with different reference categories.  For example we could include white and drop black, then all the parameter estimates would be in reference to black students. With black as the reference category, the parameter estimate for Hispanics would be the difference in the dependent variable relative to blacks at different times points. 

My understanding is that in this model, level one is the ethnicity and level two is time.  I don’t think modeling person within time makes sense for several reasons.  First, I do not believe we are interested in the estimates for individuals only for ethnic groups.  Second, modeling person people over time will make the model more complex likely increasing the number of observations that we need.  Third, when I tried running person within time, the model failed (wouldn’t run), so with this set up, it may not be possible to model person within time.
```{r}
model = lmer(y~ AI + HIS + BL + AA + (AI | time) + (HIS  | time) + (BL | time)  + (AA | time) , data = dataTest)
ranef(model)
```
Here is an example of the power analysis package that I know how to use.  There are two problems.  First, is that this package, and all others that I have tried, do not provide power for the random effects only the fixed effects, which I do not believe we are interested in, so that may not be very useful.  However, it may be the case that we do not need to estimate power for the random effects and may only need to ensure that we have enough power for the fixed effects.  

Even if we only need to estimate the power for the fixed effects, it appears as though we can only get the power for one independent variable (i.e. one ethnicity at a time) at a time.  Therefore, we may need to run a separate power analysis for each ethnicity.  This should be fine, because power analysis that we run has the full model (i.e. all the ethnicities are included) it just estimates the power for one ethnicity at a time.  Therefore, each time we run the power analysis while estimating the power for a different ethnicity, the model is the same each time.

If we do want to estimate the power for the fixed effects only, we will need to establish what reasonable parameter estimates are for each ethnicity.  Below I have assumed that each parameter estimate is .6; however, we will need to look at the literature and assess what a reasonable parameter estimate for each of the ethnicities for each of the dependent variables that we are interested in. 
```{r}
fixef(model)["AA"] <- 0.6
fixef(model)["HIS"] <-0.6
fixef(model)["BL"] <- 0.6
fixef(model)["AI"] <- 0.6
fixef(model)["AI"] <- 0.6

powerSim(model, nsim = 10)
```



