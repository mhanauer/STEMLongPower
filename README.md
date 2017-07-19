---
title: "Power Analysis in R for Multilevel Models"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
Here is an example using an artificial data set as pilot data to estimate power for a random intercepts model.  The z variable is a count dependent variable, while x is a time variable going from 1 to 10 (i.e. there are ten measurements per person), and g is the group factor splitting the data into three groups a, b, and c, which for this example we will assume are different people.  In this simulation study, we want to estimate the sample size required for detecting an effect for time (i.e. the x variable) of -.045.  That is, we want to know the number of participants (i.e. the number of groups in variable g) to estimate an average change in the dependent variable, z, over time to have a power of .8 that would be able to detect a slope of -.045 (i.e. 80% chance of not making a type two error).  

We begin by setting up the model, which is a standard multilevel model with g as the hieratical variable with three levels a, b, c (i.e. three different people) getting their own intercepts (i.e. random intercepts).  Then we change x's slope to -.045 and run the powerCurve model to estimate the number of people needed to achieve a desired level of power (usually .8).  The powerCurve function usually runs 1,000 simulations, which in real data analysis the user should use; however, because 1000 simulations can take a long time, I have set the number of simulations (nsim) to 10.  Additionally, because in this example, we are interested in understanding the number of participants (i.e. the g variable), I specify along "g".  The default is along the x variable, which, in this example, would provide the power for different amounts of measurements take across three participants needed to achieve a certain level of power.

As the reader can see the power is very low (.2) to detect a slope of -045 with three participants each receiving ten measurements.

Ok looks like everything needs some pilot data, so will need to find a similar study and get data from them to construct your outcome variable and the relationships to the independent variable.
```{r}
library(simr)

data("simdata")
head(simdata)
y = rnorm(100) 
eth = c(rep("AA", 10), rep("AI", 10), rep("H", 10), rep("W", 70))
person = rep(1:10,10)
time = rep(1:4, 25)
length(time)

dataTest = cbind(y,eth)

length(eth)

model = glmer(z~x + (1|g), family = "poisson", data = simdata)
model

fixef(model)["x"]
fixef(model)["x"] = -.045

set.seed(123)
powerCurve(model , along = "g", nsim =  10)
```
Therefore, we need to simulate increases in the sample to evaluate how many people are needed across ten measurements to achieve 80% power.  We can do this by using the extend function.  With the extend function, we set the simulation to simulate the model with participants 15 participants (i.e. 150 data points), because we specified the extension along the participant variable g with n = 15.  Then we run the powerCurve function along g again.  The results show that with ten measurements per person, we would need at least 10 people to have 80% power in our model.  
```{r}
model13 = extend(model, along = "g", n = 15)
set.seed(123)
powerCurve(model13, nsim = 10, along = "g")
```

