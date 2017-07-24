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

Notes for meeting
Problem is that there is 
```{r}
library(simr)

# This will be some type of ordinal value with and sd whatever the research finds
# Think about if we need to specify the correlation between the indepdent and dependent, because we are specificying the parameter estimates, which are the correlations.
set.seed(123)
y = rnorm(720) 
# If we want to model 80 participants with four time points then we need to have each ethnicity and stem times 4 totaling 80*4 320 data points.


eth = c(rep("HIS", 4*20), rep("BL", 4*20), rep("AI", 4*10), rep("AA", 4*20), rep("WHITE", 4*40+200+80))
eth = as.data.frame(eth)
dim(eth)
HIS = ifelse(eth == "HIS", 1, 0)
BL = ifelse(eth == "BL", 1, 0)
AI = ifelse(eth == "AI", 1,0)
AA = ifelse(eth == "AA", 1, 0)
WHITE = ifelse(eth == "WHITE", 1, 0)

# These variables need to be different from eth, because there is perfect multilinarity.  Need to make sure that person stays in the same stem field for the most part.  If everything is divisable by four you should be fine
stem = c(rep("PS", 2*20), rep("TECH", 5*20), rep("ENGIN", 6*20), rep("MATH", 10*20), rep("NON", 2*40+100+80))
PS = ifelse(stem == "PS", 1, 0)
TECH = ifelse(stem == "TECH", 1, 0)
EGNIN = ifelse(stem == "EGNIN", 1, 0)
MATH = ifelse(stem == "MATH", 1,0)
NON = ifelse(stem == "NON", 1, 0)
dim(eth)
stem = as.data.frame(stem)
dim(stem)
# We want to measure the quantiative questionnaire, because if there is enough power for the quantiative questionnaire with four time points 

# Need to make sure that the person matches their ethnicity.  For person, need to make sure the total above is divisible by 16 and by 4 for time.   

person = rep(1:4, each = 4, 45)
person = as.data.frame(person)
dim(person)

time = rep(1:4, 180)
time = as.data.frame(time)
dim(time)
dataTest = cbind(y, HIS, BL , AI , AA , WHITE , PS, TECH, EGNIN, MATH, NON, person, time)
dim(dataTest)
dataTest
names(dataTest) = c("y", "HIS", "BL", "AI" , "AA" , "WHITE" , "PS", "TECH", "EGNIN", "MATH", "NON", "person", "time")
dataTest = as.data.frame(dataTest)


# Model with no fixed effects and all random.  No EGNIN still work on that
model = lmer(y~ 1 + (time | HIS) + (time | BL) + (time | AI) + (time | AA) + (time | PS) + (time | TECH) + (time |MATH) , data = dataTest)

model
confint(model)

# Something wrong with EGNIN variable to similar to another one I think.
# + (time | EGNIN)
ranef(model)

VarCorr(model)["AI"] = .5
model

set.seed(123)
powerCurve(model , along = "person", nsim =  10)
```
Therefore, we need to simulate increases in the sample to evaluate how many people are needed across ten measurements to achieve 80% power.  We can do this by using the extend function.  With the extend function, we set the simulation to simulate the model with participants 15 participants (i.e. 150 data points), because we specified the extension along the participant variable g with n = 15.  Then we run the powerCurve function along g again.  The results show that with ten measurements per person, we would need at least 10 people to have 80% power in our model.  
```{r}
model13 = extend(model, along = "g", n = 15)
set.seed(123)
powerCurve(model13, nsim = 10, along = "g")
```
Model with all fixed effects and more data
```{r}
# Model that includes all the fixed effects
model = lmer(y~ HIS + BL + AI + AA + PS + TECH + MATH + (time | HIS) + (time | BL) + (time | AI) + (time | AA) + (time | PS) + (time | TECH) + (time |MATH) +(time | person) , data = dataTest)
```


Test
```{r}
install.packages("pamm")
library(pamm)


(fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy))
ours2 <- EAMM(numsim=10, mer.model=list(model=fm1,env="Days",random="Subject"),
VI=seq(0.3,0.5,0.1), VS=seq(0.05,0.2,0.05) )
plot(ours2, "both")
```

