---
title: "Power Analysis in R for Multilevel Models"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
I created each of the erroer terms and I also set B00 to 3, because when I got the min and max with an intercept of 3 the values were closer to a 1 through 7 value.
```{r}

n = 20000
set.seed(123)

r0s = rnorm(n,0,.2)
r1s = rnorm(n,0,.2)
B00 = 3
B10 = .5 
time = rep(1:5, 4000)
eis = rnorm(n,0,.2)
Yis = B00 + B10*(time) +eis + r0s + r1s
head(Yis)
min(Yis); max(Yis)
id = rep(1:4000, each = 5)
dataTest = as.data.frame(cbind(Yis, time, id))

```
I think, because we included an error term for the slope this model is a random slopes model, so I included the time as a random slope.  It seems like I am able recover the fixed intercept and slope.
```{r}
library(lme4)

model = lmer(Yis ~ + time + (time| id), data = dataTest)
summary(model)
```
