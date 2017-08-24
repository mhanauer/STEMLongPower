library(MASS)
library(psych)
library(reshape)
library(lme4)

###############################################################
# L1 Equation
# ITPis<-pi0s + pi1s*TIMEis + eis
# L2 Equation a
# p0s<-b00 + b01*Black+ b02*Hisp  + r0s
# L2 Equation b
# p0s<-b10 + b11*Black+ b12*Hisp  + r1s
###########################################################################
# Mixed Equation
# ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
###########################################################################

# Let's use the mixed equation
# Unknowns, but we can set the value for.
# ITPis
# b00
# b01;b02
# b10
# b11;b12
# ETHs
# r0s
# r1s
# TIMEis
##########################################################################################
### Since we have only one equation in the mixed equation, we can only have one unknowns
# Let's make the eis is the unknown
######150###################################################################################
STEMLong <- function(){
  ### High Level Parameters
  n<-150				# Sample Size
  timePoints<-5		# Number of Time Points
  
  
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.22 ; b12Val <- -.22  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(10, STEMLong())
repNum = 10
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); t_valuesT
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))
n150 = 150
# Then divide by the total to get the percentage of model that had t-values above 2, which is the power.    
Power150 = as.data.frame(t(t_valuesPower/repNum))
PowerB150 = cbind(Power150$time_black, n150)
colnames(PowerB150) = c("Power", "N")

PowerH150 = cbind(Power150$time_hisp, n150)
colnames(PowerH150) = c("Power", "N")

##############160#####################
STEMLong <- function(){
  ### High Level Parameters
  n<-160				# Sample Size
  timePoints<-5		# Number of Time Points
  
  # These are for -.2
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.22 ; b12Val <- -.22  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(1000, STEMLong())
repNum = 1000
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); t_valuesT
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))
n160 = 160
Power160 = as.data.frame(t(t_valuesPower/repNum))
PowerB160 = cbind(Power160$time_black, 160)
colnames(PowerB160) = c("Power", "N")

PowerH160 = cbind(Power160$time_hisp, 160)
colnames(PowerH160) = c("Power", "N")

#############170#######################################

STEMLong <- function(){
  ### High Level Parameters
  n<-170				# Sample Size
  timePoints<-5		# Number of Time Points
  
  # These are for -.2
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.22 ; b12Val <- -.22  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(1000, STEMLong())
repNum = 1000
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); dim(t_valuesT)
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))
n170 = 170
Power170 = as.data.frame(t(t_valuesPower/repNum))
PowerB170 = cbind(Power170$time_black, n170)
colnames(PowerB170) = c("Power", "N")

PowerH170 = cbind(Power170$time_hisp, n170)
colnames(PowerH170) = c("Power", "N")

#########180###########################################
STEMLong <- function(){
  ### High Level Parameters
  n<-180				# Sample Size
  timePoints<-5		# Number of Time Points
  
  # These are for -.2
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.22 ; b12Val <- -.22  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(1000, STEMLong())
repNum = 1000
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); t_valuesT
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))
n180 = 180
Power180 = as.data.frame(t(t_valuesPower/repNum))
PowerB180 = cbind(Power180$time_black, n180)
colnames(PowerB180) = c("Power", "N")

PowerH180 = cbind(Power180$time_hisp, n180)
colnames(PowerH180) = c("Power", "N")

####################################################
STEMLong <- function(){
  ### High Level Parameters
  n<-190				# Sample Size
  timePoints<-5		# Number of Time Points
  
  # These are for -.2
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.22 ; b12Val <- -.22  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(1000, STEMLong())
repNum = 1000
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); dim(t_valuesT)
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))

n190 = 190
Power190 = as.data.frame(t(t_valuesPower/repNum))
PowerB190 = cbind(Power190$time_black, n190)
colnames(PowerB190) = c("Power", "N")

PowerH190 = cbind(Power190$time_hisp, n190)
colnames(PowerH190) = c("Power", "N")

##########200##########################################
STEMLong <- function(){
  ### High Level Parameters
  n<-200				# Sample Size
  timePoints<-5		# Number of Time Points
  
  # These are for -.2
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.22 ; b12Val <- -.22  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(1000, STEMLong())
repNum = 1000
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); dim(t_valuesT)
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))

n200 = 200
Power200 = as.data.frame(t(t_valuesPower/repNum))
PowerB200 = cbind(Power200$time_black, n200)
colnames(PowerB200) = c("Power", "N")

PowerH200 = cbind(Power200$time_hisp, n200)
colnames(PowerH200) = c("Power", "N")

#######Power Table#############################################
# Just rbind everything and plot it 
PowerTableH_22 = rbind(PowerH150, PowerH160, PowerH170, PowerH180, PowerH190, PowerH200)
PowerTableB_22 = rbind(PowerB150, PowerB160, PowerB170, PowerB180, PowerB190, PowerB200)
PowerTableH_22
PowerTableB_22

library(MASS)
library(psych)
library(reshape)
library(lme4)

###############################################################
# L1 Equation
# ITPis<-pi0s + pi1s*TIMEis + eis
# L2 Equation a
# p0s<-b00 + b01*Black+ b02*Hisp  + r0s
# L2 Equation b
# p0s<-b10 + b11*Black+ b12*Hisp  + r1s
###########################################################################
# Mixed Equation
# ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
###########################################################################

# Let's use the mixed equation
# Unknowns, but we can set the value for.
# ITPis
# b00
# b01;b02
# b10
# b11;b12
# ETHs
# r0s
# r1s
# TIMEis
##########################################################################################
### Since we have only one equation in the mixed equation, we can only have one unknowns
# Let's make the eis is the unknown
######150###################################################################################
STEMLong <- function(){
  ### High Level Parameters
  n<-150				# Sample Size
  timePoints<-5		# Number of Time Points
  
  
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.24 ; b12Val <- -.24  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(10, STEMLong())
repNum = 10
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); t_valuesT
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))
n150 = 150
# Then divide by the total to get the percentage of model that had t-values above 2, which is the power.    
Power150 = as.data.frame(t(t_valuesPower/repNum))
PowerB150 = cbind(Power150$time_black, n150)
colnames(PowerB150) = c("Power", "N")

PowerH150 = cbind(Power150$time_hisp, n150)
colnames(PowerH150) = c("Power", "N")

##############160#####################
STEMLong <- function(){
  ### High Level Parameters
  n<-160				# Sample Size
  timePoints<-5		# Number of Time Points
  
  # These are for -.2
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.24 ; b12Val <- -.24  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(1000, STEMLong())
repNum = 1000
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); t_valuesT
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))
n160 = 160
Power160 = as.data.frame(t(t_valuesPower/repNum))
PowerB160 = cbind(Power160$time_black, 160)
colnames(PowerB160) = c("Power", "N")

PowerH160 = cbind(Power160$time_hisp, 160)
colnames(PowerH160) = c("Power", "N")

#############170#######################################

STEMLong <- function(){
  ### High Level Parameters
  n<-170				# Sample Size
  timePoints<-5		# Number of Time Points
  
  # These are for -.2
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.24 ; b12Val <- -.24  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(1000, STEMLong())
repNum = 1000
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); dim(t_valuesT)
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))
n170 = 170
Power170 = as.data.frame(t(t_valuesPower/repNum))
PowerB170 = cbind(Power170$time_black, n170)
colnames(PowerB170) = c("Power", "N")

PowerH170 = cbind(Power170$time_hisp, n170)
colnames(PowerH170) = c("Power", "N")

#########180###########################################
STEMLong <- function(){
  ### High Level Parameters
  n<-180				# Sample Size
  timePoints<-5		# Number of Time Points
  
  # These are for -.2
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.24 ; b12Val <- -.24  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(1000, STEMLong())
repNum = 1000
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); t_valuesT
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))
n180 = 180
Power180 = as.data.frame(t(t_valuesPower/repNum))
PowerB180 = cbind(Power180$time_black, n180)
colnames(PowerB180) = c("Power", "N")

PowerH180 = cbind(Power180$time_hisp, n180)
colnames(PowerH180) = c("Power", "N")

####################################################
STEMLong <- function(){
  ### High Level Parameters
  n<-190				# Sample Size
  timePoints<-5		# Number of Time Points
  
  # These are for -.2
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.24 ; b12Val <- -.24  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(1000, STEMLong())
repNum = 1000
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); dim(t_valuesT)
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))

n190 = 190
Power190 = as.data.frame(t(t_valuesPower/repNum))
PowerB190 = cbind(Power190$time_black, n190)
colnames(PowerB190) = c("Power", "N")

PowerH190 = cbind(Power190$time_hisp, n190)
colnames(PowerH190) = c("Power", "N")

##########200##########################################
STEMLong <- function(){
  ### High Level Parameters
  n<-200				# Sample Size
  timePoints<-5		# Number of Time Points
  
  # These are for -.2
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.24 ; b12Val <- -.24  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(1000, STEMLong())
repNum = 1000
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); dim(t_valuesT)
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))

n200 = 200
Power200 = as.data.frame(t(t_valuesPower/repNum))
PowerB200 = cbind(Power200$time_black, n200)
colnames(PowerB200) = c("Power", "N")

PowerH200 = cbind(Power200$time_hisp, n200)
colnames(PowerH200) = c("Power", "N")

#######Power Table#############################################
# Just rbind everything and plot it 
PowerTableH_24 = rbind(PowerH150, PowerH160, PowerH170, PowerH180, PowerH190, PowerH200)
PowerTableB_24 = rbind(PowerB150, PowerB160, PowerB170, PowerB180, PowerB190, PowerB200)
PowerTableH_24
PowerTableB_24

library(MASS)
library(psych)
library(reshape)
library(lme4)

###############################################################
# L1 Equation
# ITPis<-pi0s + pi1s*TIMEis + eis
# L2 Equation a
# p0s<-b00 + b01*Black+ b02*Hisp  + r0s
# L2 Equation b
# p0s<-b10 + b11*Black+ b12*Hisp  + r1s
###########################################################################
# Mixed Equation
# ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
###########################################################################

# Let's use the mixed equation
# Unknowns, but we can set the value for.
# ITPis
# b00
# b01;b02
# b10
# b11;b12
# ETHs
# r0s
# r1s
# TIMEis
##########################################################################################
### Since we have only one equation in the mixed equation, we can only have one unknowns
# Let's make the eis is the unknown
######150###################################################################################
STEMLong <- function(){
  ### High Level Parameters
  n<-150				# Sample Size
  timePoints<-5		# Number of Time Points
  
  
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.26 ; b12Val <- -.26  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(10, STEMLong())
repNum = 10
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); t_valuesT
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))
n150 = 150
# Then divide by the total to get the percentage of model that had t-values above 2, which is the power.    
Power150 = as.data.frame(t(t_valuesPower/repNum))
PowerB150 = cbind(Power150$time_black, n150)
colnames(PowerB150) = c("Power", "N")

PowerH150 = cbind(Power150$time_hisp, n150)
colnames(PowerH150) = c("Power", "N")

##############160#####################
STEMLong <- function(){
  ### High Level Parameters
  n<-160				# Sample Size
  timePoints<-5		# Number of Time Points
  
  # These are for -.2
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.26 ; b12Val <- -.26  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(1000, STEMLong())
repNum = 1000
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); t_valuesT
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))
n160 = 160
Power160 = as.data.frame(t(t_valuesPower/repNum))
PowerB160 = cbind(Power160$time_black, 160)
colnames(PowerB160) = c("Power", "N")

PowerH160 = cbind(Power160$time_hisp, 160)
colnames(PowerH160) = c("Power", "N")

#############170#######################################

STEMLong <- function(){
  ### High Level Parameters
  n<-170				# Sample Size
  timePoints<-5		# Number of Time Points
  
  # These are for -.2
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.26 ; b12Val <- -.26  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(1000, STEMLong())
repNum = 1000
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); dim(t_valuesT)
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))
n170 = 170
Power170 = as.data.frame(t(t_valuesPower/repNum))
PowerB170 = cbind(Power170$time_black, n170)
colnames(PowerB170) = c("Power", "N")

PowerH170 = cbind(Power170$time_hisp, n170)
colnames(PowerH170) = c("Power", "N")

#########180###########################################
STEMLong <- function(){
  ### High Level Parameters
  n<-180				# Sample Size
  timePoints<-5		# Number of Time Points
  
  # These are for -.2
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.26 ; b12Val <- -.26  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(1000, STEMLong())
repNum = 1000
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); t_valuesT
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))
n180 = 180
Power180 = as.data.frame(t(t_valuesPower/repNum))
PowerB180 = cbind(Power180$time_black, n180)
colnames(PowerB180) = c("Power", "N")

PowerH180 = cbind(Power180$time_hisp, n180)
colnames(PowerH180) = c("Power", "N")

####################################################
STEMLong <- function(){
  ### High Level Parameters
  n<-190				# Sample Size
  timePoints<-5		# Number of Time Points
  
  # These are for -.2
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.26 ; b12Val <- -.26  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(1000, STEMLong())
repNum = 1000
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); dim(t_valuesT)
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))

n190 = 190
Power190 = as.data.frame(t(t_valuesPower/repNum))
PowerB190 = cbind(Power190$time_black, n190)
colnames(PowerB190) = c("Power", "N")

PowerH190 = cbind(Power190$time_hisp, n190)
colnames(PowerH190) = c("Power", "N")

##########200##########################################
STEMLong <- function(){
  ### High Level Parameters
  n<-200				# Sample Size
  timePoints<-5		# Number of Time Points
  
  # These are for -.2
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.26 ; b12Val <- -.26  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(1000, STEMLong())
repNum = 1000
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); dim(t_valuesT)
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))

n200 = 200
Power200 = as.data.frame(t(t_valuesPower/repNum))
PowerB200 = cbind(Power200$time_black, n200)
colnames(PowerB200) = c("Power", "N")

PowerH200 = cbind(Power200$time_hisp, n200)
colnames(PowerH200) = c("Power", "N")

#######Power Table#############################################
# Just rbind everything and plot it 
PowerTableH_26 = rbind(PowerH150, PowerH160, PowerH170, PowerH180, PowerH190, PowerH200)
PowerTableB_26 = rbind(PowerB150, PowerB160, PowerB170, PowerB180, PowerB190, PowerB200)
PowerTableH_26
PowerTableB_26

library(MASS)
library(psych)
library(reshape)
library(lme4)

###############################################################
# L1 Equation
# ITPis<-pi0s + pi1s*TIMEis + eis
# L2 Equation a
# p0s<-b00 + b01*Black+ b02*Hisp  + r0s
# L2 Equation b
# p0s<-b10 + b11*Black+ b12*Hisp  + r1s
###########################################################################
# Mixed Equation
# ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
###########################################################################

# Let's use the mixed equation
# Unknowns, but we can set the value for.
# ITPis
# b00
# b01;b02
# b10
# b11;b12
# ETHs
# r0s
# r1s
# TIMEis
##########################################################################################
### Since we have only one equation in the mixed equation, we can only have one unknowns
# Let's make the eis is the unknown
######150###################################################################################
STEMLong <- function(){
  ### High Level Parameters
  n<-150				# Sample Size
  timePoints<-5		# Number of Time Points
  
  
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.28 ; b12Val <- -.28  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(10, STEMLong())
repNum = 10
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); t_valuesT
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))
n150 = 150
# Then divide by the total to get the percentage of model that had t-values above 2, which is the power.    
Power150 = as.data.frame(t(t_valuesPower/repNum))
PowerB150 = cbind(Power150$time_black, n150)
colnames(PowerB150) = c("Power", "N")

PowerH150 = cbind(Power150$time_hisp, n150)
colnames(PowerH150) = c("Power", "N")

##############160#####################
STEMLong <- function(){
  ### High Level Parameters
  n<-160				# Sample Size
  timePoints<-5		# Number of Time Points
  
  # These are for -.2
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.28 ; b12Val <- -.28  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(1000, STEMLong())
repNum = 1000
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); t_valuesT
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))
n160 = 160
Power160 = as.data.frame(t(t_valuesPower/repNum))
PowerB160 = cbind(Power160$time_black, 160)
colnames(PowerB160) = c("Power", "N")

PowerH160 = cbind(Power160$time_hisp, 160)
colnames(PowerH160) = c("Power", "N")

#############170#######################################

STEMLong <- function(){
  ### High Level Parameters
  n<-170				# Sample Size
  timePoints<-5		# Number of Time Points
  
  # These are for -.2
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.28 ; b12Val <- -.28  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(1000, STEMLong())
repNum = 1000
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); dim(t_valuesT)
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))
n170 = 170
Power170 = as.data.frame(t(t_valuesPower/repNum))
PowerB170 = cbind(Power170$time_black, n170)
colnames(PowerB170) = c("Power", "N")

PowerH170 = cbind(Power170$time_hisp, n170)
colnames(PowerH170) = c("Power", "N")

#########180###########################################
STEMLong <- function(){
  ### High Level Parameters
  n<-180				# Sample Size
  timePoints<-5		# Number of Time Points
  
  # These are for -.2
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.28 ; b12Val <- -.28  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(1000, STEMLong())
repNum = 1000
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); t_valuesT
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))
n180 = 180
Power180 = as.data.frame(t(t_valuesPower/repNum))
PowerB180 = cbind(Power180$time_black, n180)
colnames(PowerB180) = c("Power", "N")

PowerH180 = cbind(Power180$time_hisp, n180)
colnames(PowerH180) = c("Power", "N")

####################################################
STEMLong <- function(){
  ### High Level Parameters
  n<-190				# Sample Size
  timePoints<-5		# Number of Time Points
  
  # These are for -.2
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.28 ; b12Val <- -.28  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(1000, STEMLong())
repNum = 1000
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); dim(t_valuesT)
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))

n190 = 190
Power190 = as.data.frame(t(t_valuesPower/repNum))
PowerB190 = cbind(Power190$time_black, n190)
colnames(PowerB190) = c("Power", "N")

PowerH190 = cbind(Power190$time_hisp, n190)
colnames(PowerH190) = c("Power", "N")

##########200##########################################
STEMLong <- function(){
  ### High Level Parameters
  n<-200				# Sample Size
  timePoints<-5		# Number of Time Points
  
  # These are for -.2
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.28 ; b12Val <- -.28  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(1000, STEMLong())
repNum = 1000
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); dim(t_valuesT)
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))

n200 = 200
Power200 = as.data.frame(t(t_valuesPower/repNum))
PowerB200 = cbind(Power200$time_black, n200)
colnames(PowerB200) = c("Power", "N")

PowerH200 = cbind(Power200$time_hisp, n200)
colnames(PowerH200) = c("Power", "N")

#######Power Table#############################################
# Just rbind everything and plot it 
PowerTableH_28 = rbind(PowerH150, PowerH160, PowerH170, PowerH180, PowerH190, PowerH200)
PowerTableB_28 = rbind(PowerB150, PowerB160, PowerB170, PowerB180, PowerB190, PowerB200)
PowerTableH_28
PowerTableB_28

library(MASS)
library(psych)
library(reshape)
library(lme4)

###############################################################
# L1 Equation
# ITPis<-pi0s + pi1s*TIMEis + eis
# L2 Equation a
# p0s<-b00 + b01*Black+ b02*Hisp  + r0s
# L2 Equation b
# p0s<-b10 + b11*Black+ b12*Hisp  + r1s
###########################################################################
# Mixed Equation
# ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
###########################################################################

# Let's use the mixed equation
# Unknowns, but we can set the value for.
# ITPis
# b00
# b01;b02
# b10
# b11;b12
# ETHs
# r0s
# r1s
# TIMEis
##########################################################################################
### Since we have only one equation in the mixed equation, we can only have one unknowns
# Let's make the eis is the unknown
######150###################################################################################
STEMLong <- function(){
  ### High Level Parameters
  n<-150				# Sample Size
  timePoints<-5		# Number of Time Points
  
  
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.3 ; b12Val <- -.3  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(10, STEMLong())
repNum = 10
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); t_valuesT
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))
n150 = 150
# Then divide by the total to get the percentage of model that had t-values above 2, which is the power.    
Power150 = as.data.frame(t(t_valuesPower/repNum))
PowerB150 = cbind(Power150$time_black, n150)
colnames(PowerB150) = c("Power", "N")

PowerH150 = cbind(Power150$time_hisp, n150)
colnames(PowerH150) = c("Power", "N")

##############160#####################
STEMLong <- function(){
  ### High Level Parameters
  n<-160				# Sample Size
  timePoints<-5		# Number of Time Points
  
  # These are for -.2
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.3 ; b12Val <- -.3  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(1000, STEMLong())
repNum = 1000
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); t_valuesT
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))
n160 = 160
Power160 = as.data.frame(t(t_valuesPower/repNum))
PowerB160 = cbind(Power160$time_black, 160)
colnames(PowerB160) = c("Power", "N")

PowerH160 = cbind(Power160$time_hisp, 160)
colnames(PowerH160) = c("Power", "N")

#############170#######################################

STEMLong <- function(){
  ### High Level Parameters
  n<-170				# Sample Size
  timePoints<-5		# Number of Time Points
  
  # These are for -.2
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.3 ; b12Val <- -.3  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(1000, STEMLong())
repNum = 1000
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); dim(t_valuesT)
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))
n170 = 170
Power170 = as.data.frame(t(t_valuesPower/repNum))
PowerB170 = cbind(Power170$time_black, n170)
colnames(PowerB170) = c("Power", "N")

PowerH170 = cbind(Power170$time_hisp, n170)
colnames(PowerH170) = c("Power", "N")

#########180###########################################
STEMLong <- function(){
  ### High Level Parameters
  n<-180				# Sample Size
  timePoints<-5		# Number of Time Points
  
  # These are for -.2
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.3 ; b12Val <- -.3  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(1000, STEMLong())
repNum = 1000
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); t_valuesT
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))
n180 = 180
Power180 = as.data.frame(t(t_valuesPower/repNum))
PowerB180 = cbind(Power180$time_black, n180)
colnames(PowerB180) = c("Power", "N")

PowerH180 = cbind(Power180$time_hisp, n180)
colnames(PowerH180) = c("Power", "N")

####################################################
STEMLong <- function(){
  ### High Level Parameters
  n<-190				# Sample Size
  timePoints<-5		# Number of Time Points
  
  # These are for -.2
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.3 ; b12Val <- -.3  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(1000, STEMLong())
repNum = 1000
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); dim(t_valuesT)
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))

n190 = 190
Power190 = as.data.frame(t(t_valuesPower/repNum))
PowerB190 = cbind(Power190$time_black, n190)
colnames(PowerB190) = c("Power", "N")

PowerH190 = cbind(Power190$time_hisp, n190)
colnames(PowerH190) = c("Power", "N")

##########200##########################################
STEMLong <- function(){
  ### High Level Parameters
  n<-200				# Sample Size
  timePoints<-5		# Number of Time Points
  
  # These are for -.2
  b00Val <-6.5
  
  b01Val <- -.1 ; b02Val <- -.1 
  
  # This should be low, because we are not expecting alot of change among the white owmen
  b10Val <- .05
  # These values result in scores of 4.5 for black and Hispanic women
  b11Val <- -.3 ; b12Val <- -.3  # Coefficients for the different ethnicities
  ETHsCategories<-c(1,2,3)
  ETHsLabel<-c("W","B","H")
  ETHsProb<- c(0.50,0.25,.25) # highest prob is whites
  
  # ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites
  
  # Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
  ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
  head(ETHs)
  # If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
  # with particular probabilites.
  Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
  Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
  
  
  l1mean<-0
  l1sd<-.5
  l2errorMean<-0
  # Randomly select correlation between the random intecept and slope error terms
  L2corVal<-runif(1,min=.1,max=0.3)
  (l2errorCor<-matrix (c(1,L2corVal,L2corVal,1),byrow=TRUE,ncol=2,nrow=2))
  # Why take the square root?
  sqrt(l2errorCor)
  
  # What is this doing?  It seems to be doing the same thing as above, but adding an eigenvalue check?
  if(1==0){
    repeat {
      l2errorCor<-matrix(runif(4, min = .1, max = .2), ncol = 2,nrow=2)
      l2errorCor[2,1]<-l2errorCor[1,2]
      diag(l2errorCor)<-1
      l2errorCor
      eigenValueCheck<-sum(eigen(l2errorCor)$values>0)
      if(eigenValueCheck==2) {
        break
      }
    }
    l2errorCor
    sqrt(l2errorCor)
  }
  
  ###############################################################
  # L1 error term is known, ITP is unknown
  ###############################################################
  # error term # This is always random.  
  eis<-as.data.frame(matrix(rnorm(n*timePoints,mean=l1mean,sd=l1sd),nrow=n,ncol=timePoints))
  head(eis)
  
  
  #Create the coefficient matrix for b00
  b00 <-as.data.frame(matrix(rep(b00Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b00)
  
  #Create the coefficient matrix for b01;b02;b03 # Coefficient for Black, Hispanic, and Asian
  b01 <-as.data.frame(matrix(rep(b01Val, n*timePoints),nrow=n, ncol=timePoints))
  b02 <-as.data.frame(matrix(rep(b02Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b01);head(b02)
  #Create the coefficient matrix for b10 # Average slope
  b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
  head(b10)
  
  #Create the coefficient matrix for b11
  b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
  b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
  
  head(b11);head(b12)
  
  # Create the TIME variable # Sort puts each column into order so 1,2,3,4,5 in this case.
  TIMEis <- as.data.frame(matrix(sort(rep(1:timePoints, n)),byrow=FALSE,nrow=n, ncol=timePoints))
  head(TIMEis)
  
  # Create the Level 2 error terms #   
  l2error<-mvrnorm(n, mu = rep(l2errorMean,2), Sigma = sqrt(l2errorCor))
  l2error #
  
  r0s <- as.data.frame(matrix(rep(l2error[,1],timePoints),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r0s)
  
  
  r1s <- as.data.frame(matrix(sort(rep(l2error[,2],timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  head(r1s)
  
  # Mixed Equation
  # ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
  
  ITPis<-(b00 + b01*Black+ b02*Hisp  + r0s) + (b10 + b11*Black+ b12*Hisp  + r1s)*TIMEis + eis
  head(ITPis)
  # Create and ID Variable
  idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
  dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(ITPis)[,-2]))
  names(dat)<-c("id","time","black","hisp","itp")
  head(dat)
  
  #Growth model
  m2<-lmer(itp ~ time * (black+hisp) + (time | id), data=dat)
  # Here I am grabbing the t-values
  t_value = coef(summary(m2))[,"t value"] 
  t_value
}

# Replicate the function above
t_values <- replicate(1000, STEMLong())
repNum = 1000
# Tranpose so I can grab the t-values for the variables of interest
t_valuesT = as.data.frame(t(t_values)); dim(t_valuesT)
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))

n200 = 200
Power200 = as.data.frame(t(t_valuesPower/repNum))
PowerB200 = cbind(Power200$time_black, n200)
colnames(PowerB200) = c("Power", "N")

PowerH200 = cbind(Power200$time_hisp, n200)
colnames(PowerH200) = c("Power", "N")

#######Power Table#############################################
# Just rbind everything and plot it 
PowerTableH_3 = rbind(PowerH150, PowerH160, PowerH170, PowerH180, PowerH190, PowerH200)
PowerTableB_3 = rbind(PowerB150, PowerB160, PowerB170, PowerB180, PowerB190, PowerB200)




