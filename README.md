library(MASS)
library(psych)
library(reshape)
library(lme4)

###############################################################
# L1 Equation
# ITPis<-pi0s + pi1s*TIMEis + eis
# L2 Equation a
# p0s<-b00 + b01*Black+ b02*Hisp+ b03*Asian  + r0s
# L2 Equation b
# p0s<-b10 + b11*Black+ b12*Hisp+ b13*Asian  + r1s
###########################################################################
# Mixed Equation
# ITPis<-b00+b01*ETHs+ r0s + b10*TIMEis+b11*ETHs*TIMEis+r1s*TIMEis+eis
###########################################################################

# Let's use the mixed equation
# Unknowns, but we can set the value for.
# ITPis
# b00
# b01;b02;b03
# b10
# b11;b12;b13
# ETHs
# r0s
# r1s
# TIMEis
##########################################################################################
### Since we have only one equation in the mixed equation, we can only have one unknowns
# Let's make the eis is the unknown
#########################################################################################
STEMLong <- function(){
### High Level Parameters
n<-150				# Sample Size
timePoints<-5		# Number of Time Points

# Here is where I would change the values.  I would change the reflect the differences for each ethnicity.
b00Val <- 2
b01Val <- -.4 ; b02Val <- -.5 ; b03Val <- .6

b10Val <- .2
b11Val <- -.3 ; b12Val <- -.3 ; b13Val <- .3 # Coefficients for the different ethnicities
ETHsCategories<-c(1,2,3,4)
ETHsLabel<-c("W","B","H","A")
ETHsProb<- c(0.55,0.15,.25,.05) # highest prob is whites

# ETHsProb<- c(0.50,0.30,.20) # W,B, H highest prob is whites

	# Create the ETHs variable # Sampled from the ETHS categories with the probs specied above five times
	ETHs<-as.data.frame(matrix(rep(sample(ETHsCategories,n,replace=TRUE,prob=ETHsProb),timePoints),nrow=n,ncol=timePoints,byrow=FALSE))
	head(ETHs)
	# If a row equals a particular value then that is one.  Colmeans should be close to the probabilities, because they are randomly
	# with particular probabilites.
	Black<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==2,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Black)
	Hisp<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==3,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Hisp)
	Asian<-as.data.frame(matrix(rep(ifelse(ETHs[,1]==4,1,0),timePoints),nrow=n,ncol=timePoints,byrow=FALSE));colMeans(Asian)


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
b03 <-as.data.frame(matrix(rep(b03Val, n*timePoints),nrow=n, ncol=timePoints))
head(b01);head(b02);head(b03)
#Create the coefficient matrix for b10 # Average slope
b10 <-as.data.frame(matrix(rep(b10Val, n*timePoints),nrow=n, ncol=timePoints))
head(b10)

#Create the coefficient matrix for b11
b11 <-as.data.frame(matrix(rep(b11Val, n*timePoints),nrow=n, ncol=timePoints))
b12 <-as.data.frame(matrix(rep(b12Val, n*timePoints),nrow=n, ncol=timePoints))
b13 <-as.data.frame(matrix(rep(b13Val, n*timePoints),nrow=n, ncol=timePoints))

head(b11);head(b12);head(b13)

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

ITPis<-(b00 + b01*Black+ b02*Hisp+ b03*Asian  + r0s) + (b10 + b11*Black+ b12*Hisp+ b13*Asian  + r1s)*TIMEis + eis
head(ITPis)
# Create and ID Variable
idVar<-as.data.frame(matrix(sort(rep(1:n, timePoints)),byrow=TRUE,nrow=n, ncol=timePoints))
dat<-as.data.frame(cbind(stack(idVar)[,-2],stack(TIMEis)[,-2],stack(Black)[,-2],stack(Hisp)[,-2],stack(Asian)[,-2],stack(ITPis)[,-2]))
names(dat)<-c("id","time","black","hisp","asian","itp")
head(dat)

#Growth model
m2<-lmer(itp ~ time * (black+hisp+asian) + (time | id), data=dat)
t_value = coef(summary(m2))[,"t value"] 
t_value
}

t_values <- replicate(10, STEMLong())
t_valuesT = as.data.frame(t(t_values)); t_valuesT
# I think I can use apply somehow
t_valuesTimeEth = as.data.frame(cbind(t_valuesT$`time:black`, t_valuesT$`time:hisp`, t_valuesT$`time:asian`))
colnames(t_valuesTimeEth)  = c("time_black", "time_hisp", "time_asian")

# Changing the t- values that are less than or equal -2, because we are expecting statistically significantly lower ITP scores.
t_valuesOnes = as.data.frame(apply(t_valuesTimeEth, 2, function(x){ifelse(x <= -2, 1, 0)}))
t_valuesPower = as.data.frame(apply(t_valuesOnes, 2, sum))
t_valuesPower
t_valuesTimeEth[t_values <= 2]
# Something like for t-values
p_vals[p_vals < 0.05]


