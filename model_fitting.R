### check if have required libraries and install if not ###

if(!"deSolve" %in% rownames(installed.packages())){
  install.packages("deSolve", repos = "http://cran.us.r-project.org")
}

if(!"tidyverse" %in% rownames(installed.packages())){
  install.packages("ggplot2" , repos = "http://cran.us.r-project.org")
}

if(!"data.table" %in% rownames(installed.packages())){
  install.packages("data.table", repos = "http://cran.us.r-project.org")
}

if(!"lhs" %in% rownames(installed.packages())){
  install.packages("lhs", repos = "http://cran.us.r-project.org")
}

if(!"bbmle" %in% rownames(installed.packages())){
  install.packages("bbmle", repos = "http://cran.us.r-project.org")
}

require(deSolve)
require(tidyverse)
require(data.table)
require(lhs)
require(bbmle)

########################
### set up SIR model ###
########################

### model equations ###
source("Agestructure.r")



#####################################
### input data to use for fitting ###
#####################################

flu <- data.frame(read.delim("outbreak.txt"))
d <- ggplot(flu, aes(x=time, y=obs))+
  geom_point() +
  labs(x="Time (weeks)", y="Reported cases"); d

#############################
### Least squares fitting ###
#############################

### set up objective function ###

ss.sir <- function(params0, data) {
  data <- data  #data set used for fitting
  t <- data[,"time"] #time points in data set
  cases <- data[,"obs"] #number of cases in data set
  R0 <- exp(params0[1]) 
  dur.inf <- exp(params0[2])
  theta <- c(R0, dur.inf)
  out <- as.data.frame(ode(init.state, times=t, sir.model, 
                           theta, method='euler'))
  out <- mutate(out, Inc=c(0,diff(Inc))) #calculate inc from cumulative inc
  sse <- sum((out$Inc-cases)^2) #sum of squared errors
}

### define starting values and run optimizer ###

params0 <- log(c(1.3,1))
fit1 <- optim(params0,ss.sir, data=flu)
ss.theta<- exp(fit1$par) 
exp(fit1$par)*c(1,7) #print the best-fit params & convert inf.dur to days

### run model with parameter values from least squares fitting ###

model.pred <- as.data.frame(ode(init.state, times=flu[,"time"], sir.model, 
                                exp(fit1$par), method='euler'))
model.pred <- mutate(model.pred, Inc=c(0,diff(Inc)))
pred.traj <- pivot_longer(model.pred, cols=-time, names_to="state", values_to="value")
ss<- d+
  geom_line(data=subset(pred.traj, state=="Inc"), aes(x=time, y=value)); ss

################################
### Latin Hypercube Sampling ###
################################

### set up latin hypercube sample ###
h <- 2000 #choose the number of points to sample
lhs <- maximinLHS(h,2) #draw 'h' samples from uniform dist'n U(0,1) using LH design

### rescale sample based on prior estimates for each parameter ###
R0.min <- 1  # define minimum and maximum values for each of the parameters
R0.max <- 3
dur.inf.min <- 1/7
dur.inf.max <- 3
## create a parameter set by rescaling our simulated latin hypercube sample
params.set = cbind(R0=lhs[,1]*(R0.max-R0.min)+R0.min,
                   dur.inf=lhs[,2]*(dur.inf.max-dur.inf.min)+dur.inf.min)

### run the model for each LHS parameter set and save calibration targets to a data frame ###

lhs.data <- data.frame(matrix(rep(NA,h*5),nrow=h))
for(i in 1:h){
  params <- (params.set[i,])   #select the parameters to use for a particular model run
  times=flu[,"time"]
  out <- as.data.frame(ode(init.state, times, sir.model, 
                           params, method='euler')) #run the model
  peak<-max(diff(out$Inc))  #maximum number of cases that occur during the outbreak
  cum.inc <- out$Inc[max(times)]  #total number of cases that occur over the time period
  peak.week <-which.max(diff(out$Inc)) #week with the largest number of cases
  lhs.data[i, 1:2]<-params[1:2]
  lhs.data[i,3] <- peak
  lhs.data[i,4] <-cum.inc
  lhs.data[i,5] <-peak.week
}
names(lhs.data) <- c(names(params), 'peak', 'cum.inc', 'peak.week')  #data set with param values and model outputs

### calculate calibration targets from the flu data set we're using for fitting ###

peak.data <- max(flu$obs)  # maximum number of cases reported in the data set
cum.inc.data <- sum(flu$obs) # cumulative cases reported in the data set
peak.week.data <- which.max(flu$obs) #week with highest case reports

### define the acceptable target ranges and select parameter sets that fall within these ranges ###

peak.min <- peak.data*0.8
peak.max <- peak.data*1.2
cum.inc.min <- cum.inc.data*0.8
cum.inc.max <- cum.inc.data*1.2
peak.week.min <- peak.week.data-3
peak.week.max <- peak.week.data + 3

lhs.data$isFittingData=0
## select parameter sets that result in outputs that fall within specified calibration targets 
lhs.data[(lhs.data$peak>=peak.min) & (lhs.data$peak<=peak.max) & 
           (lhs.data$cum.inc>=cum.inc.min) & (lhs.data$cum.inc<=cum.inc.max) &   
           (lhs.data$peak.week>=peak.week.min) & (lhs.data$peak.week<=peak.week.max) & 
           !is.na(lhs.data$peak | lhs.data$cum.inc),]$isFittingData=1

table(lhs.data$isFittingData) #list of param sets that fit the data

### plot parameter priors and posteriors ###

par(mfrow=c(1,2))
boxplot(lhs.data$R0,subset(lhs.data,isFittingData==1)$R0,col=c(grey(0.6),2),
        names=c("Prior","Posterior"), main="LHS for 'R0':\nprior/posteriors distributions")

boxplot(lhs.data$dur.inf,subset(lhs.data,isFittingData==1)$dur.inf,col=c(grey(0.6),2),
        names=c("Prior","Posterior"), main="LHS for 'Infection duration (weeks)':
      \nprior/posteriors distributions")


### plot the epidemic trajectories for accepted parameter sets and compare to data ###

lhs.params=subset(lhs.data,isFittingData==1)

lhs.out <- data.frame(matrix(rep(NA,nrow(flu)*nrow(lhs.params)),nrow=nrow(flu)))
for(i in 1:nrow(lhs.params)){   #run the model for the selected parameter sets
  params <- unlist((lhs.params[i,1:2]))
  out <- as.data.frame(ode(init.state, times=flu[,"time"], sir.model, params, method='euler'))
  out <- mutate(out, Inc=c(0,diff(Inc)))
  lhs.out[, i]<-out$Inc
}  


pred.traj <- lhs.out %>%
  mutate(time = row_number()) %>%
  pivot_longer(cols=-time, values_to="value", names_to = "run")


lhs.plot<- d+
  geom_line(data=pred.traj, aes(x=time, y=value, group=run), color="grey"); lhs.plot

#####################################
### Maximum Likelihood Estimation ###
#####################################

### plot binomial distribution for serosurvey example ###

p <- 0.15
n <- 50
k <- seq(0,n,by=1)
prob <- dbinom(x=k, size=n, prob=p)
plot(k, prob,type='h',lwd=5,lend=1, ylab="Probability")

### plot likelihood distribution for serosurvey example ###
k=21
n <- 50
p <- seq(0,1,by=0.001)
plot(p, dbinom(x=k, size=n, prob=p, log=TRUE), ylim=c(-10,-2),ylab="Log-likelihood", type='l')
abline(h=dbinom(x=k,size=n, prob=k/n,log=TRUE)-0.5*qchisq(p=0.95,df=1),col='red')
abline(v=k/n, col='blue')

###########################
### MLE using SIR model ###
###########################

### set up function to output model incidence ###
prediction <- function(params, times) {
  out <- as.data.frame(ode(init.state, times=times, sir.model, params, method='euler'))
  out <- mutate(out, Inc=c(0,diff(Inc))) #calculate incidence from cumulative incidence
  out[,"Inc"] #return incident cases only
} 

### set up likelihood function ###
##
poisson.loglik <- function(params, data) {
  times <- data$time
  pred <- prediction(params, times)
  sum(dpois(x=data$obs[-1], lambda=pred[-1], log=TRUE))
}

### MLE objective function ###
f <- function(log.R0, log.dur.inf) {
  par <- params
  par[c("R0", "dur.inf")] <- exp(c(log.R0,log.dur.inf))
  -poisson.loglik(par,dat)
}

### define data set and initial guesses for parameters ###

dat <- flu  # data set used for fitting
params <- c(R0=NA, dur.inf=NA)
guess <- list(log.R0=log(1.3),log.dur.inf=log(1))  #starting values for R0 and dur.inf

### run MLE ###

fit0 <- mle2(f,start=guess ); fit0 #run the optimizer
fit <- mle2(f,start=as.list(coef(fit0))); fit
summary(fit)
mle.theta <- exp(coef(fit))  #get the MLE point estimate
ci<-suppressWarnings(exp(confint(fit, quietly=TRUE))) #get the confidence intervals
rownames(ci) <- c("R0", "dur.inf")
print(list("best-fit MLE"=mle.theta*c(1,7), "CI"=ci*c(1,7))) #print MLE estimate & convert dur.inf to d

### compare epidemic trajectories using MLE estimates to least squares ###

model.pred.mle <- as.data.frame(ode(init.state, times=flu[,"time"], sir.model, mle.theta, method='euler'))
model.pred.mle <- mutate(model.pred.mle, Inc=c(0,diff(Inc)))
pred.traj.mle <- pivot_longer(model.pred.mle, cols=-time, names_to = "state", values_to = "value")
model.pred.lci <- as.data.frame(ode(init.state, times=flu[,"time"], sir.model, ci[,1], method='euler'))
model.pred.lci <- mutate(model.pred.lci, Inc=c(0,diff(Inc)))
pred.traj.lci <- pivot_longer(model.pred.lci, cols=-time, names_to="state", values_to = "value")
model.pred.uci <- as.data.frame(ode(init.state, times=flu[,"time"], sir.model, ci[,2], method='euler'))
model.pred.uci <- mutate(model.pred.uci, Inc=c(0,diff(Inc)))
pred.traj.uci <- pivot_longer(model.pred.uci, cols=-time, names_to="state", values_to = "value")
ss+
  geom_line(data=subset(pred.traj.mle, state=="Inc"), aes(x=time, y=value), colour="purple", lty=2)+
  geom_line(data=subset(pred.traj.lci, state=="Inc"), aes(x=time, y=value), colour="grey", lty=1)+
  geom_line(data=subset(pred.traj.uci, state=="Inc"), aes(x=time, y=value), colour="grey", lty=1)

