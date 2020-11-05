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
source("age_structure.r")



#####################################
### input data to use for fitting ###
#####################################

covid <- data.frame(read.csv("covidtesting.csv"))[seq(5, 272, 7),c(1, 6:8)]
age.specific <- read.csv("conposcovidloc.csv")
age.categ.inc <-age.specific %>%
  group_by(Age_Group, Case_Reported_Date)%>%
  summarise(n=n())%>%
  pivot_wider(names_from = Age_Group, values_from= )
  
d <- ggplot(covid, aes(x=Reported.Date, y=Total.Cases))+
  geom_point() 
d

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

