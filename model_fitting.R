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




#####################################
### input data to use for fitting ###
#####################################

covid <- data.frame(read.csv("covidtesting.csv"))[seq(5, 272, 7),c(1, 6:8)]
covid$Reported.Date <- 1:39
colnames(covid) <- c("time", colnames(covid)[2:ncol(covid)])
colnames(covid) <- c(colnames(covid)[1:2], "obs", colnames(covid)[4])


age.specific <- read.csv("conposcovidloc.csv")
age.categ.inc <-age.specific %>%
  group_by(Age_Group, Case_Reported_Date)%>%
  summarise(n=n())%>%
  pivot_wider(names_from = Age_Group, values_from= )
  
d <- ggplot(covid, aes(x=time, y=obs))+
  geom_point() 
d




vaccine <- F

source("age_structure.R")
init.state <- yinit
### set up objective function ###

ss.sir <- function(params0, data) {
  data <- data  #data set used for fitting
  t <- data[,"time"] #time points in data set
  cases <- data[,"obs"] #number of cases in data set
  p <- params0 
 theta$p <- p
  out <- as.data.frame(ode(yinit, times=t, sir.model, 
                           theta, method='rk4'))
  
 out <- mutate(out, D= apply(out[, 52:61], 1, sum))
  out <- mutate(out, D=c(0,diff(D))) #calculate inc from cumulative inc
  sse <- sum((out$D-cases)^2) #sum of squared errors
return(sse)
  }

### define starting values and run optimizer ###



##if you put in a large R0, you need to 
#have a large duration of infection
#if duration of infection is smaller
#compared to R0, you get wrong values of the 
#parameters
fit1 <- optimize(interval = c(0.01, 0.2),f= ss.sir, data=covid)
ss.theta<- theta
ss.theta$p <- fit1$minimum


### run model with parameter values from least squares fitting ###

model.pred <- as.data.frame(ode(init.state, times=covid[,"time"], sir.model, 
                                ss.theta, method='rk4'))
model.pred <- mutate(model.pred,  D= apply(model.pred[, 52:61], 1, sum))
model.pred <- mutate(model.pred, D=c(0,diff(D)))
pred.traj <- pivot_longer(model.pred, cols=-time, names_to="state", values_to="value")
ss<- d+
  geom_line(data=subset(pred.traj, state=="D"), aes(x=time, y=value)); ss
