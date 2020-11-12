library(tidyverse)
library(tidyselect)
library(deSolve)
library(socialmixr)


# Age structure model
##
sir.model <- function (times, x, parms) { #SIR model equations
  
  S <- x[sindex]
  
  E <- x[eindex]
  I <- x[iindex]
  R <- x[rindex]
  Inc <-x[incindex]
  D <- x[dindex]
  
  N <- S +E+ I + R-D
  
  
  
  beta <- parms[["p"]]*parms[["c"]] #beta = p*c
  gamma <- 1/parms[["dur.inf"]] #gamma = 1/dur.inf
  aging <- parms[["aging"]]
  births <- parms[["births"]]
  mu <- diag(parms[["mu"]])
  
  #social distancing
for(j in 1: length(social.distancing)){
  
    
if (times %in% sdbreaks[j]:sdbreaks[j+1]) {
   beta <- beta *social.distancing[j]
  }
  
}
  
  
  sigma <- 1/parms[["latent"]]
  
  lambda <-  beta %*% (I/N)#divide by N because of age -groups. freq vs. density dependent
  #this is frequency dependent
  
  #By 10 age groups
  dS <- -lambda*S + aging%*%S + births%*%N- diag(rho)%*%S
  dE <- aging %*%E + lambda*S -sigma*E
  dI <- sigma*E - gamma*I + aging%*%I -mu%*%I #gamma*i(1-ifr)
  dR <- gamma*I + aging%*%R +diag(rho)%*%S
  dInc <- lambda*S
  dD <- mu%*%I #gamma*I*ifr
  
  

  
  return(list(c(dS,dE,  dI, dR, dInc, dD   ) ) )
}

#Next we need to define the starting conditions and parameters:

#state.names <-c("S", "I", "R", "Inc")

# numebr of age groups
n.i <- 10 # number of age groups

#indices
sindex <- 1:n.i  #these indices help sort out what's saved where
eindex <- seq(from=max(sindex)+1, by=1, length.out=n.i)
iindex <- seq(from=max(eindex)+1, by=1, length.out=n.i)

rindex <- seq(from=max(iindex)+1, by=1, length.out=n.i)
incindex <- seq(from=max(rindex)+1, by=1, length.out=n.i)
dindex <- seq(from=max(incindex)+1, by=1, length.out=n.i)

#age categories, births, ageing
age.categories <- seq(0, 90, 10)

ages <- rep(10, 10)
aging <- diag(-1/ages)
aging[row(aging)-col(aging)==1] <- 1/head(ages,-1)
births <- matrix(0, nrow=n.i, ncol=n.i)
births[1, n.i] <- -aging[n.i, n.i]

#model parameters
m<-contact_matrix(polymod, age.limits = age.categories, 
                  symmetric = F, split=F,
                  missing.participant.age = "remove", 
                  missing.contact.age = "sample") 

c <- m$matrix # number of contacts
p <- 0.01 # probability of transmission given contact
#https://www.medrxiv.org/content/10.1101/2020.03.03.20028423v3.full.pdf

times <- 1:60 # time step is in weeks
if (!vaccine){times <-1:39}
dur.inf <- 10 # in days
death.time <- 60 #in days
mu <-(1/death.time)*( c(0.00002, 0.00007, 0.0031, 0.00084, 0.00161, 0.00595, 0.0193, 0.0428, 0.078, 0.078) )
#1. Public Health Ontario. COVID-19 Case Fatality, Case Identification, and Attack Rates in Ontario. 2020 May 20;5. 

#Keep everything in days.
# cumulative  deaths 
latent <- 4.6 # in days
#social distancing params
sdbreaks <- c(1, 6, 26,max(times))
social.distancing <- c(1, 0.1, 0.8) #social distancing

#vaccination
efficacy<- 0.8
prop.vacc <- rep(0.9, 10)
rate <- 1/10
#don't model the delay. Just model it as something that happens instantaneously.
#proprotions instead of rates for death

rho <- efficacy*prop.vacc*rate
if (!vaccine){rho <- rep(0, 10)}
theta <- list(dur.inf = dur.inf, 
              births = births, 
              aging = aging, 
              c = c, 
              p = p, 
              mu=mu, 
              latent =latent,
              rho =rho
              ) #combining parameters in list



#initializing the model
ont.pop <-unlist(read_csv("population.csv", col_names = F)[,2])


ont.pop <- c(sum(ont.pop[1:4]), sum(ont.pop[5:8]),sum(ont.pop[9:12]),sum(ont.pop[13:16]), sum(ont.pop[17:21]))

ontpop <-c()
for (j in 1:((nrow(ont.pop)-1)/2)){
  ontpop[j] <- unlist(ont.pop[j,2] +ont.pop[j+1, 2])
  
  
}
pop.init <- ontpop

init.sir <- read.csv("conposcovidloc.csv")



sirs <-init.sir %>%
  filter(Outcome1== "Not Resolved")%>%
  group_by(Age_Group, Outcome1 )%>%
  summarise(n= n())
init.inf <- 10*unlist(sirs[,3])

sirs <-init.sir %>%
  filter(Outcome1== "Resolved")%>%
  group_by(Age_Group, Outcome1 )%>%
  summarise(n= n())

rec.init <- 10*unlist(sirs[1:10, 3])


sirs <-init.sir %>%
  filter(Outcome1== "Fatal")%>%
  group_by(Age_Group, Outcome1 )%>%
  summarise(n= n())


d.init <- c(0, unlist(sirs[,3]))

hh <- read.csv("covidtesting.csv")
Inc <- hh[nrow(hh), 5]

yinit <- c(
  S = ontpop - init.inf-rec.init-d.init, 
  E= init.inf,
  I = init.inf, 
  R = rec.init,
  Inc = rep(10, n.i),
  D= d.init
)
#https://www150.statcan.gc.ca/t1/tbl1/en/tv.action?pid=1710000501

if (F){if (!vaccine){yinit <- c(
  S= ontpop-1, 
  E= rep(0, 10), 
  I= rep(1,10),
  R= rep(0,10),
  Inc = rep(0, 10),
  D= rep(0, 10)
)}

}


#solving equations, storing values in traj data frame
traj <- as.data.frame(ode(y=yinit, times=times, func=sir.model,
                          parms=theta, method="rk4"))


#removing negative values from traj
k <-1  
for (j in 2:nrow(traj)) {
  if ( all(traj [j, ] >0) ){
    k  <- j
  }
}

traj <- traj[1:k,]

#combining age groups
new_vari <- str_c(rep(c("S", "E", "I", "R", "D", "Inc"), each= 3), rep(c("1", "2", "3"), 6) )
 for (k in 1:6){
   
   traj [, 59+3*k] <- traj[,10*k-8 ]+traj[,10*k-7 ]+traj[,10*k-6 ]
   traj [, 60+3*k] <- traj[,10*k-5 ]+traj[,10*k-4 ]+traj[,10*k-3 ]
   traj [, 61+3*k] <- traj[,10*k-2 ]+traj[,10*k-1 ]+traj[,10*k ]+traj[, 10*k+1]
 }



traj2 <-traj[,c(1,62:79 )]
colnames(traj2)<-c("time",new_vari)




#Plot the model outputs:
traj2%>%
  mutate_at(vars(contains("inc")), function(x) x - lag(x)) %>%
  gather(key="group", value="values", -c(time), na.rm = T) %>%
  separate(group, into=c("metric", "index"), sep = "(?<=[a-zA-Z])\\s*(?=[0-9])") %>%
  ggplot(aes(x=time, y=values, color=index)) +
  geom_line(lwd=1) +
  facet_wrap(metric~., scales="free_y") +
  scale_color_brewer(type="qual", palette=2) +
  labs(x = "Time", y="Number of people", color = "Age group")






