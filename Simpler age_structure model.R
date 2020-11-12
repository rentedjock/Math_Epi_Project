library(tidyverse)
library(tidyselect)
library(deSolve)
library(socialmixr)


#########################
#### Age structure model
#########################

sir.model <- function (times, x, parms) { #SIR model equations
  
  S <- x[sindex]
  I <- x[iindex]
  R <- x[rindex]
  Inc <-x[incindex]
  D <- x[dindex]
  
  N <- S + I + R 
  
  beta <- parms[["p"]]*parms[["c"]] #beta = p*c
  #social distancing
  for(j in 1: length(social.distancing)){
    
    
    if (times %in% sdbreaks[j]:sdbreaks[j+1]) {
      beta <- beta *social.distancing[j]
    }
    
  }
  
  
  
  gamma <- 1/parms[["dur.inf"]] #gamma = 1/dur.inf
  aging <- parms[["aging"]]
  births <- parms[["births"]]
  mu <- diag(parms[["mu"]])
  lambda <-  beta %*% (I/N)
  
  #By 10 age groups
  dS <- -lambda*S + aging%*%S + births%*%N
  dI <- lambda*S - gamma*I + aging%*%I - mu%*%I
  dR <- gamma*I + aging%*%R 
 
  dInc <- lambda*S
  dD <- mu%*%I
  
  return(list(c(dS, dI, dR, dInc, dD)))
}

###########################################
####  Model Parameters
###########################################

#state.names <-c("S", "I", "R", "Inc")

# numebr of age groups
n.i <- 4 # number of age groups

#indices
sindex <- 1:n.i  #these indices help sort out what's saved where
sindex

iindex <- seq(from=max(sindex)+1, by=1, length.out=n.i)
iindex

rindex <- seq(from=max(iindex)+1, by=1, length.out=n.i)
rindex

incindex <- seq(from=max(rindex)+1, by=1, length.out=n.i)
incindex

dindex <- seq(from=max(incindex)+1, by=1, length.out=n.i)
dindex

#age categories, births, aging
age.categories <- c(0,20,40,60)
age.categories

ages <- c(20,20,20,40)
ages

aging <- diag(-1/ages)
aging[row(aging)-col(aging)==1] <- 1/head(ages,-1)
aging

#moving this to days
aging <- aging/365
aging

births <- matrix(0, nrow=n.i, ncol=n.i)
births[1, n.i] <- -aging[n.i, n.i]
births

#Model Parameters

m<-contact_matrix(polymod, countries="United Kingdom", age.limits = age.categories, symmetric = F, split=F)
c <- m$matrix # number of contacts
c

p <- 0.01 # probability of transmission given contact
#https://www.medrxiv.org/content/10.1101/2020.03.03.20028423v3.full.pdf

times <- 1:365 # time step is in days

dur.inf <- 10 # in days

#proportion of cases who die due to COVID per age group
mu <- c(0.004579237, 0.057218699, 0.388159787, 3.777182709)/100
#divided by hundred to convert percentage to proprotions

#social distancing params
sdbreaks <- c(1, 6, 26,max(times))
social.distancing <- c(0.5, 0.5, 0.5) #social distancing



#combining parameters in list
theta <- list(dur.inf = dur.inf, 
              births = births, 
              aging = aging, 
              c = c, 
              p = p, 
              mu=mu)

###########################################
####  Starting Conditions
###########################################

#age group specific underestimation factors from ontario data
under.estimate <- c(3.572082088, 3.371225087,3.767713269,4.111897981)


#Regrouping function to regroup age groups from Ontario data
regroup <- function(x){
  x <- c(sum(x[1:2]),sum(x[3:4]),sum(x[5:6]),sum(x[7:10]) )
  
}


#Ontario population
ont.pop <-unlist(read_csv("population.csv", col_names = F)[,2])


ont.pop <- c(sum(ont.pop[1:4]), sum(ont.pop[5:8]),sum(ont.pop[9:12]),sum(ont.pop[13:21]) )


#infected compartment
init.sir <- read.csv("conposcovidloc.csv")


sirs <-init.sir %>%
  filter(Outcome1== "Not Resolved")%>%
  group_by(Age_Group, Outcome1 )%>%
  summarise(n= n())
init.inf <- regroup(unlist(sirs[,3]))*under.estimate

#Recovered compartment
sirs <-init.sir %>%
  filter(Outcome1== "Resolved")%>%
  group_by(Age_Group, Outcome1 )%>%
  summarise(n= n())

rec.init <- regroup(unlist(sirs[1:10, 3]))*under.estimate

#Death compartment
sirs <-init.sir %>%
  filter(Outcome1== "Fatal")%>%
  group_by(Age_Group, Outcome1 )%>%
  summarise(n= n())


d.init <- regroup(c(0, unlist(sirs[,3])) )


#Susceptible compartment
init.s <- ont.pop- rec.init-d.init-init.inf

###########################################
####  Shock Vaccinate
###########################################
vaccinate <- function(props =rep(0.5, 4)){
  
  if (length(props)!= 4) {
    warning("vaccination vector length not equal to 4, smaller vactor recycled")
  }
  
  if (any(props >1)){ stop("Proportion vaccinated can't be greater than 1")}
  
  init.s <- init.s*(1-props)
  rec.init <- rec.init+init.s*props
  
}


vaccinate()
###########################################
####  Running Model
###########################################

yinit <- c(
  S = init.s,
  I = init.inf,
  R = rec.init,
  Inc = rep(0, 4),
  D= d.init
)

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

traj %>%
  mutate_at(vars(contains("inc")), function(x) x - lag(x)) %>%
  gather(key="group", value="values", -c(time)) %>%
  separate(group, into=c("metric", "index"), sep = "(?<=[a-zA-Z])\\s*(?=[0-9])") %>%
  ggplot(aes(x=time, y=values, color=index)) +
  geom_line(lwd=1) +
  facet_wrap(metric~., scales="free_y") +
  scale_color_brewer(type="qual", palette=2) +
  labs(x = "Time", y="Number of people", color = "Age group")




