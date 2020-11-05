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
  N <- S + I + R-D
  
  
  
  beta <- parms[["p"]]*parms[["c"]] #beta = p*c
  gamma <- 1/parms[["dur.inf"]] #gamma = 1/dur.inf
  aging <- parms[["aging"]]
  births <- parms[["births"]]
  mu <- diag(parms[["mu"]])
  for(j in 1: length(social.distancing)){
    
    
  if (times %in% sdbreaks[j]:sdbreaks[j+1]) {
    beta <- beta *social.distancing[j]
    }
  
   }
  sigma <- 1/parms[["latent"]]
  
  lambda <-  beta %*% (I/N)
  
  dS <- -lambda*S + aging%*%S + births%*%N
  dE <- aging %*%E + lambda*S -sigma*E
  dI <- sigma*E - gamma*I + aging%*%I -mu%*%I
  dR <- gamma*I + aging%*%R
  dInc <- lambda*S
  dD <- mu%*%I
  return(list(c(dS,dE,  dI, dR, dInc, dD)))
}

#Next we need to define the starting conditions and parameters:

#state.names <-c("S", "I", "R", "Inc")
n.i <- 4 # number of age groups
sindex <- 1:n.i  #these indices help sort out what's saved where
eindex <- seq(from=max(sindex)+1, by=1, length.out=n.i)
iindex <- seq(from=max(eindex)+1, by=1, length.out=n.i)
rindex <- seq(from=max(iindex)+1, by=1, length.out=n.i)
incindex <- seq(from=max(rindex)+1, by=1, length.out=n.i)
dindex <- seq(from=max(incindex)+1, by=1, length.out=n.i)
age.categories <- c(0,5, 20, 60)
ages <- c(5,20, 40, 30)
aging <- diag(-1/ages)
aging[row(aging)-col(aging)==1] <- 1/head(ages,-1)
births <- matrix(0, nrow=n.i, ncol=n.i)
births[1, n.i] <- -aging[n.i, n.i]


m<-contact_matrix(polymod, age.limits = age.categories, symmetric = T, split=T) # countries="United Kingdom",
c <- m$matrix # number of contacts
p <- 0.1 # probability of transmission given contact
times <- 1:39 # time step is in weeks
dur.inf <- 10 # in days
mu <-1/30*( c(0.0001, 0.001, 0.03, 0.1))
latent <- 10 # in days
sdbreaks <- c(1, 20, 30, max(times))
social.distancing <- c(1, 0.5, 1)
  
theta <- list(dur.inf = dur.inf, 
              births = births, 
              aging = aging, 
              c = c, 
              p = p, 
              mu=mu, 
              latent =latent)

pop.init <- 1e6*(ages/sum(ages))
init.inf <- 1
rec.init <-0
yinit <- c(
  S = pop.init - init.inf, 
  E= rep(init.inf, n.i),
  I = rep(0, n.i), 
  R = rep(0, n.i),
  Inc = rep(0, n.i),
  D= rep(0, n.i)
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

#Plot the model outputs:
traj%>%
  mutate_at(vars(contains("inc")), function(x) x - lag(x)) %>%
  gather(key="group", value="values", -c(time), na.rm = T) %>%
  separate(group, into=c("metric", "index"), sep = "(?<=[a-zA-Z])\\s*(?=[0-9])") %>%
  ggplot(aes(x=time, y=values, color=index)) +
  geom_line(lwd=1) +
  facet_wrap(metric~., scales="free_y") +
  scale_color_brewer(type="qual", palette=2) +
  labs(x = "Time", y="Number of people", color = "Age group")





