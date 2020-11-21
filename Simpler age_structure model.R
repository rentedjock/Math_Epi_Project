library(tidyverse)
library(tidyselect)
library(deSolve)
library(socialmixr)
#for surface graphs
library(plotly)


#########################
#### Age structure model
#########################



###########################################
####  Starting Conditions
###########################################

#age group specific underestimation factors from ontario data
under.estimate <- c(3.574562704, 3.581050923,4.37442076)
#see associated excel file for data from PHO


#Regrouping function to regroup age groups from Ontario data
regroup <- function(x){
  x <- c(sum(x[1:3]),sum(x[4:6]),sum(x[7:10]) )
  
}



####model Initializing values


#Ontario population
ont.pop <-unlist(read_csv("population.csv", col_names = F)[,2])


ont.pop <- c(sum(ont.pop[1:6]),sum(ont.pop[7:12]),sum(ont.pop[13:21]) )


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


model <- function(props, dose){
sir.model <- function (times, x, parms) { #SIR model equations
  
  S <- x[sindex]
  I <- x[iindex]
  R <- x[rindex]
  Inc <-x[incindex]
  D <- x[dindex]
  
  N <- S + I + R 
  
  beta <- parms[["p"]]*parms[["c"]] #beta = p*c
  
  
  
  
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
n.i <- 3 # number of age groups

#indices
sindex <- 1:n.i  #these indices help sort out what's saved where


iindex <- seq(from=max(sindex)+1, by=1, length.out=n.i)


rindex <- seq(from=max(iindex)+1, by=1, length.out=n.i)


incindex <- seq(from=max(rindex)+1, by=1, length.out=n.i)


dindex <- seq(from=max(incindex)+1, by=1, length.out=n.i)


#age categories, births, aging
age.categories <- c(0,30,60)


ages <- c(30,30,50)


aging <- diag(-1/ages)
aging[row(aging)-col(aging)==1] <- 1/head(ages,-1)


#moving this to days
aging <- aging/365


births <- matrix(0, nrow=n.i, ncol=n.i)
births[1, n.i] <- -aging[n.i, n.i]


#Model Parameters

m<-contact_matrix(polymod, countries="United Kingdom", age.limits = age.categories, symmetric = F, split=F)
c <- m$matrix # number of contacts


p <- 0.01 # probability of transmission given contact
#https://www.medrxiv.org/content/10.1101/2020.03.03.20028423v3.full.pdf

times <- 1:365 # time step is in days

dur.inf <- 10 # in days

#proportion of cases who die due to COVID per age group
#from PHO data, see associated excel file
mu <- c(0.01517539, 0.283819849, 3.777182709)/100
#divided by hundred to convert percentage to proprotions




#combining parameters in list
theta <- list(dur.inf = dur.inf, 
              births = births, 
              aging = aging, 
              c = c, 
              p = p, 
              mu=mu)




###########################################
####  Shock Vaccinate
###########################################

  
  if (length(props)!= 2) {
    warning("vaccination vector length not equal to 2, smaller vactor recycled")
  }
  
  if ( sum(props) >1 ){ stop("Proportion can't be greater than 1")}
  
  for (j in seq_along(props)){
    
    init.s[j] <- init.s[j]- props[j]*dose
    rec.init[j] <- rec.init[j] + dose*props[j]
    
  }
  
  init.s[j+1] <- init.s[j+1]- (1-sum(props))*dose
  rec.init[j+1] <- rec.init[j+1]+ (1-sum(props))*dose
  

###########################################
####  Running Model
###########################################

yinit <- c(
  S = init.s,
  I = init.inf,
  R = rec.init,
  Inc = rep(0, 3),
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
return(traj)

}


#Model function returns the trajectory data frame for a given
# number of doses and the given proportion of doses across age groups



#ploting function
plot_my_traj <- function(traj){
  traj %>%
    mutate_at(vars(contains("inc")), function(x) x - lag(x)) %>%
    gather(key="group", value="values", -c(time)) %>%
    separate(group, into=c("metric", "index"), sep = "(?<=[a-zA-Z])\\s*(?=[0-9])") %>%
    ggplot(aes(x=time, y=values, color=index)) +
    geom_line(lwd=1) +
    facet_wrap(metric~., scales="free_y") +
    scale_color_brewer(type="qual", palette=2) +
    labs(x = "Time", y="Number of people", color = "Age group")
  }


#sample plots
traj<- model(props = c(0.2, 0.5), dose = 1e6)
plot_my_traj(traj)


#deaths fucntion gives number of 
#deaths in age group and total number of deaths

deaths<- function(x){
  death <- unlist(x[365, 14:16])
  
  death <-c(death, sum(death))
  return(death)
}

#initializing data frame
results <- data.frame()
counter <-1

#vector of proportionns from 0-100%
#by steps of 1%
#this is proportion of vaccine going to the youngest age group
prop <- seq(0, 1, 0.01)

for ( i in seq_along(prop) ) {
  #this is the proportion of vaccines going to the middle age group
  #the "model" function assumes that the proportion to
  #oldest age group is 1-sum if otehr two proportions
  prop2 <- seq(0, 1- prop[i], by=0.01)
  
  for (j in  seq_along(prop2)){
    
    
    #calculating deaths in each age group for each set of proportions
    # fed into the model
    h <- model(props = c(prop[i], prop2[j]), dose = 1e6)%>%
      deaths()
    
    #storing the proportions in the data frame
    results[counter, 1] <- prop[i]
    results[counter, 2] <- prop2[j]
   
     for (k in 1:4){
       #storing age specific deaths
       #and overall deaths for each specific proportion 
       #in results data frame
       results[counter, k+2] <- h[k]
     } 
    counter <- counter +1
   print(counter)     
  }
  
  
  
}


#converting the results dataframe 
#to matrix format for plotting
ll <-results[, c(1, 2,6)] %>%
  pivot_wider(names_from = V1, values_from = V6 )
ll <-ll[1:100, ]
rownames(ll) <-unlist( ll[, 1])
ll <- ll[, -1]


ll <- matrix(unlist(ll), nrow=100, ncol=101)

#plotting
p <-plot_ly( z= ~ll )%>% add_surface()
 
p
