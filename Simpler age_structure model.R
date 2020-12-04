#######################################################################################################################
#### Course:      CHL5425: Mathematical Epidemiology of Communicable Diseases
####  Full Project available at https://github.com/rentedjock/Math_Epi_Project
#### Project    
####  Objective:  Ontario has been given 3 millions doses of a COVID-19 vaccine. What is the optimal vaccination strategy                 
####                  that will minimize the number of deaths due to COVID-19 in Ontario?
#### Authors:     Abhinav Thakral and Mackenzie A. Hamilton
####   
#### Date:        23rd November 2020                                                                      
#######################################################################################################################

library(tidyverse)
library(tidyselect)
library(deSolve)
library(socialmixr)
library(plotly)

### Table of Contents
### Part 1: Initializing starting conditions
### Part 2: Writing a function that runs a simple deterministic SIR model with age-assortative mixing
### Part 3: Functions that plot the SIR trajectories and print out key metrics from the trajectories
### Part 4: Running age-assortative COVID-19 SIR model with shock vaccination for all possible vaccination strategies
### Part 5: Plotting surfance plots of SIR model output for all possible vaccination strategies 
### Part 6: Output SIR model metrics for extreme vaccination strategies 

#############################################
#### Part 1: Initializing Starting Conditions
#############################################

### Setting up Ontario population data ###

#read in Ontario population data from 
ont.pop <-unlist(read_csv("population.csv", col_names = F)[,2])

#regroup age groups from the Ontario population data  into age-groups in the format :<30, 30-59, 60+
ont.pop <- c(sum(ont.pop[1:6]),sum(ont.pop[7:12]),sum(ont.pop[13:21]) )



### Setting up initial number of infected, recovered and susceptible individuals in the 
### Ontario population                                                                  

#creating a function that will group number of infected, susceptible and recovered into age groups <30, 30-59 and 60+ 
regroup <- function(x){
  x <- c(sum(x[1:3]),sum(x[4:6]),sum(x[7:10]) )
}

#age-group specific underestimation factors for true number of infected and recovered individuals in Ontario 
#under estimation fator= true cases in Ontario/reported cases in Ontario
under.estimate <- c(3.574562704, 3.581050923,4.37442076)

## Infected compartment ##
init.sir <- read.csv("conposcovidloc.csv")

sirs <-init.sir %>%
  filter(Outcome1== "Not Resolved")%>%
  group_by(Age_Group, Outcome1 )%>%
  summarise(n= n())

#infecteds regrouped according to age-groups and multiplied by underestimation factor
init.inf <- regroup(unlist(sirs[,3]))*under.estimate

## Recovered compartment ##
sirs <-init.sir %>%
  filter(Outcome1== "Resolved")%>%
  group_by(Age_Group, Outcome1 )%>%
  summarise(n= n())
#Recovered regrouped according to age-groups and multiplied by underestimation factor
rec.init <- regroup(unlist(sirs[1:10, 3]))*under.estimate

## Death compartment ##
sirs <-init.sir %>%
  filter(Outcome1== "Fatal")%>%
  group_by(Age_Group, Outcome1 )%>%
  summarise(n= n())

#Deaths regrouped according to age-groups
d.init <- regroup(c(0, unlist(sirs[,3])) )

## Susceptible compartment ##
init.s <- ont.pop- rec.init-d.init-init.inf

#############################################################################################
#### Part 2: Function that runs as simple deterministic SIR Model with age-assortative mixing
#############################################################################################

#This function models the trajectory of COVID-19 in the Ontario given a fixed number 
# of doses of vaccine have been allocated across 3 age groups (<30, 30-59, 60+). 

#dose: numeric vector of length 1. 
#   Number of people for whom the doses are sufficient
#props: numeric vector of length 2. 
#   The first element is the proportion of vaccines given to <30 years
#   The second element is the proportion of vaccines given to 30-59 years
#   The proportion of vaccine allocted to >60 years is by default 1-sum of the other two proportions,

vaccination.model <- function(props, dose){
  
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
    
    dS <- -lambda*S + aging%*%S + births%*%N
    dI <- lambda*S - gamma*I + aging%*%I - mu%*%I
    dR <- gamma*I + aging%*%R 
    
    dInc <- lambda*S
    dD <- mu%*%I
    
    return(list(c(dS, dI, dR, dInc, dD)))
  }
  
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
  ages <- c(30,30,25)
  aging <- diag(-1/ages)
  aging[row(aging)-col(aging)==1] <- 1/head(ages,-1)
  aging <- aging/365
  
  births <- matrix(0, nrow=n.i, ncol=n.i)
  births[1, n.i] <- -aging[n.i, n.i]
  
  # Model Parameters #
  # contact matrix
  m <-contact_matrix(polymod, countries="United Kingdom", age.limits = age.categories, symmetric = F, split=F)
  c <- m$matrix # number of contacts
  p <- 0.03 # probability of transmission given contact
  dur.inf <- 10 # in days
  mu <- c(0.01517539, 0.283819849, 3.777182709)/100 #proportion of cases who die due to COVID per age group from PHO data
  
  times <- 1:365 # time step is in days
  
  # combining model parameters in a list titled theta
  theta <- list(dur.inf = dur.inf, 
                births = births, 
                aging = aging, 
                c = c, 
                p = p, 
                mu=mu)
  
  # Shock vaccinate
  if (length(props)!= 2) {
    warning("vaccination vector length not equal to 2, smaller vactor recycled")
  }
  
  if ( sum(props) >1 ){ stop("Proportion can't be greater than 1")}
  
  init.s2 <- init.s
  rec.init2 <- rec.init
  
  for (j in seq_along(props)){
    #assuming 95% efficacy in preventing infection as well
    init.s2[j] <- init.s2[j]- props[j]*dose*0.95
    rec.init2[j] <- rec.init2[j] + dose*props[j]*0.95
  }
  
  init.s2[j+1] <- init.s2[j+1] - (1-sum(props))*dose*0.95
  rec.init2[j+1] <- rec.init2[j+1] + (1-sum(props))*dose*0.95
  
  #Run model after applying the shock vaccination to the initial susceptible and recovered compartments
  yinit <- c(
    S = init.s2,
    I = init.inf,
    R = rec.init2,
    Inc = rep(0,3),
    D= d.init
  )
  
  traj <- as.data.frame(ode(y=yinit, times=times, func=sir.model,
                            parms=theta, method="rk4"))
  
  #removing negative values from traj
  k <- 1  
  for (j in 2:nrow(traj)) {
    if ( all(traj [j, ] >0) ){
      k  <- j
    }
  }
  
  traj <- traj[1:k,]
  return(traj)
}

######################################################################################################
###### Part 3: Functions that will plot our COVID-19 trajectories under various vacciantion strategies
######################################################################################################

## Creating a graph labeller
#Next few lines are for changing facet labels

variable.names <- list(
  'S' ="Susceptibles",
  'I' ="Infecteds",
  'R' ="Recovereds",
  'Inc' ="Incidence",
  'D' ="Deaths"
)

graph_labeller <- function(variable, value)(
  return(variable.names[value])
) 

vaccine.names <- list(
  '1' ="All vaccines given to individuals 0-29",
  '2' ="All vaccines given to individuals 30-59",
  '3' ="All vaccines given to individuals 60+"
)  

vacc_labeller <- function(variable, value)(
  return(vaccine.names[value])
)

# Creating a function that will plot the trajectory
plot_my_traj <- function(traj){
  traj %>%
    mutate_at(vars(contains("inc")), function(x) x - lag(x)) %>%
    gather(key="group", value="values", -c(time)) %>%
    separate(group, into=c("metric", "index"), sep = "(?<=[a-zA-Z])\\s*(?=[0-9])") %>%
    mutate(metric= factor(metric, 
                          levels= c("S", "I", "R", "Inc", "D") ))%>%
    ggplot(aes(x=time, y=values, color=index)) +
    geom_line(lwd=1) +
    facet_wrap(metric~., scales="free_y", labeller = graph_labeller) +
    scale_color_brewer(type="qual", palette=2,
                       name= "Age Group", 
                       labels= c("<30 years", "30-59 years", "60+ years")) +
    labs(x = "Time", y="Number of People")
}


#Function to plot graphs by specific metrics (S/I/R/D/Inc)

plot_by_metric <- function(ctraj){
  met <- 1:5
  
  for (k in seq_along(met) ){
    if (k==4){
      ktraj <-ctraj %>%
        
        select(time, (3*k-1):(3*k+1), vaccine )%>%
        mutate_at(2:4, function(x){
          
          c(x[1:365]-lag(x[1:365]),x[366:730]-lag(x[366:730]),x[731:1095]-lag(x[731:1095]))
        } #coding this was so hard!!!
        
        )%>%
        
        mutate(total= rowSums (.[2:4]))%>%
        pivot_longer(cols = c(2:4,6), names_to = "index",values_to ="value") }
    
    else {
      ktraj <-ctraj %>%
        select(time, (3*k-1):(3*k+1), vaccine )%>%
        mutate(total= rowSums (.[2:4]))%>%
        pivot_longer(cols = c(2:4,6), names_to = "index",values_to ="value")}
    
    print(ggplot(ktraj, aes(x=time, y=value, color=index)) +
            geom_line(lwd=1) +
            facet_wrap(vaccine~.,  labeller = vacc_labeller) +
            scale_color_brewer(type="qual", palette=2,
                               name= "Age Group", 
                               labels= c("<30 years", "30-59 years", "60+ years", "Total")) +
            labs(x = "Time", y="Number of People"))
  }
  
}

#Function that provides the cumulative number of deaths in age group and overall 
deaths<- function(x){
  death <- unlist(x[365, 14:16])
  
  death <-c(death, sum(death))
  return(death)
}

#Function that provides the peak incidence of infections per age group and overall 
infecteds<- function(x){ 
  infe <- c(max(x[, 5]),
            max(x[, 6]),
            max(x[, 7]), 
            max(x[, 5]+x[, 6]+x[, 7] )  
  )
  return(infe)
}

################################################################################
###### Part 4: Running age-assortative COVID-19 SIR model with shock vaccination
################################################################################

#initializing data frame to store results
results <- data.frame()
counter <- 1

#vector of proportionns from 0-100% by steps of 1%
# this will be the proportion of vaccine going to the youngest age group
prop.young <- seq(0, 1, 0.01)

for ( i in seq_along(prop.young) ) {
  #this is the proportion of vaccines going to the middle age group
  #the "model" function assumes that the proportion to
  #oldest age group is 1-sum if otehr two proportions
  prop.mid <- seq(0, 1- prop.young[i], by=0.01)
  
  for (j in  seq_along(prop.mid)){
    
    #calculating deaths in each age group for each set of proportions fed into the model
    model.output <- vaccination.model(props = c(prop.young[i], prop.mid[j]), dose = 3e6)
    death.toll <- model.output%>%
      deaths()
    peak.inf.toll <- model.output%>%
      infecteds()
    
    #storing the proportions in the data frame
    results[counter, 1] <- prop.young[i]
    results[counter, 2] <- prop.mid[j]
    
    for (k in 1:4){
      #storing age specific deaths
      #and overall deaths for each specific proportion 
      #in results data frame
      results[counter, k+2] <- death.toll[k]
    }
    
    #storing the infecteds
    for (k in 1:4){
      #storing age specific deaths
      #and overall deaths for each specific proportion 
      #in results data frame
      results[counter, k+6] <- peak.inf.toll[k]
    }
    
    counter <- counter +1
    print(counter)   #to see progress  
  }
  
}

write.csv(results, "deaths.csv")

########################################################################################
#### Part 5: Surface plots for SIR model output from all possible vaccination strategies
########################################################################################

### Peak incidence of infection plot ###

peak.inf.plot <- read.csv("deaths.csv")[,  c(1:3, 8:11)]
peak.inf.plot <- peak.inf.plot[,-1]

#converting the results dataframe to matrix format for plotting
peak.inf.plot <-peak.inf.plot[, c(1, 2,6)] %>%
  pivot_wider(names_from = V1, values_from = V10 )
peak.inf.plot <-peak.inf.plot[1:100, ]
rownames(peak.inf.plot) <-unlist(peak.inf.plot[, 1])
peak.inf.plot <- peak.inf.plot[, -1]

peak.inf.plot <- matrix(peak.inf.plot(ll), nrow=100, ncol=101)

Infections <- peak.inf.plot

#plotting
final.infection.plot <-plot_ly(z = ~Infections, colors= "Greys")%>% 
  add_surface()%>%
  layout (
    title ="Peak COVID-19 infection incidence by \n vaccine allocation strategy",
    scene =list( xaxis =list(title="% vaccine to 0-29"), 
                 yaxis = list(title = "% vaccine to 30-59"), 
                 zaxis = list(title = "Peak COVID-19 Incidence"))
  )

final.infection.plot

### Cumulative COIVD-19 Mortality Plot ###

cum.death.plot <-read.csv("deaths.csv")[,  c(1:7)]
cum.death.plot <- cum.death.plot[,-1]

#converting the results dataframe to matrix format for plotting
cum.death.plot <-cum.death.plot[, c(1, 2,6)] %>%
  pivot_wider(names_from = V1, values_from = V6 )
cum.death.plot <- cum.death.plot[1:100, ]
rownames(cum.death.plot) <-unlist( cum.death.plot[, 1])
cum.death.plot <- cum.death.plot[, -1]

cum.death.plot <- matrix(unlist(cum.death.plot), nrow=100, ncol=101)

Deaths <- cum.death.plot

#plotting
final.death.plot <-plot_ly( z= ~Deaths, colors= "Greys")%>% 
  add_surface()%>%
  layout (
    title ="Cumulative COVID-19 deaths by \n vaccine allocation strategies",
    scene =list( xaxis =list(title="% vaccine to 0-29"), 
                 yaxis = list(title = "% vaccine to 30-59", tickvals =c(0, 20,40,60,80)), 
                 zaxis = list(title = "Cumulative Deaths"))
  )
final.death.plot

#plotting the S/I/R/Inc/D metric graphs  
ctraj <- rbind(vaccination.model(c(1, 0), 3e6), vaccination.model(c(0, 1), 3e6), vaccination.model(c(0, 0), 3e6))
ctraj <- ctraj%>%
  mutate (vaccine= rep(c(1,2,3), each =365) )


plot_by_metric(ctraj)
#########################################################################
#### Part 6: Output cumulative metrics for extreme vaccination strategies
#########################################################################

d1 <-vaccination.model(c(1, 0), 3e6)%>%deaths()
d2 <-vaccination.model(c(0, 1), 3e6)%>%deaths()
d3 <-vaccination.model(c(0, 0), 3e6)%>%deaths()

write.csv(data.frame(d1, d2, d3), "agedeaths.csv")

i1 <-vaccination.model(c(1, 0), 3e6)%>%infecteds()
i2 <-vaccination.model(c(0, 1), 3e6)%>%infecteds()
i3 <-vaccination.model(c(0, 0), 3e6)%>%infecteds()

write.csv(data.frame(i1, i2, i3), "maxinfect.csv")


