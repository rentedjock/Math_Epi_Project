###########################################################
#Mathematical Epidemiology Of Communicable Diseases
#
#
#
#Final Assignment
#
#
#Objective:Ontario has been given 3 millions doses of a COVID-19 vaccine.
#What is the optimal vaccination strategy for Ontario that will minimize the number of deaths due to COVID-19?
#
#
#Authors:
#   Abhinav Thakral
#   Mackenzie Anne Hamilton
#
#
#Date: 23rd November 2020
###########################################################
library(tidyverse)
library(tidyselect)
library(deSolve)
library(socialmixr)
library(plotly)


#########################
#### Age structure model
#########################



###########################################
####  Starting Conditions
###########################################

#age group specific underestimation factors from ontario data 
#for each age-group
under.estimate <- c(3.574562704, 3.581050923,4.37442076)
#see associated excel file for data from PHO
#
# under estimation fator= true cases in Ontario/reported cases in Ontario
#1. Public Health Ontario. COVID-19 Case Fatality, Case Identification, and Attack Rates in Ontario. 2020 May 20;5. 
#
#Regrouping function to regroup age groups from Ontario data into age-groups in the format :<30, 30-59, 60+
regroup <- function(x){
  x <- c(sum(x[1:3]),sum(x[4:6]),sum(x[7:10]) )
  
}



####model Initializing values


#Ontario population
ont.pop <-unlist(read_csv("population.csv", col_names = F)[,2])

##regroup age groups from this Ontario data file into age-groups in the format :<30, 30-59, 60+
ont.pop <- c(sum(ont.pop[1:6]),sum(ont.pop[7:12]),sum(ont.pop[13:21]) )


#infected compartment
#file from https://covid-19.ontario.ca/data
init.sir <- read.csv("conposcovidloc.csv")


sirs <-init.sir %>%
  filter(Outcome1== "Not Resolved")%>%
  group_by(Age_Group, Outcome1 )%>%
  summarise(n= n())
#infecteds regrouped according to age-groups and multiplied by underestimation factor
init.inf <- regroup(unlist(sirs[,3]))*under.estimate

#Recovered compartment
sirs <-init.sir %>%
  filter(Outcome1== "Resolved")%>%
  group_by(Age_Group, Outcome1 )%>%
  summarise(n= n())
#Recovered regrouped according to age-groups and multiplied by underestimation factor
rec.init <- regroup(unlist(sirs[1:10, 3]))*under.estimate

#Death compartment
sirs <-init.sir %>%
  filter(Outcome1== "Fatal")%>%
  group_by(Age_Group, Outcome1 )%>%
  summarise(n= n())

#Deaths regrouped according to age-groups
d.init <- regroup(c(0, unlist(sirs[,3])) )


#Susceptible compartment
init.s <- ont.pop- rec.init-d.init-init.inf

#this function models the trajectory of cOVID-19
#given a fixed number of doses of vaccine
# and how we divide them across age-groups
#dose------------ numeric vector of length 1. Number of people for whom the doses are sufficient
#props-----------numeric vector of length 2. The first element is the proportion of vaccines given to <30 years
                                            #the second element is teh proportion of vaccines given to 30-59 years
                                            # the proportion to >60 years is by default 1-sum of other two proportions
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


p <- 0.03 # probability of transmission given contact
#from https://nccid.ca/phac-seir-model-on-covid-19/

times <- 1:365 # time step is in days

dur.inf <- 10 # in days

#proportion of cases who die due to COVID per age group
#from PHO data, see associated excel file
mu <- c(0.01517539, 0.283819849, 3.777182709)/100
#divided by hundred to convert percentage to proprotions
#data from
#1. Public Health Ontario. COVID-19 Case Fatality, Case Identification, and Attack Rates in Ontario. 2020 May 20;5.


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
    #assuming 95% efficacy in preventing infection as well
    init.s[j] <- init.s[j]- props[j]*dose*0.95
    rec.init[j] <- rec.init[j] + dose*props[j]*0.95
    
  }
  
  init.s[j+1] <- init.s[j+1]- (1-sum(props))*dose*0.95
  rec.init[j+1] <- rec.init[j+1]+ (1-sum(props))*dose*0.95
  

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

##########################
######PLOTTING Functions
##########################


# code inspired from https://stackoverflow.com/questions/3472980/how-to-change-facet-labels

#graph labeller

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

#ploting function
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


#sample plots
traj<- model(props = c(0.2, 0.5), dose = 1e6)
plot_my_traj(traj)


#####To plot graphs by metrics (S/I/R/D/Inc)
# code inspired from https://stackoverflow.com/questions/3472980/how-to-change-facet-labels

#Next few lines are for changing facet labels

vaccine.names <- list(
  '1' ="Vaccinating Under 30 years",
  '2' ="Vaccinating 30-59 years",
  '3' ="Vaccinating 60+ years"
)


vacc_labeller <- function(variable, value)(
  return(vaccine.names[value])
)

#To plot graphs by metrics (S/I/R/D/Inc)
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

##########################
######Metric Functions
##########################




#deaths fucntion gives number of 
#deaths in age group and total number of deaths

deaths<- function(x){
  death <- unlist(x[365, 14:16])
  
  death <-c(death, sum(death))
  return(death)
}


#infected fucntion gives peak number of 
#infecteds in age group and peak total number of infecteds

infecteds<- function(x){ 
  infe <- c(max(x[, 5]),
            max(x[, 6]),
            max(x[, 7]), 
            max(x[, 5]+x[, 6]+x[, 7] )  
            )
  return(infe)
}

#end of functions
##############################################################################################################

##getting the results
#initializing data frame to store results
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
    jj <- model(props = c(prop[i], prop2[j]), dose = 3e6)
    h <- jj%>%
      deaths()
    b <- jj%>%
      infecteds()
    
    #storing the proportions in the data frame
    results[counter, 1] <- prop[i]
    results[counter, 2] <- prop2[j]
   
     for (k in 1:4){
       #storing age specific deaths
       #and overall deaths for each specific proportion 
       #in results data frame
       results[counter, k+2] <- h[k]
     }
    
    #storing the infecteds
    for (k in 1:4){
      #storing age specific deaths
      #and overall deaths for each specific proportion 
      #in results data frame
      results[counter, k+6] <- b[k]
    }
    
    counter <- counter +1
   print(counter)   #to see progress  
  }
  
  
  
}
### total number of iterations = 101*51

#storing results because 5151 iterations 
#take a lot of time


write.csv(results, "deaths.csv")
####################################
####making surface plot
####################################



########Infection plot##########

ll<-read.csv("deaths.csv")[,  c(1:3, 8:11)]
ll <- ll[,-1]

#converting the results dataframe 
#to matrix format for plotting
ll <-ll[, c(1, 2,6)] %>%
  pivot_wider(names_from = V1, values_from = V10 )
ll <-ll[1:100, ]
rownames(ll) <-unlist( ll[, 1])
ll <- ll[, -1]


ll <- matrix(unlist(ll), nrow=100, ncol=101)

Infections <-ll
#plotting
p <-plot_ly( z= ~Infections)%>% 
  add_surface()%>%
  layout (
    title ="Infections by Vaccine Distribution Across Ages",
    scene =list( xaxis =list(title="Under 30 years"), 
                 yaxis = list(title = "30-59 years"), 
                 zaxis = list(title = "Infections"))
  )

p





########mortality plot##########

ll<-read.csv("deaths.csv")[,  c(1:7)]
ll <- ll[,-1]

#converting the results dataframe 
#to matrix format for plotting
ll <-ll[, c(1, 2,6)] %>%
  pivot_wider(names_from = V1, values_from = V6 )
ll <-ll[1:100, ]
rownames(ll) <-unlist( ll[, 1])
ll <- ll[, -1]


ll <- matrix(unlist(ll), nrow=100, ncol=101)

Deaths <-ll
#plotting
p <-plot_ly( z= ~Deaths, colors= "YlOrRd")%>% 
  add_surface()%>%
  layout (
    title ="Deaths by Vaccine Distribution Across Ages",
   scene =list( xaxis =list(title="Under 30 years"), 
    yaxis = list(title = "30-59 years"), 
    zaxis = list(title = "Deaths"))
  )
 
p


####################################
####other  plots
####################################

model(c(0, 1), 3e6)%>%plot_my_traj()
model(c(1, 0), 3e6)%>%plot_my_traj()
model(c(0, 0), 3e6)%>%plot_my_traj()

ctraj <- rbind(model(c(1, 0), 3e6),model(c(0, 1), 3e6), model(c(0, 0), 3e6) )
ctraj <- ctraj%>%
  mutate (vaccine= rep(c(1,2,3), each =365) )


plot_by_metric(ctraj)


####################################
####Tables
####################################

a1 <-model(c(1, 0), 3e6)%>%deaths()
a2 <-model(c(0, 1), 3e6)%>%deaths()
a3 <-model(c(0, 0), 3e6)%>%deaths()


write.csv(data.frame(a1, a2, a3), "agedeaths.csv")



i1 <-model(c(1, 0), 3e6)%>%infecteds()
i2 <-model(c(0, 1), 3e6)%>%infecteds()
i3 <-model(c(0, 0), 3e6)%>%infecteds()

write.csv(data.frame(i1, i2, i3), "maxinfect.csv")




