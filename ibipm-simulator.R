# This file performs individual based simulations 
# of individual-based IPMs based on the Salguero et al. (2012)
# data. This simulations are used to estimate extinction probabilities
# numerically the "hard" way. 


# run the R file which provides 
# the GLMs for the Salguero et al. (2012) data
source("ipm.R")


# the simulator function which starts with a single individual 
# of size x.start, runs the simulation Tf time steps or until 
# extinction (which ever comes first), repeats the simulation
# reps number of times. Output are 
# state  - final state of the population in each simulation (a list)
# n - number of individuals at the the end of simulation (a vector)
# time - time at which simulation stopped. Less than Tf only if 
# extinction occurred during the simulation. 

simulator=function(x.start=100,Tf=10,reps=10){

# store the final population state and population size
# and time

final.state=list()
final.n=numeric(reps)
final.time=numeric(reps)

# run simulations

for(i in 1:reps){
t=0
x=x.start
while((t<Tf)&&(length(x)>0)){
  t=t+1
  # determine whether each individual survived: 0 - died; 1 -survived
  survivors=rbinom(length(x),size=1,prob=Survival(x))
  # determine growth of all individuals 
  sizes=rnorm(length(x),mean = predict(model.growth, data.frame(size = x), type = "response"),
        sd = sd(model.growth$residuals))
  # find who flowered
  temp=rbinom(length(x),size=1,prob=as.numeric(Flowering(x)))
  flowered=which(temp==1)
  # create an empty vector to hold all the offsprings
  kids=c()
  # determine the number of kids and their size
  if(length(flowered)>0){
    lambda=sum(Fecundity(x[flowered])) 
    number.kids=rpois(1,lambda=lambda)
    if(number.kids>0){
      kids=rgamma(number.kids, shape = coef(recruit.size)[1], rate = coef(recruit.size)[2])
    }
  }
  
  
  # pull out the dead and negative sized adults
  temp=which((sizes>0)&(survivors==1))
  x=sizes[temp]
  x=c(x,kids)
}
# store the final state, population size, and simulation time
final.state[[i]]=x
final.n[i]=length(x)
final.time[i]=t
}
return(list(state=final.state,n=final.n,time=final.time))
}


# The following code uses the simulator 
# to estimate the exinction probabilities in Figure 2. 
# To get the values in Figure 2b, the simulator
# is adjusted so that lambda-> lambda*3

Tfs=c(1, 10, 20)
l=length(Tfs)
k=5 # number of size classes
xs=seq(alpha,beta,length=k) # equally spaced sizes
extinct.prob=matrix(0,k,l) # for holding the extinction probabilities
reps=10000
for(i in 1:k){
  out=simulator(x.start=xs[i],Tf=max(Tfs)+1,reps=reps)
  for(j in 1:l){
  extinct.prob[i,j]=length(which(out$time<=Tfs[j]))/reps
  }
  print(i)
}

# save data into a file (sims take a while)

save(file="extinction.Rdata",extinct.prob,xs.temp)
