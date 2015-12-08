# run the R file which fits GLMs to the Salguero et al. (2012) data
source("ipm.R")

# min and maximum sizes
x.min=alpha
x.max=2*beta

# initial population state (written as a vector of sizes)

x=c(beta,beta,beta,beta)

# length of run 

Tf=4

# number of reps

reps=1

# store the final population sizes

final=numeric(reps)

# run simulations

for(i in 1:reps){
t=-1
while((t<Tf)&&(length(x)>0)){
  t=t+1
  # determine whether each individual survived: 0 - died; 1 -survived
  survivors=rbinom(length(x),size=1,prob=Survival(x))
  # determine growth of all individuals 
  sizes=rnorm(length(x),mean = predict(model.growth, data.frame(size = x), type = "response"),
        sd = sd(model.growth$residuals))
  # kill off individuals that grew too large and non-survivors
  temp=which(sizes>x.max)
  sizes[temp]=0
  sizes=sizes*survivors
  # find who flowered
  flowered=rbinom(length(x),size=1,prob=Flowering(x))
  # create an empty vector to hold all the offspring
  kids=c()
  
  
  # pull out the non-zero sized individuals
  temp=which(sizes>0)
  x=sizes[temp]
  print(x)
  }
}

