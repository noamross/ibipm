# run the R file which fits GLMs to the Salguero et al. (2012) data
source("ipm.R")


simulator=function(x.start=100,Tf=10,reps=10){
# min and maximum sizes
x.min=alpha*0
x.max=3*beta

# store the final population state and population size

final.state=list()
final.n=numeric(reps)

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
  flowered=rbinom(length(x),size=1,prob=Flowering(x))
  # create an empty vector to hold all the offspring
  kids=c()
  # determine the number of kids and their size
  if(length(flowered)>0){
    lambda=sum(Fecundity(x[flowered]))*3
    number.kids=rpois(1,lambda=lambda)
    if(number.kids>0){
      kids=rgamma(number.kids, shape = coef(recruit.size)[1], rate = coef(recruit.size)[2])
    }
  }
  
  
  # pull out the dead, too small, too large
  temp=which((sizes>x.min)&(sizes<x.max)&(survivors==1))
  x=sizes[temp]
  temp=which((kids>x.min)&(kids<x.max))
  x=c(x,kids[temp])
}
final.state[[i]]=x
final.n[i]=length(x)
}
return(list(state=final.state,n=final.n))
}


# trying this out to recreate some of the curves from Figure 2. 

Tfs=c(1,10,20)
l=length(Tfs)
k=5 # number of size classes
xs=seq(alpha,beta,length=k)
extinct.prob=matrix(0,k,l)
reps=10000
for(j in 1:l){
for(i in 1:k){
  out=simulator(x=xs[i],Tf=Tfs[j],reps=reps)
  extinct.prob[i,j]=length(which(out$n==0))/reps
  print(k*(j-1)+i)
}
}
matplot(xs,extinct.prob,ylim=c(0,1))
xs.temp=xs

save(file="extinction-v2.Rdata",extinct.prob,xs.temp)
