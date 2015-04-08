#===============================================================================
#-------------------------------------------------------------------------------
#-------------------------Perennial IPM---------------------------------
#-------------------------Individual based model --------------------------
#===============================================================================

# clear variables etc
rm (list = ls ())
# load the IPM build from the IPMs.R file

load("IPMs.Rdata")

###########################
# KEY FUNCTIONS
# s(x) survival function: size dependent
# g(y,x) growth kernel: probability of growing from size x to y
# f(x) fecundity kernel: mean number of offspring from poisson fit times est. probability
# os(y) offspring size kernel: probabily density of offspring sizes
#--------------------------
s=function(x)Survival(x)
g=function(y,x)Growth(y,x)
f=function(x)f.fecundity(x)*establishment.prob
os=function(y)dnorm(y,mean=recruit.size.mean,sd=recruit.size.sd)

# generating function for fecundity.
phi=function(x,s)exp(f(x)*(s-1))


# The probability generating functional for this
# branching process is given by
#
# Psi(h)(x)=phi(x,int os(y) h(y) dy)
#						+s(x)*int g(y,x)*h(y)dy





#######################################
# The generating functionals
# Key observation: except for fecundity all is linear!
#---------------------------------------

# matrix for survival; columns are different seasons
s.mat=s(xs)

# growth matrices corresponding to integration with midpoint rule; already has dx
g.mat=outer(xs,xs,g)*dx

# need to make sure to get rid of the eviction
# in the growth kernel by killing off the evicted
# do this by normalizing g.mat
# and by introducing survival for growth matrix
# new survival growth model is s.evict*g.evict
# which agrees with the orginal kernel s*g
# but has a column stochastic matrix g.evict now

q.evict=1-colSums(g.mat) # note: fair amount eviction
g.evict=g.mat%*%diag(1/(1-q.evict))
s.evict=(1-q.evict)*s.mat

# create offspring size vector. needs to sum to one

os.mat=os(xs)/sum(os(xs))  # without normalization get a sum of 1.010046

# create the generating functional
# need to take transpose to take integrals down columns
Psi=function(h){
	h.new=(s.evict*(t(g.evict)%*%h)+1-s.evict)*phi(xs,os.mat%*%h)
	return(h.new)
}

# iterating any number of years

Psi.iterate=function(T){
	hs=matrix(0,n,T+1)
	for(t in 1:T){
		hs[,t+1]=Psi(hs[,t])
		}
	return(hs)
}

# function to compute probabilities of extinction (not cdf)

probs=function(hs,N){
	hs=hs^N
	T=dim(hs)[2]-1
	pr=matrix(0,n,T)
	for(i in 1:T){pr[,i]=hs[,i+1]-hs[,i]}
	return(pr)
}


# illustrating how things work

# estimating the probability of extinction

out=Psi.iterate(100)
plot(xs,out[,100],type="l",xlab="size",ylab="probability of extinction",lwd=2,bty="n")

# probability of extinction in a given year for certain size classes

ps=probs(out,1)
matplot(t(ps[c(50,60,75),]),lwd=3,lty=1,typ="l",xlab="year",ylab="probabilty",bty="n")


# eigenvalue check
A=g.evict%*%diag(s.evict)+outer(xs,xs,Fecundity)*dx
eigen(A)$values[1]
