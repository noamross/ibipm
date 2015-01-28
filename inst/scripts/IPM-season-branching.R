#===============================================================================
#-------------------------------------------------------------------------------
#-------------------------  SEASONAL WALLFLOWER IPM   --------------------------
#------------------------------Individual based model --------------------------
#===============================================================================


# change.the.parameters are passed to this code 
# for creating bifurcation diagrams
# if the first coordinate (survival.factor) equals 0 then
# no changes are made


# TURN ON FOR THE BASE RUN. TURN OFF FOR BIF RUNS
rm (list = ls ())
change.the.parameters=c(0,0)

# TURN OFF FOR THE BASE RUN. TURN ON FOR BIF RUNS
survival.factor=change.the.parameters[1];
rust.percent=change.the.parameters[2]

# TURN ON FOR THE BASE RUN. TURN OFF FOR BIF RUNS
load("data.Rdata")
load("IPM-season-stats.Rdata")
load("IPM-season-build.Rdata")

###########################
# KEY FUNCTIONS
# s survival function: size and season dependent
# p probability of flowering (only winter): size dependent
# g growth kernel: size and season dependent
# phi offspring generating function for reproducing individuals
# e size of germinating indvidual in fall density
#--------------------------

s=function(x,season)FuncSurv(x,season)
p=function(x)FuncFlower(x)
g=function(x,y,season)FuncGrowth(x,y,season)
e=function(x)FuncSeedling(x)

disease.factor=1
if(survival.factor>0){
         temp=predict(model.disease,data.frame(X.temp=rust.percent),type='response')
         disease.factor=temp/model.seedcount
         }
            
phi=function(s)exp(disease.factor*FuncFecund()*(s-1))

#################################
# Basic set up for discretization
#--------------------------------

xs=bin.mids # using bin.mids from determinsitic model (500)
dx=bin.mids[2]-bin.mids[1]

a=min(bin.mids) # left endpoint
b=max(bin.mids) # right endpoint
n=length(bin.mids) # number of bins


#######################################
# The generating functionals
# coordinates 1:n is for the continuous part i.e. the function h
# coordinate n+1 is for the seeds i.e. the variable z
# 
# Key observation: except for fecundity all is linear!
#---------------------------------------

# matrix for survival; columns are different seasons
s.mat=matrix(0,n,4)
for(i in 1:4)s.mat[,i]=s(xs,i)
if(survival.factor>0)s.mat[,2]=pmin(1,survival.factor*s.mat[,2])

# vector flowering
p.mat=matrix(0,n,1)
p.mat=p(xs)

# growth matrices corresponding to integration with midpoint rule; already has dx

g.mat=kern.growth[-(n+1),-(n+1),]

# need to make sure to get rid of the eviction
# in the growth kernel by killing off the evicted
# do this by normalizing g.mat
# and by introducing survival for growth matrix
# new survival growth model is s.evict*g.evict
# which agrees with the orginal kernel s*g
# but has a column stochastic matrix g.evict now

g.evict=array(0,dim=c(n,n,4))
s.evict=array(0,dim=c(n,4))
q.evict=array(0,dim=c(n,4))
for(i in 1:4){
	q.evict[,i]=1-colSums(g.mat[,,i])
	g.evict[,,i]=g.mat[,,i]%*%diag(1/(1-q.evict[,i]))
	s.evict[,i]=s.mat[,i]*(1-q.evict[,i])
}

# emergence distribution. normalized to add up to one

e.mat=e(xs)*dx
e.mat=e.mat/sum(e.mat)


# spring and summer generating functionals
# need to take transpose to take integrals down columns
Gsp=function(hz){
	h=hz[1:n]
	z=hz[n+1]
	h.new=1-s.evict[,1]+s.evict[,1]*(t(g.evict[,,1])%*%h)
	z.new=z
	return(c(h.new,z.new))
	}
	
Gsu=function(hz){
	h=hz[1:n]
	z=hz[n+1]
	h.new=1-s.evict[,2]+s.evict[,2]*(t(g.evict[,,2])%*%h)
	z.new=z
	return(c(h.new,z.new))
	}
	
# winter generating functional

Gwi=function(hz){
	h=hz[1:n]
	z=hz[n+1]
	h.new=1-s.evict[,4]+s.evict[,4]*(1-p.mat)*(t(g.evict[,,4])%*%h)+s.evict[,4]*p.mat*phi(z)
	z.new=1
	return(c(h.new,z.new))
	}

# fall generating functional

Gfa=function(hz){
	h=hz[1:n]
	z=hz[n+1]
	h.new=1-s.evict[,3]+s.evict[,3]*(t(g.evict[,,3])%*%h)
	z.new=sum(h*e.mat)
	return(c(h.new,z.new))
	}

# put them all together

G.quarter=list(sp=Gsp,su=Gsu,fa=Gfa,wi=Gwi)

# yearly starting in winter

G.year=function(hz){	
	years=c(3,2,1,4)
	for(i in 1:4){
		hz=G.quarter[[years[i]]](hz)
	}
	return(hz)
}

# iterating any number of years

G.iterate=function(T){
	Wi=matrix(0,n+1,T+1)
	Sp=matrix(0,n+1,T+1)
	Su=matrix(0,n+1,T+1)
	Fa=matrix(0,n+1,T+1)
	Wi[,1]=G.quarter[[4]](rep(0,n+1))
	temp=G.quarter[[1]](rep(0,n+1))
	Sp[,1]=G.quarter[[4]](temp)
	temp=G.quarter[[2]](rep(0,n+1))
	temp=G.quarter[[1]](temp)
	Su[,1]=G.quarter[[4]](temp)
	temp=G.quarter[[3]](rep(0,n+1))
	temp=G.quarter[[2]](temp)
	temp=G.quarter[[1]](temp)
	Fa[,1]=G.quarter[[4]](temp)
	for(t in 1:T){
		Wi[,t+1]=G.year(Wi[,t])
		Sp[,t+1]=G.year(Sp[,t])		
		Su[,t+1]=G.year(Su[,t])
		Fa[,t+1]=G.year(Fa[,t])
	}
	return(list(Wi=Wi,Sp=Sp,Su=Su,Fa=Fa))
}

# function to compute the mean time to extinction

mean.time=function(out){
	Wi=out$Wi
	Sp=out$Sp
	Su=out$Su
	Fa=out$Fa
	temp=colSums(1-t(Wi))+colSums(1-t(Sp))+colSums(1-t(Su))+colSums(1-t(Fa))
	temp=temp/4
	return(temp)
}

# function to compute probabilities of extinction (not cdf)

probs=function(out,N){
	Wi=out$Wi^N
	Sp=out$Sp^N
	Su=out$Su^N
	Fa=out$Fa^N
	T=dim(Fa)[2]-1
	pr=matrix(0,n+1,4*T)
	count=1
	cum=matrix(0,n+1,1)
	for(i in 1:T){
		pr[,4*(i-1)+1]=Wi[,i]-cum
		cum=Wi[,i]
		pr[,4*(i-1)+2]=Sp[,i]-cum
		cum=Sp[,i]
		pr[,4*(i-1)+3]=Su[,i]-cum
		cum=Su[,i]
		pr[,4*(i-1)+4]=Fa[,i]-cum
		cum=Fa[,i]
	}
	return(pr)
	
}

# sensitivity analysis for extinction probability with respect to survival
# INPUT: the extinction probability function q as a vector of length n+1
# for winter
# be sure to compute the derivate function for the reproductive pgf
disease.factor=1
phi.prime=function(s)disease.factor*FuncFecund()*exp(disease.factor*FuncFecund()*(s-1))

Sensitivity=function(q){
	DG=array(0,dim=c(n+1,n+1,4))
	DG[1:n,1:n,1]=diag(s.evict[,1])%*%t(g.mat[,,1])
	DG[n+1,n+1,1]=1
	DG[1:n,1:n,2]=diag(s.evict[,2])%*%t(g.mat[,,2])
	DG[n+1,n+1,2]=1
	DG[1:n,1:n,3]=diag(s.evict[,3])%*%t(g.mat[,,3])
	DG[n+1,1:n,3]=e.mat
	DG[1:n,1:n,4]=diag(s.evict[,4]*(1-p.mat))%*%t(g.mat[,,4])
	DG[1:n,n+1,4]=s.evict[,4]*p.mat*phi.prime(q[n+1])
	DG[n+1,n+1,4]=0
	
	dGeps=array(0,dim=c(n+1,n+1,4))
	dGeps[1:n,1:n,3]=diag(-1+t(g.mat[,,3])%*%q[1:n])
	q.temp=Gfa(q)
	dGeps[1:n,1:n,2]=diag(-1+t(g.mat[,,3])%*%q.temp[1:n])
	q.temp=Gsu(q.temp)
	dGeps[1:n,1:n,1]=diag(-1+t(g.mat[,,3])%*%q.temp[1:n])
	q.temp=Gsu(q.temp)
	dGeps[1:n,1:n,4]=diag(-1+diag(1-p.mat)%*%t(g.mat[,,i]))%*%q.temp[1:n]+p.mat*phi.prime(q.temp[n+1])
	
	A=solve(diag(n+1)-DG[,,4]%*%DG[,,1]%*%DG[,,2]%*%DG[,,3])
	S=array(0,dim=c(n+1,n+1,4))
	S[,,4]=A%*%dGeps[,,4]
	S[,,1]=A%*%(DG[,,4]%*%dGeps[,,1])
	S[,,2]=A%*%(DG[,,4]%*%DG[,,1]%*%dGeps[,,2])
	S[,,3]=A%*%(DG[,,4]%*%DG[,,1]%*%DG[,,2]%*%dGeps[,,3])
	return(S)
}