#=========================
# Integral Project Models
#-------------------------
# adapted from code of Jessica Metcalf, Sean McMahon, Rob Salguero-Gomez, and Eelke Jongenans from an ESA workshop in 2012.
# Their code, in turn, is adapted from code by Johan Dahlgren from the publication
# Nicole et al. 2011. Interdependent effects of habitat quality and climate on population growth of an endangered plant. J Ecol 99 1211-1218.
#
# The IPM is of the form
# N[y,t+1]=Integal of s[x]*G[y,x]+F[y,x]dx
# where s[x] is the probability of surviving for an individual of size x
# G[y,x] is the probability density for growing to size y given individual is of size x and has survived
# F[y,x] is the density of offsrping of size y produced by individuals of size x

#=======================
# Prep work
#-----------------------

# clear variables etc
rm (list = ls ())
require(IPMpack)
# read in the data from the J Ecology paper. This data includes information about year to year size changes, survivorship, and fecundity.

data("dataIPMpackCryptantha")
data=dataIPMpackCryptantha
#====================
# Estimate survival kernel
#--------------------

# Use a generalized linear model with a binomial family (i.e. logistic regression)
# to determine the probability of suriving as function of size

model.survival=glm(surv~size,data=data,family=binomial)

summary(model.survival) # look at the statistics

# create the survival function

Survival=function(x)predict(model.survival,data.frame('size'=x),type='response')

#plot the survivorship data and the model; using the jitter option to make the values clearer
alpha=min(data$size[which(data$size>0)],data$sizeNext[which(data$sizeNext>0)]) # the minimum size of an individual
beta=max(data$size[which(data$size>0)],data$sizeNext[which(data$sizeNext>0)]) # the maximum size of an individual
xs=seq(alpha,beta,length=100) # a sequence of sizes from min size to max size

#pdf("Week5-A.pdf")
plot(data$size,jitter(data$surv,factor=0.1),pch=21,bg=rgb(1,0,0,0.5))
lines(xs,Survival(xs),lwd=4,col="blue")
#dev.off()()

######################################
# Estimating the growth kernel G(x,y)
######################################

# estimate relationship between mean size today and tomorrow where we assume it is a linear relationship
model.growth=lm(sizeNext~size,data=data)
summary(model.growth) # provides the summary statistics

# plotting the residuals below shows that they are approximately normal.
#pdf("Week5-B.pdf")
hist(model.growth$residuals)
#dev.off()()

# hence, we will assume that one can model the variation in growth as a mean size based trend with gaussian variation around this mean.
b=sd(model.growth$residuals) # the SD of the gaussian
Growth=function(y,x)dnorm(y,mean=predict(model.growth,data.frame('size'=x),type='response'),sd=b) # the growth kernel

# plot the data and the fit
#pdf("Week5-C.pdf")
par(cex.lab=1.5,cex.axis=1.5)
image(xs,xs,t(outer(xs,xs,Growth)),xlab=expression(size[t]),ylab=expression(size[t+1]))
points(data$size,data$sizeNext,pch=21,bg=rgb(0,0,1,0.5))
#dev.off()()

#######################
# Estimating the fecundity kernel
#######################

# use a general linear model assuming poisson distributed
# Should be doing a zero-inflated poisson but keeping it simple for now.
# See http://www.ats.ucla.edu/stat/r/dae/zipoisson.htm

# model.fecundity=glm(fec.seed~size,data=data,family=poisson())
# f.fecundity=function(x)predict(model.fecundity,data.frame('size'=x),type="response")
#
# #plot the seed data and the fit
#
# #pdf("Week5-D.pdf")
# plot(data$size,data$fec.seed,pch=21,bg=rgb(1,0,0,0.5),bty="n",xlab="size",ylab="fecundity")
# lines(xs,f.fecundity(xs),lwd=2)
#dev.off()()

# alternative model with flowering probability & poisson for flowering individuals
model.flowering=glm(fec0~size,data=data,family=binomial)
flowering=function(x)predict(model.flowering,data.frame('size'=x),type="response")

# plotting

plot(data$size,jitter(data$fec0,factor=0.1),pch=21,bg=rgb(1,0,0,0.5))
lines(xs,flowering(xs),lwd=4,col="blue")

# poisson on individuals that flowered

who=which(data$fec0==1)
seed.temp=data$fec1[who]
size.temp=data$size[who]
model.fecundity=glm(seed.temp~size.temp,family=poisson)

fecundity=function(x)exp(model.fecundity$coefficients[1]+model.fecundity$coefficients[2]*x)

#plotting

plot(size.temp,seed.temp,pch=21,bg=rgb(1,0,0,0.5))
lines(xs,fecundity(xs),lwd=4,col="blue")



# estimate probability of seed establishment
establishment.prob=sum(is.na(data$size))/sum(data$fec1,na.rm=TRUE)
# estimate mean recruit size and variation around that mean
require(fitdistrplus)
recruit.size=fitdist(data$sizeNext[is.na(data$size)],"gamma")
denscomp(recruit.size)
dgamma(1:23,shape=coef(recruit.size)[1],rate=coef(recruit.size)[2])





save.image("IPMs-v2.Rdata")
