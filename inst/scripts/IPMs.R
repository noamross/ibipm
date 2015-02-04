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

# read in the data from the J Ecology paper. This data includes information about year to year size changes, survivorship, and fecundity.

data=as.data.frame(read.csv("data/IPM-data.csv"))

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

pdf("Week5-A.pdf")
plot(data$size,jitter(data$surv,factor=0.1),pch=21,bg=rgb(1,0,0,0.5))
lines(xs,Survival(xs),lwd=4,col="blue")
dev.off()

######################################
# Estimating the growth kernel G(x,y)
######################################

# estimate relationship between mean size today and tomorrow where we assume it is a linear relationship
model.growth=lm(sizeNext~size,data=data)
summary(model.growth) # provides the summary statistics

# plotting the residuals below shows that they are approximately normal.
pdf("Week5-B.pdf")
hist(model.growth$residuals)
dev.off()

# hence, we will assume that one can model the variation in growth as a mean size based trend with gaussian variation around this mean.
b=sd(model.growth$residuals) # the SD of the gaussian
Growth=function(y,x)dnorm(y,mean=predict(model.growth,data.frame('size'=x),type='response'),sd=b) # the growth kernel

# plot the data and the fit
pdf("Week5-C.pdf")
par(cex.lab=1.5,cex.axis=1.5)
image(xs,xs,t(outer(xs,xs,Growth)),xlab=expression(size[t]),ylab=expression(size[t+1]))
points(data$size,data$sizeNext,pch=21,bg=rgb(0,0,1,0.5))
dev.off()

#######################
# Estimating the fecundity kernel
#######################

# use a general linear model assuming poisson distributed
# Should be doing a zero-inflated poisson but keeping it simple for now.
# See http://www.ats.ucla.edu/stat/r/dae/zipoisson.htm

model.fecundity=glm(fec.seed~size,data=data,family=poisson())
f.fecundity=function(x)predict(model.fecundity,data.frame('size'=x),type="response")

#plot the seed data and the fit

pdf("Week5-D.pdf")
plot(data$size,data$fec.seed,pch=21,bg=rgb(1,0,0,0.5),bty="n",xlab="size",ylab="fecundity")
lines(xs,f.fecundity(xs),lwd=2)
dev.off()

# estimate probability of seed establishment
establishment.prob=sum(is.na(data$size))/sum(data$fec.seed,na.rm=TRUE)
# estimate mean recruit size and variation around that mean
recruit.size.mean=mean(data$sizeNext[is.na(data$size)])
recruit.size.sd=sd(data$sizeNext[is.na(data$size)])
# put together the fecundity function
Fecundity=function(y,x){ establishment.prob*dnorm(y,mean=recruit.size.mean,sd=recruit.size.sd)*f.fecundity(x)
  }



###############
# The Entire Kernel
################

K=function(y,x)Survival(x)*Growth(y,x)+Fecundity(y,x)

#####################################
# Discretizing the IPM
#######################################

# parameters

n=101   # subdivisions


#  create discretized size vector
# create a Gauss quadrature version
# http://www.damtp.cam.ac.uk/lab/people/sd/lectures/nummeth98/integration.htm

xs=seq(alpha,beta,length=n) # recall alpha and beta are min/max sizes
dx=xs[2]-xs[1] # dx increment in sizes
xs.left=xs[-n]+(1/2-sqrt(3)/6)*dx
xs.right=xs[-n]+(1/2+sqrt(3)/6)*dx
IPM.left=(outer(xs.left,xs.left,K)*dx)
IPM.right=(outer(xs.right,xs.right,K)*dx)
IPM=(IPM.left+IPM.right)/2


#plot kernel
pdf("Week5-E.pdf")
image(xs[-n],xs[-n],IPM)
dev.off()


#############################################
# analyze IPM with the usual matrix machinery
#############################################
out=eigen(IPM)
lambda=Re(out$values[1])
v=Re(out$vectors[,1])
v=v/sum(v)
w=Re(eigen(t(IPM))$vectors[,1])
w=w/sum(v*w)

pdf("Week5-F.pdf")
par(cex.axis=1.5,cex.lab=1.5,mar=c(4,4,1,4))
plot(xs[-n],v,type="l",ylab="stable size distribution",xlab="size",lwd=3,col="black")
par(new=TRUE)
plot(xs[-n],w,type="l",yaxt="n",xlab="",ylab="",lwd=3,col="red")
axis(4)
mtext("reproductive values",side=4,line=3,cex=1.5)
legend("top",c("stable size","reproductive values"),col=c("black","red"),lwd=3,bty="n",cex=1.5)
dev.off()


S=w%o%v
E=S*IPM/lambda

pdf("Week5-G.pdf")
image(xs[-n],xs[-n],E/dx^2,xlab="size at t",ylab="size at t+1")
contour(xs[-n],xs[-n],E/dx^2,add=TRUE)
dev.off()
