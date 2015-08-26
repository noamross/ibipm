# Individual-based Integral Projection Models: The role of size-structure on extinction risk and establishment success
Sebastian Schreiber and Noam Ross  
April 8, 2015  


$\def\N{\mathbb N}$
$\def\R{\mathbb R}$
$\def\P{\mathbb P}$
$\def\N{\mathbb N}$
$\def\S{\mathcal S}$
$\def\n{\mathbf n}$
$\def\s{\mathbb s}$
$\def\x{\mathbf x}$
$\def\y{\mathbf y}$
$\def\o{\%*\%}$

##Introduction



Computations of extinction probabilities or likelihoods of establishment success lie on opposing sides of a theoretician's coin and are used to address questions in conservation biology, restoration ecology, biological invasions and population genetics. Risks of extinction or establishment failure stem from populations consisting of a finite number of individuals, each of which faces a non-zero risk of mortality on any given day. These risks are shaped, in part, by the size and composition of a population in which individuals may differ in age, size, geographical location, or other important characteristics influencing demography. When population structure is finite-dimensional (i.e. a finite number of age classes, stages, geographical locations), multi-type branching processes serve as the stochastic counterpart of matrix models (A & Ney 19XX, Caswell 2001). These models allow one to examine how extinction risk or establishment success depend on the population size and its distribution across a finite number of states (e.g. ages, size classes, patches). 

When collecting demographic data, one commonly makes continuous measurements (e.g. mass, length) about individuals. Consequently, since their introduction in the early 2000s, there has been a surge of interest in integral projection models (IPMs) that account for continuous population structure (REFs). As IPMs correspond to infinite-dimensional matrix models which are nummerically approximated by finite-dimensional matrix models, many of the standard demographic concepts and methods for matrix models (e.g. stable state distributions, reproductive values, life table response experiements) extend to IPMs (REFs).  

Here, we describe individual-based counterparts of IPMs by making use of the theory of continuous-state branching processes (Harris XX). We present new analytical methods for computing extinction probabilities for these individual-based IPMs. These general methods allow one to examine how extinction probabilities depend on any mixture of  continuous and discrete population structure. We also introduce methods for examining the sensitivity of the extinction probabilities to model parameters. To illustrate the implementation of these methods, we introduce and analyze an individual-based IPM of XXX, which is a stochastic counterpart of the IPM developed by XXX. 



##The General Models and Methods

### The Individual Based IPM
We consider an individual-based model where the set of all possible individual states (e.g. age, size, etc.) lies in a compact metric space $X$. For example, for the standard size-structured IPM, $X=[a,b]$ corresponds to the range of sizes measured in the field. For models with a mixture of age and size stucture, $X$ could be given by $\{1,\dots,T\}\times [a,b]$ where $T$ corresponds to the maximal age of an individual.

As we consider finite populations, the state of the population at any point of the time is characterized by the different sizes $(x_1,x_2,\dots, x_k)$  of individuals within the population and how many individuals there are of each these sizes $(n_1,n_2,\dots, n_k)$. Specifically, if there $n_1$ individuals of size $x_1$, $n_2$ individuals of size $x_2$,..., $n_k$ individuals of size $x_k$, then the state of the population is given by $(n_1,n_2,\dots, n_k; x_1,x_2,\dots, x_k)$. Following Harris (1963), the set of all possible population states is
\[
\S= \{(\n,\x)=(n_1,\dots, n_k, x_1\dots , x_k): k,n_i\in \N, x_i \in X\} \cup \{0\}
\]
where $\N=\{1,2,3,\dots\}$ denotes the natural numbers and $0$ is the extinction state corresponding to no individuals in the population.

Let  $\s_t=( \n, \x)\in \S$ be the population state at time $t$. The dynamics of $\s_t$ are determined by a set of probabilistic rules that determine the contribution of each individual in the population to the population in next time step $t+1$. These "contributions" may correspond to an individual surviving and changing state (e.g. growing in size), or having offspring.  Consistent with standard branching process theory, each individual in the population updates independent of all other individuals. The update rule for an individual of size $x$ is given by a probability measure $m(x,d\s)$ on the state space $\S$. Specically, the probability an individual in state $x$ contributes $\s$ individuals to the population in the next time step where $\s$ lies in a subset $A\subset \S$ equals
\[
\P[\s_1\in A| \s_0 =(1,x)]=\int_A m(x,d\s)
\]
where the left hand side reads "the probability the population state lies $A$ at time $t=1$ after initially having only one individual in state $x$." 

If the population state is currently $\s_t=(n_1,\dots,n_k,x_1,\dots,x_k)$, then the state $\s_{t+1}$ is determined as follows:

1. for each of the $n_1$ individuals in state $x_1$, randomly and indepedently choose the number of replacement individuals from distribution $m(x_1,d\,\s)$.
2. repeat step 1. for the types $x_2,\dots, x_k$.
3. determine the new population state by identifying the states of all individuals and counting the total number of individuals in each of these states. 

This iterative algorithm can be used to create individual based simulations of the individual based IPM. As with any branching process, stochastic realizations of this process, with probability one, either go to extinction in finite time, or the population size grows without bound. 

### Probability generating functionals and extinction probabilities

As with multiple-type branching processes, we can characterize the probabilistic state of the system using probability generating functionals $\Psi$ (pgfs). Our approach differs from Harris (1963) who used moment generating functionals instead. Unlike mgfs, the pgfs allow us to compute how extinction probabilties change in time as well as compute asymptotic extinction probabilities. 

To define the pgf $\Psi$, we introduce the following notation: given a continuous function $h:X\to\R$ and $\s\in \S$, let 
\[
h^\s=
\left\{
\begin{array}{cc}
\prod_{i=1}^k h(x_i)^{n_i}&\mbox{ if }\s=(n_1,\dots,n_k,x_1,\dots,x_k)\\
1 &\mbox{ if }\s=0.
\end{array}
\right.
\]
We define the pgf by 
\[
\Psi(h)(x)=\int h^\s \,m(x,d\s).\tag{*}
\]
In words, $\Psi(h)(x)$ corresponds to the expected value of $h^\s$ due to the contributions of individuals in the next time step from an individual in state $x$. This requires integrating over all possible populations states in the next time step.  

The utility of $\Psi$ for computing extinction probabilities follows from two facts. The first fact stems immediately from the definition: if $h_0$ is the zero function (i.e. $h_0(x)=0$ for all $x$), then 
\[
\Psi(h_0)(x)=\int 1 \,m(d\s,x)= \P[\s_1=0|\s_0=(1,x)]
\]
is the probability the population goes extinct in one time step given that initially it consisted of one individual is state $x$. For the second fact, we define 
\[
\Psi^t(h)=\underbrace{\Psi(\Psi(\dots \Psi}_{\mbox{$t$ times}}(h)\dots))
\]
to be the $t$-fold composition of $\Psi$ with itself. We claim that 
\[
\Psi^t(h_0)(x)=\P[\s_t=0|\s_0=(1,x)]=:E_t(x)\tag{**}
\]
is the probability of extinction by time $t$ given that the population initially consisted of one individual in state $x$. To see why this fact is true, we argue by induction. The first fact proves that the statement holds for $t=1$. Now suppose the assertion holds for $t$, we will show it holds for $t+1$.  On the event that $\s_1=(n_1,\dots,n_k,x_1,\dots,x_k)$ is the population state at time $1$, extinction occurs by time $t+1$ only if each of the lineages of the $n_1+\dots +n_k$ individuals go extinct in the next $t$ time steps. As the fates of these lineages are independent of one another, it follows that 
\[
\begin{aligned}
\P[s_{t+1}=0|\s_0=(1,x),\s_1=(n_1,\dots,n_k,x_1,\dots,x_k)]&= \prod_{i=1}^k \P[\s_t=0|\s_0=(1,x_i)]^{n_i}\\
&=\prod_{i=1}^k\left(\Psi^t(h_0)(x_i)\right)^{n_i}=\Psi^t(h_0)^{\s_1}
\end{aligned}
\]
where the second equality follows from our inductive hypothesis. By the law of total probability
\[
\P[s_{t+1}=0|\s_0=(1,x)]=\int [\Psi^t(h_0)]^\s \,m(d\s,x)
\]
which, by definition, equals $\Psi^{t+1}(0)(x)$ as claimed. 

Hence, equation (**) can be used to compute extinction probabilities iteratively. Furthremore, as individuals update independent of one another, the probability of the population going extinct by time $t$ for any initial condition $\s=(n_1,\dots,n_k,x_1,\dots,x_k)$ equals 
\[
\P[\s_t=0|\s_0=\s]=E_t^\s.
\]
These analytic expressions allow us to efficiently compute extinction probabilities by constructing a numerical approximation of the pgf and iterating it with an initial condition of a zero vector which corresponds to the numerical approximation of the zero function. 

As $E_0(x)\le E_1(x),\le E_2(x)\dots$ for any $x\in X$ and $E_t(x)\le 1$ for all $t$, there is a well defined limit corresponding to the probability of eventual extinction:
\[
E_\infty(x):=\lim_{t\to\infty} E_t(x).
\]
Using moment generating functionals and under suitable technical hypotheses, Harris showed that $E_\infty(x)<1$ for all $x$ if the dominant eigenvalue of the mean-field IPM is greater than one, and $E_\infty(x)=1$ for all $x$ otherwise.



### Sensitivity Formulas

**I corrected the formula below, but it will be much harder to explain and illustrate how to generally implement. Consequently, I am planning on removing this section**

To examine the sensitivity of the extinction probability $E_t(x)$ to a parameter $\alpha$ in an individual-based IPM such as the slope or intercept of the logistic regression for a survival function, the chain rule implies that 
\[
\frac{\partial E_{t+1} }{\partial \alpha}(x)=\frac{\partial\Psi}{\partial \alpha}(E_t)(x)+\frac{\partial \Psi}{\partial h}(E_t)(x)\frac{\partial E_t}{\partial \alpha}(x)
\]
We can use this equation to inductively find these sensitivities as $\frac{\partial E_0}{\partial \alpha}(x)=0$ for all $x$. 

## An Illustration with an Endangered Plant Species

To illustrate how these general methods can be applied to a data set, we develop an individual-based IPM based on **add details about the endangered plant species and the work of Nicols et al. Mention reducing survival by 10% for illustrative purposes**. 

For the mean-field IPM, Nicols et al. had a model of the form 
\[
N_{t+1}(x)=\int_a^b s(y)G(y,x)+e(x)f(y)N_t(y)\,dy
\]

where $s(x)$ is the probability of surviving to the next year for individuals of size $y$, $G(y,x)dy$ is the infinitesimal probability that a surviving individual of size $x$ is size $y$ in the next year, $f(x)$ is the mean number of offspring produced by an individual of size $x$, and $e(x)$ is the infinitesimal probability of an offspring is of size $x$ at the time of the annual census. *PROVIDE brief description of the model selection choices*. 



```r
# The following code loads the functions for the mean field IPM: 

rm (list = ls ())
load("IPMs.Rdata")
g=function(y,x)Growth(y,x)
s=function(x)Survival(x)
f=function(x)f.fecundity(x)*establishment.prob
pf=function(x)1
e=function(y)dnorm(y,mean=recruit.size.mean,sd=recruit.size.sd) 
```

The kernels $s(x)$, $G(y,x)$, and $e(x)$ provide us with all the information that the individual-based IPMs requires for probabilistic updating individuals for survival, growth, and size at first census. The fecundity kernel $f(y)$, however, only specifies the mean number of offspring produced, but for an individual-based IPM, we need the distribution of the number of offspring produced by an individual. Fortunately, when using GLMs for the fecundity data, this extra information comes for free. For example, the fecundity data was modeled using a Poisson family for the GLM.  In which case, the mean number $f(x)$ of offspring also specifies the distribution. More generaly, one might use multi-parameter distributions such as a zero inflated poisson or a negative binomial in which case parameters beyond the mean are needed to specify the distribution of offspring number. 


### Deriving the pgf $\Psi$

To define $\Psi$, we observe that the contributions of an individual of size $x$ to the population in the next time step involves the sum of two independent random variables:  the contribution due to survival and growth and the contribution due to reproduction. We will identify the pgfs, $\Psi_g$ and $\Psi_f$, for each of these processes separately. Then, we make use of a fundamental property of pgfs: 

**The pgf for a sum of independent random variables is the product of the pgfs of these random variables.**

to get that 
\[
\Psi=\Psi_g \times \Psi_f.
\]

For the pgf $\Psi_g$ for survival and growth, $\Psi_g(h)(x)$ corresponds to integrating $h^\s$ over all possible contributions $\s$ from an individual of size $x$ surviving and growing. These contributions are of two types: $\s=\emptyset$ when the individual ies,  and $\s=(1,y)$ when the individual survives and grows to size $y$. The first event occurs with probability $1-s(x)$ and the infinitesimal probability of the second event is $s(x)G(y,x)dy$. As $h^\s=1$ when $\s=\emptyset$ and $h^\s=h(y)$ when $\s=(1,y)$, integrating over all possible contributions due to survival and growth yields 
\[
\Psi_{g}(h)(x)=(1-s(x))\times 1 + s(x)\int h(y)G(y,x)dy.
\]

For the pgf $\Psi_f$ corresponding to contributions due to fecundity, $\Psi_f(h)(x)$ is given by integrating $h^\s$ over all possible states $\s$ corresponding to the offspring produced by an individual of size $x$.  To write this down, we can take advantage of the fact that the reproductive contribution of an individual of size $x$ corresponds to a Poisson number of offspring with mean $f(x)$ whose sizes are drawn *independently* from the same offspring distribution. To make use of this observation, we can condition on the number of offspring an individual has. Lets say our individual of size $x$ has $N$ offspring. Then its contribution is a sum of $N$ independent random variables taking values of the form $\s=(1,y)$. The pgf $\Psi_{kid}$ associated with one of these contributions is given by 
\[
\Psi_{kid}(h)(x)=\int h(y)e(y) \, dy
\]
as $h^{{1,y}}(x)=h(y)$ and $e(y) \, dy$ is the infinitesimal probability of an offspring being of size $y$. Since there are $N$ offspring independent of one another, the pgf for this sum is given by the product of the pgfs:
\[
\Psi_{N\,kids}(h)(x)=\left(\int h(y)e(y) \, dy\right)^N.
\]
Finally, to get $\Psi_f$, we can sum over all possible number of kids weighted by their probabilities:
\[
\Psi_f(h)(x)=\sum \frac{\exp(-f(x))f(x)^N}{N!} \Psi_{kid}(h)(x)^N.
\]
As the pgf for a Poisson distribution with mean $f(x)$ is
\[
\phi_x(s)=\sum_N \exp(-f(x))\frac{(f(x)s)^N}{N!}=\exp(f(x)(s-1))
\]
we get  
\[
\Psi_f(h)(x)=\phi(x, \Psi_{kid}(h)(x))=\phi(x,\int h(y)e(y)dy).
\]
Putting together the two pieces, we get the desired pgf $\Psi$:
\[
\Psi(h)(x)=\left((1-s(x))\times 1 + s(x)\int h(y)G(y,x) \, dy\right)\phi(x,\int h(y)e(y) \, dy).
\]
We note that this argument works did not rely on the form of the offspring pgf $\phi$. Hence, we can use the same expression for alternative offspring distributions as illustrate later. 


### Numerically Implementing and Utilizing $\Psi$

To create this probability generating functional numerically, we discretize the size interval $[\alpha,\beta]$ using $n$ equal sized intervals of width $dx=(\beta-\alpha)/n$. To use the midpoint rule to approximate the integrals, we created a vector $xs$ corresponding to the midpoints of these intervales. Using this vector we discretized the survival function as a vector $s.vec=s(xs)$, the growth kernel as a matrix using the outer product $g.mat=outer(xs,xs,g)$, the fecudity function as a vector $f.vec=f(xs)$, and the offspring size distribution as a vector $e.vec=e(xs)$. 

For the previously described methods to work, it is critical that column sums for growth and the sum of offspring size distribution vector equal one i.e. the columns are probability vectors. For most IPMs, this will not occur automatically due to individuals being evicted from the size interval $[\alpha,\beta]$. There are a variety of ways to handle this issue (Ref). As the  offspring size vectors nearly summed to one, we simply renormalized it so that it summed to one. For the growth matrix, we treated eviction as mortality to be consistent with the model developed by Nicoles et al. To do this, we took one minus the column sums and substracted them from the survival vector and then normalized the column sums so that they added to one. When taking the product of survival and growth, the mean-field IPM is unaffected by this change. 

Using these discretized demographic components and the pgf $\phi$ for fecunidy, we get the discretized pgf $\Psi$ is given by 
\[
\Psi_{discrete}(h)=s.vec\circ (g.mat^T\, h)+1-s.vec)\circ phi(xs,e.vec^T\,h)
\]
where $h$ corresponds to a discretized function i.e. a vector of length $n$, $^T$ denotes the transpose of a matrix or vector,  and $\circ$ denotes element by element multiplication. 


```r
n=100 # number of bins
xs=seq(alpha,2*beta,length=n)
dx=xs[2]-xs[1]
s.vec=s(xs)*0.9
e.vec=e(xs)*dx
e.vec=e.vec/sum(e.vec)
g.mat=outer(xs,xs,g)*dx
q.evict=1-colSums(g.mat) 
g.mat=g.mat%*%diag(1/(1-q.evict))
# Note: most of the eviction is from negative sized individuals
s.vec=(1-q.evict)*s.vec
phi=function(x,s)exp(f(x)*(s-1))
Psi=function(h){
	h.new=(s.vec*(t(g.mat)%*%h)+1-s.vec)*(phi(xs,e.vec%*%h))
	return(h.new)
}
Psi.iterate=function(T){
	hs=matrix(0,T+1,n)
	for(t in 1:T){
		hs[t+1,]=Psi(hs[t,])
		}
	return(hs)
}

# function for estimating asymptotic extinction probability
AE=function(tol=0.0001){
  h=rep(0,n)
  diff=1
  while(diff>tol){
    h.old=h
    h=Psi(h)
    diff=sum(abs(h.old-h))
  }
  return(h)
}
```

The figure below illustrates how extinction probabilities $E_t(x)$ vary with size over a $100$ year time frame. Intuitively, this figure illustrates that the probability of extinction decreases with the size of individual initially founding the population, and extinction probabilities increase over time. Furthermore, this figure illustrates that $E_t(x)$ are approaching limiting extinction probabilities $E_{\infy}(x)$ which are less than one. This stems from the fact that the dominant eigenvalue of the mean-field IPM is greater than one. To estimate this asymptotic probability of extinction one can  iterate $\Psi$ until a tolerance condition is met e.g. $|E_{t+1}(x)-E_t(x)|\le \epsilon$ for all $x$.  


```r
T=100
out=Psi.iterate(T)
times=c(1,5,10,25,30,40)+1
matplot(xs,t(out[times,]),xlim=c(alpha,beta),type="l",lty=1,xlab="size",ylab="extinction probabilities",bty="n",lwd=2)
legend("topright",c(as.character(times-1),expression(Inf)),col=c(1:length(times),"black"),lty=1,lwd=c(rep(2,length(times)),5),bty="n")
# adding the asymptotic extinction probability
h=AE()
lines(xs,h,type="l",xlab="size",ylab="asymptotic extinction probability",lwd=5,xlim=c(alpha,beta))
```

![](IBIPM_files/figure-html/unnamed-chunk-3-1.png) 

To scale things up to an entire population, recall that for a population initially in state $\s_0$, the extinction probability at time $t$ is $E_t^{\s_0}$. For example, the figure below illustrates how the extinction probabilities over time vary for a population with initially $100$ individuals of the smallest size $\alpha$ versus a population with $5$ or $6$ individuals of the largest size $\beta$. This figure suggests that, from the extinction risk perspective, about $5$ or $6$ large individuals are equivalent to $100$ small individuals. 


```r
#Inputs: N vector of length $n$ which specifies the initial number of individuals in each size bin, and T the length of the time to look things over.
#Output:a vector of length T+1 corresponding to the probabilities of extinction by time t or sooner
PVA=function(N,T){
	hs=Psi.iterate(T)
	output=numeric(T+1)
	for(t in 1:(T+1)){
		output[t]=prod(hs[t,]^N)
	}
	return(output)
}
T=75
N=rep(0,n)
N[1]=100
out1=PVA(N,T)
N=rep(0,n)
N[50]=6
out2=PVA(N,T)
N=rep(0,n)
N[50]=5
out3=PVA(N,T)
out=cbind(out1,out2,out3)
matplot(0:T,out,typ="l",lwd=3,bty="n",xlab="time t",ylab="probability of extinction by year t")
legend("topleft",c("100 small","6 large","5 large"),lty=1:3,col=1:3,lwd=3,bty="n")
```

![](IBIPM_files/figure-html/unnamed-chunk-4-1.png) 