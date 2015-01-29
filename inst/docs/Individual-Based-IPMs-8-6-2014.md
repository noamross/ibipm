---
title: Establishment of size-structured populations]{Establishment and extinction risk of size-structured populations
author: Sebastian Schreiber
---
Introduction
============

General models and methods
==========================

We consider an individual based model where the set of all possible
individual states is a compact set $X$. Possibilities for this state
space include $X=[a,b]$ for a size-structured population where the
minimum size of individuals is $a$ and the maximum size is $b$,
$X=\{1,2,\dots, k\}$ for a population with $k$ discrete stages,
$X=\{0\}\cup [a,b]$ for a size-structured plant population with seeds
represented by state $0$, or $X=\{1,2,\dots,k\} \times [a,b]$ for
populations structured by age and size. As we only consider populations
with a finite number of individuals, the state of the population at any
point of the time is characterized by the individual states of
individuals within the population and how many individuals are within
each of these states $(n_1,n_2,\dots, n_k)$. For example, if there $n_1$
individuals in states $x_1$, $n_2$ individuals in state $x_2$,..., $n_k$
individuals in state $x_k$, then the state of the population is given by
$(n_1,n_2,\dots, n_k; x_1,x_2,\dots, x_k)$. The set of all possible
population states is
$$\S= \{(\n,\x)=(n_1,\dots, n_k, x_1\dots , x_k): k,n_i\in \N, x_i \in X\} \cup \{0\}$$
where $\N=\{1,2,3,\dots\}$ denotes the natural numbers and $0$ is the
extinction state i.e. no individuals in the population.

Let $\s_t=( \n, \x)\in \S$ be the population state at time $t$. The
dynamics of $\s_t$ are determined by a set of probabilistic rules that
determine the contribution of each individual in the population to the
population in next time step $t+1$. Given let $p_t^x(\n,\y)$ be the
probability that an individual in state $x$ at time $t$ contributes
$n_1$ individuals in state $y_1$, $n_2$ individuals in state $y_2$, etc.
to the population in time $t+1$. These “contributions” include
individuals changing state, having offspring, or dying in which case
there is no contribution with probability $p_t^x(0)$. We assume that
individual updates are independent of one other. Under these
assumptions, $\s_t$ is what is known as a continuous-state branching
process in the probability literature. These branching processes can be
viewed as the demographic stochastic counterparts to integral projection
models (IPMs). However, unlike the classical IPMs, these generalized
branching process provide insights into how demographic stochasticity
influences extinction risk and establishment success. Harris developed a
general theory for these branching processes. Here, we discuss some of
the highlights of that theory and illustrate the theory with an
application to assessing extinction risk for an endangered perennial
plant, Menzies wallower. Since seasonality plays an important role in
this population, we extend the theory of Harris to fluctuating
environments.

Need some assumptions on $p_t^x(\n,\y)$....erg...need to define this
more carefully wrt to a measure

To study the dynamics of $\s_t$ one can either perform individual-based
simulations using the probabilistic update rules or analytical methods
if one is interested in how the moments of the population and the
probability of extinction changes in time. To understand these
quantities, it is useful to introduce the probability generating maps
(PGMs) associated with the probabilistic update rules. These maps take
continuous functions $h:\S\to [0,1]$ to continuous functions of the same
type
$$\Phi^x_t(h) = \sum_{(\n, \y)\in \S} p_t^x(\n,\y) \prod_i h(y_i)^{n_i}$$
I think these have the nice property that they generate all the moments
$$D\Phi_t^x (1)h(y) =$$

The Model
---------

Since Menzies wallflower is a perennial monocarp, we focus on models for
semelparous populations. However, in the Appendices, we describe how to
develop models for iteroparous populations. Furthermore, as the data for
Menzies wallflower is seasonal, we develop the models for periodic
environments. Furthermore, as this population has a mixture of a
discrete state (seeds) and continuos states (maximal leaf lengths of
individual plants) and we want to write down the simplest representation
of the dynamics, we assume that there is a positive, finite measure
$\mu$ on $X$ for defining the demographic kernels. For instance, in the
case of a purely size structured state space $X=[a,b]$, $\mu$ would be
the usual Lebesgue measure on $[a,b]$. Alternatively, as we discuss
further later one, if $X=\{0\} \cup [a,b]$ where $\{0\}$ corresponds to
the seed class, we choose $\mu$ to be a measure that places a weight of
one on $\{0\}$ and Lebesgue measure on $[a,b]$.

For semelparous populations, individuals do one of three things: die,
survive and transition to another individual state, or survive and
reproduce. The events are mutually exclusive. Let $s_t(x)$ be the
probability that an individual in state $x$ survives from time $t$ to
time $t+1$. A surviving individual in state $x$ either reproduces with
probability $r_t(x)$ or grows with the complementary probability
$1-r_t(x)$. For non-reproducing individuals, let $G_t(y,x)$ be the
transition or “growth” kernel that describes the probability of
transitioning from state $x$ to state $y$ for a surviving,
non-reproducing individual at time $t$. Roughly $G_t(x,y)dy$ corresponds
to the probability of an individual of state $x$ transitioning to a
state in the interval $[y,y+dy]$. To describe reproduction, we assume
that the distribution of offspring produced by individual depends on
their state $x$ and offspring sizes are draw independently from a
distribution that only depends on the parents state. For Humboldt Bay
wallflower, this is a reasonable assumption as there are not significant
differences in seed sizes among individuals (a lll are quite small) and,
consequently, it is likely that the size of an offspring in the next
census is independent of the sizes of other offspring from the same
parent (i.e. sizes are primarily environmentally determined). To
describe this reproductive process, let $f_k(x)$ be the probability an
individual in state $x$ has $k$ offspring and $K_F(x,y)$ be the kernel
that describes the probability of that an offspring of an individual in
state $x$ is born into state $y$.

Let $\s_t=( n_1,\dots, n_k, x_1\dots , x_k)\in \S$ be the population
state at time $t$ which updates stochastically as follows

1.  each individual of type $x_i$ dies with probability
    $1-p_T(x_i)-p_F(x_i)$.

2.  each surviving individual of type $x_i$ reproduces with probability
    $\frac{p_F(x_i)}{p_F(x_i)+p_T(x_i)}$ or transitions with the
    complementary probability,

3.  each transitioning individual of type $x_i$ transitions to a state
    in set $A\subset X$ with probability $\int_A K_T(x_i,y)dy$

4.  each reproducing individual of type $x_i$ produces $n$ offspring
    whose states lie in the sets $A_1,A_2,\dots A_n$ with probability
    $f_n(x) \int_{A_1} K_F(x,y)dy\dots \int_{A_n} K_F(x,y)dy$

For models of this type, the dynamics of extinction can be described
explicitly. These dynamics depend o the probability generating
functional for the stochastic process. To describe this probability
generating functional, let $$\phi_x(s) = \sum_{k=0}^\infty f_k s^k$$ be
the probability generating function for the offspring number of an
individual in state $x$. Then the probability generating functional
$\Psi$ for stochastic process $\s_t$ takes continuous functions
$h:X\to [0,1]$ to continuous functions from $X$ to $[0,1]$. This
functional is defined by (see Appendix)
$$\Psi (h) (x) =1-p_T(x)-p_F(x)+  p_T(x) \int_X K_T(x,y) h(y) dy + p_F(x)\phi_x\left( \int_X K_F(x,y) h(y)dy\right)$$
By work of Harris, the probability of being extinct at time $t$ when
there is initially one individual of type $x$ is given by
$$\Psi^t(\mathbf 0)(x):=\underbrace{\Psi \circ \Psi \circ \dots \circ \Psi}_{t\,\,\rm{ times}}({\mathbf 0})(x)$$
where $\mathbf 0$ denotes the zero function i.e. $\mathbf 0(x)=0$ for
all $x$. In particular, the probability of eventual extinction is given
by $$q(x)=\lim_{t\to\infty} \Psi^t (\mathbf 0)(x)$$ It can be shown that
the extinction function $q:X\to [0,1]$ is the smallest solution to the
equation $$\Psi (q) = q$$ The first characterization of $q$ is useful
for numerical implementation. The second characterization is useful for
developing sensitivity and elasticity formulas.

A useful result about $q$ is a threshold theorem due to Harris. Namely,
consider the mean-field IPM for this stochastic process:
$$n_{t+1}(x)=\int M(y,x) n_t(y) dy$$ where
$$M(x,y)=  p_T(x) K_T(x,y)  + p_F(x) \phi_x'(1) K_F(x,y)$$ If the
dominant eigenvalue associated with this IPM is greater than one, then
(under suitable technical conditions) the asymptotic probability of
extinction is strictly less than one i.e. $q(x)<1$ for all $x$.
Alternatively, if the dominant eigenvalue is $\le 1$, then extinction in
inevitable i.e. $q(x)=1$ for all $x$.

Fluctuating environments
------------------------

These ideas are easily extended to fluctuating environments (see
Appendices for details). For the semelparous model, we have the
probabilities $p_T(x,t),p_F(x,t)$, the kernels $K_F(x,y,t),K_T(x,y,t)$,
and the offspring probability generation function $\phi_x(s,t)$ are now
time dependent. These time dependences result in the time-dependent
probability generating functional for the stochastic process:
$$\begin{aligned}
\Psi_t(h)(x)
=&1-p_T(x,t)-p_F(x,t)+  p_T(x,t) \int_X K_T(x,y,t) h(y) dy \\ &+ p_F(x,t)\phi_x\left( \int_X K_F(x,y,t) h(y)dy,t\right)\\
\end{aligned}$$ The only subtlety in using these time-dependent
functionals to determine extinction likelihoods is that one needs to
iterate them backwards. Namely,
$$\Psi_1 \circ \Psi_2 \circ \dots \Psi_{t-1}\circ \Psi_t({\mathbf 0})(x)=\Psi_1\left(\Psi_2 \left( \dots \Psi_{t-1}\left( \Psi_t({\mathbf 0})\right)\right)\right)(x)$$
describes the probability of the process going extinct in $t$ time steps
given there is initially one individual of type $x$. The probability of
eventual extinction is given by
$$\lim_{t\to\infty} \Psi_1 \circ \Psi_2 \circ \dots \Psi_{t-1}\circ \Psi_t({\mathbf 0})(x)$$
Unfortunately, in general, this asymptotic extinction probability can
not be expressed implicitly as the solution of a fixed point equation as
in the case of constant environments. However, in the case of periodic
environments, i.e. $\Psi_{t+k}=\Psi_t$ for all $t$ for some period $k$,
the probability of ultimate extinction must satisfy
$$\Psi_1\circ \Psi_2\circ \dots\circ\Psi_{k-1}\circ \Psi_k (q)(x) = q(x) \mbox{ for all }x$$
For this periodic case, the Appendix proves there is a threshold
theorem. Namely, consider the mean field dynamic
$$n_{t+1}(x)=\int M_t(y,x) n_t(y) dy$$ where
$$M_t(x,y)=  p_T(x,t) K_T(x,y,t)  + p_F(x,t)\frac{\partial  \phi_x}{\partial s}(1,t) K_F(x,y,t)$$
If the dominant eigenvalue associated with the time $k$ mapping is
greater than one, then $q(x)<1$ for all $x$, else $q(x)=1$ for all $x$.
It is natural to conjecture that a similar threshold theorem can be
stated in the case of stationary environments. A statement of this
conjecture is provided in the Appendices.

Numerical implementation
========================

For this part restrict the discussion to $X=[a,b]$ e.g. continuous size
structure. We describe the procedure for approximating a single
numerical functional $\Psi$. In the case of fluctuating environments,
one needs to approximate each time-dependent functional $\Psi_t$
separately. In the case of the seasonal model, this requires estimating
four functionals corresponding to the four seasons.

To approximate a nonlinear functional $\Psi$ , do the following

1.  Discretize the interval $[a,b]$ into $2n$ subintervals of equal
    width $\Delta x = \frac{b-a}{2n}$. This means there will $2n+1$ end
    points $x_i = a+i \Delta x$ for $0\le i \le 2n$. Call the vector of
    these endpoints xs.

2.  For each of the integral operators
    $h(x)\mapsto \int_X K_\star(x,y) h(y)dy$ where $\star = T$ or $F$
    create the matrix approximation using Simpson’s rule. Namely create
    the $2n+1 \times 2n+1$ matrices $\rm{K\star}$ whose $i+1$-th row is
    $$\frac{\Delta x}{6}
    (K_\star(x_i,x_0),4K_\star(x_i,x_1),2K_\star(x_i,x_2),4K_\star(x_i,x_3),\dots, 2K_\star(x_i,x_{2n-2}),4K_\star(x_i,x_{2n-1}),K_\star(x_i,x_{2n}))$$
    If there is no “eviction” in the model, the sums of these rows
    should equal $1$. However, this may not happen. To deal with this
    issue, define $\star$.evict to be vectors that correspond to row
    sums of $\rm{K}\star$ and renormalize $\rm{K}\star$ by multiplying
    it on the left by a diagonal matrix whose diagonal entries are
    $1/\star.\rm{evict}$.

3.  for each of the $p_\star$ functions create the vectors
    $\rm{p\star}=p_\star(xs)*\star\rm{.evict}$.

4.  create the function $\rm{phi}(x,y)=\phi_x(y)$ which is vector
    friendly i.e. given vectors x and y it returns a vector of the same
    length.

5.  create the discretized version of $\Psi$ which takes vectors of
    length $2n+1$ and returns vectors of length $2n+1$. This function is
    defined by
    $$\rm{Psi}(h) = 1 - \rm{pT}-\rm{pF}+\rm{pT}*(AT\o h)+\rm{pF}*\rm{phi}(\rm{xs},AF\o h)$$

To approximate the extinction probabilities up to time $t$, run a for
loop with $q_0$ as a vector of zeros and $q_{t+1}=\rm{Psi}(q_t,t-1)$.
For stationary (e.g. constant or periodic) environments, one gets a good
approximation of the asymptotic probability of extinction for large
enough $t$.

Illustration
============

The two figures below illustrate output from the wallflower model with
disease. Since the eigenvalue for the mean field model is $<1$, eventual
extinction occurs with probability one. Figure [fig:extinction-disease]
illustrates probabilities of extinction depend on the season of the
founding individual as well as the size of the founding individual.
Clearly, initiating a population with larger sized individuals (in
general) result in smaller likelihoods of extinction over finite time
intervals. Also initiating populations in Fall or Winter seems to result
in lower likelihoods of extinction. Figure [fig:histogram-disease]
predicts that populations generated by the first cohort is most likely
to go extinct after 16 years (i.e. 2007) but has a 5% chance of
persisting at least 26 years (i.e. 2017).

![The probabilities of extinction are plotted as a function of the size
of an initial founder. The green curves correspond to the probability of
extinction in the first year through twentieth year. The blue curve is
the probability of eventual extinction. ](Figs/extinction-disease)

[fig:extinction-disease]

![The probabilities of the time to extinction for the lineage produced
by the first cohort of the data set. The solid red curve corresponds to
the median extinction time and the dashed red lines correspond to the 5%
and 95% for extinction times. ](Figs/histogram-disease)

[fig:extinction-disease][fig:histogram-disease]

The next two figures illustrate output from the wallflower model without
disease. Since the eigenvalue for the mean field model is $>1$,
long-term persistence occurs with positive probability.
Figure [fig:extinction-nodisease] illustrates a similar trend to the
disease case, with extinction being now less likely.
Figure [fig:extinction-nodisease2] predicts the necessary population
size to ensure long-term persistence. Populations started with seedlings
need on the order of 10,000 individuals to persist, while populations of
larger sized individuals only need on the order of 100 to persist. The
effect seasonality is most pronounced for populations of seedlings.

![The probabilities of extinction are plotted as a function of the size
of an initial founder. The green curves correspond to the probability of
extinction in the first year through twentieth year. The blue curve is
the probability of eventual extinction. ](Figs/extinction-nodisease)

[fig:extinction-nodisease]

![The necessary population size to ensure long-term persistence as a
function of the size of individuals and the
season](Figs/extinction-nodisease-2)

[fig:extinction-nodisease2]

Appendices
==========

For iteroparous populations, $p_T(x)$ is the probability of surviving
but not reproducing and $p_F(x)$ the probability of surviving and
reproducing. Then $$\begin{aligned}
\Psi(h)(x)&=1-p_T(x)-p_F(x)\\
&+\int_X K_T(x,y) h(y) dy \left(p_T(x) + p_F(x)\phi_x\left( \int_X K_F(x,y) h(y)dy\right)\right)\\
\end{aligned}$$

