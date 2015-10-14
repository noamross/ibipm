## ----components, cache=FALSE, fig.width=7, fig.height=5, echo=FALSE, fig.cap="The demographic kernels for Plateau Yellow Miner's Candle with the corresponding data from Salguero et al. 2012."----
# run the R file which fits GLMs to the Salguero et al. (2012) data
source("ipm.R")

# Establish a grid based on the range of sizes in the data
alpha <- min(cryptantha$size, cryptantha$sizeNext, na.rm = TRUE)
beta <- max(cryptantha$size, cryptantha$sizeNext, na.rm = TRUE)
xs <- seq(alpha, beta, length = 100)

library(viridis)
# Plot the kernels against the data
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 1, 1), cex.lab = 1.5, cex.axis = 1.5)
# the survival kernel
plot(cryptantha$size, jitter(cryptantha$surv, factor = 0.1), pch = 21,
     bg = rgb(1, 0, 0, 0.5), xlab = "size x", ylab = "survival probability")
lines(xs, Survival(xs), lwd = 4, col = "blue")
# the growth kernel
image(xs, xs, t(outer(xs, xs, Growth)), xlab = expression(size[t]),
      ylab = expression(size[t + 1]), col = inferno(24)[8:24], useRaster=TRUE)
points(cryptantha$size, cryptantha$sizeNext, pch = 21,
       bg = rgb(0.2, 0.5, 1, 0.9))
# the flowering kernel
plot(cryptantha$size, jitter(cryptantha$fec0, factor = 0.1),
     pch = 21, bg = rgb(1, 0, 0, 0.5), xlab = "size", ylab = "flowering")
lines(xs, Flowering(xs), lwd = 4, col = "blue")
# the fecundity data
plot(cryptantha$size[cryptantha$fec0 == 1],
     cryptantha$fec1[cryptantha$fec0 == 1], pch = 21, bg = rgb(1, 0, 0, 0.5),
     bty = "n", xlab = "size", ylab = "fecundity")
lines(xs, Fecundity_flowered(xs), lwd=2)

## ----main.IBIPM.code, cache=FALSE, echo=FALSE----------------------------
#---------
# Creating the discretized probability generating functional
#---------

# the discretization choices
n <- 500  # number of bins
xs <- seq(alpha, 2 * beta, length = n) # the x values
dx <- xs[2] - xs[1] # delta x width

# discretizing the survival, offspring size, and flowering functions
s.vec <- Survival(xs)
e.vec <- dgamma(xs, shape = coef(recruit.size)[1], rate = coef(recruit.size)[2]) * dx
e.vec <- e.vec/sum(e.vec) # normalizing to add up to one
flowering.vec <- as.numeric(Flowering(xs))

# discretizing the growth kernel
# need to make sure the columns add
# to one, i.e. no eviction
g.mat <- outer(xs, xs, Growth) * dx
q.evict <- 1 - colSums(g.mat)
g.mat <- g.mat %*% diag(1/(1 - q.evict))
# Note: most of the eviction is from negative-sized individuals
s.vec <- (1 - q.evict) * s.vec


# defining the probability generating function
# for size-dependent offspring number distribution
phi <- function(x, s) exp(Fecundity(x) * (s - 1))


# putting all the pieces together to create
# the discretized pgf `Phi` which takes
# in discretized functions (i.e. vectors)
# taking values beween [0, 1]
# and returns the same.

Psi <- function(h) {
    h.new <- (s.vec * (t(g.mat) %*% h) + 1 - s.vec) * (1 - flowering.vec + flowering.vec * phi(xs, e.vec %*% h))
    return(h.new)
}


# The following function iterates
# `Psi` for `Tf` time steps on the
# zero vector which yields
# size-dependent extinction probabilities
# in a matrix where rows correspond
# to time and columns correspond to size

Psi.iterate <- function(Tf) {
    hs <- matrix(0, Tf + 1, n)
    for (t in 1:Tf) {
        hs[t + 1, ] <- Psi(hs[t, ])
    }
    return(hs)
}

# A function for estimating asymptotic extinction probability (`AE`)
# by iterating the pgf `Psi` on the zero function until
# changes between iterates have an L1 norm of less than the
# specified tolerance `tol`

AE <- function(tol = 1e-04) {
    h <- rep(0, n)
    diff <- 1
    while (diff > tol) {
        h.old <- h
        h <- Psi(h)
        diff <- sum(abs(h.old - h))
    }
    return(h)
}

# Inputs: `N` vector of length $n$ which specifies the
# initial number of individuals in each size bin, and `Tf`, the
# length of the time to look things over.
# Output: a vector of length T+1 corresponding to the probabilities
# of extinction by time `Tf` or sooner
PVA <- function(N, Tf) {
    hs <- Psi.iterate(Tf)
    output <- numeric(Tf + 1)
    for (t in 1:(Tf + 1)) {
        output[t] <- prod(hs[t, ]^N)
    }
    return(output)
}

## ----extinction, cache=FALSE, fig.width=7, fig.height=4, echo=FALSE, fig.cap="Extinction probabilities for a population initially consisting of a single individual of size $x$. Different curves correspond to extinction occuring in $1, 5, 10, 15$ or $20$ years. Asymptotic extinction probabilities are shown by the thick gray curve. In A, the extinction curves for the baseline individual-based IPM. In B, extinction curves are shown when seed surival is increased by a factor of three."----

# Using the Psi.iterate function to get size- and year-dependent
# extinction probabilities for populations starting with a single
# individual

Tf <- 100
out <- Psi.iterate(Tf)
times <- c(1, 5, 10, 15, 20) + 1

par(cex.lab = 1.25, cex.axis = 1.25, mar = c(4.5, 4.5, 1, 1))
layout(matrix(c(1, 2), 1, byrow = TRUE), widths = c(4, 4))
matplot(xs, t(out[times, ]), xlim = c(alpha, beta), type = "l", lty = 1, xlab = "size x", ylab = "extinction probabilities",
    bty = "n", lwd = 2, col=rev(viridis(7)[2:6]))
h <- AE()
lines(xs, h, type = "l", lwd = 5, xlim = c(alpha, beta), col = viridis(7)[1])

# The following code changes the mean fecundity
# by a factor of 3,  defines the
# appropriate new `Phi` pgf, and computes
# size and year dependent extintion probabilities
# for comparative purposes.

phi.plus = function(x, s) {
  exp(3 * Fecundity(x) * (s - 1))
}

Psi.plus = function(h) {
    h.new <- (s.vec * (t(g.mat) %*% h) + 1 - s.vec) *
      (1 - flowering.vec + flowering.vec * phi.plus(xs, e.vec %*% h))
    return(h.new)
}

Psi.plus.iterate <- function(Tf){
  hs <- matrix(0, Tf+1, n)
  for(t in 1:Tf){
    hs[t + 1, ]=Psi.plus(hs[t, ])
  }
  return(hs)
}

# function for estimating asymptotic extinction probability
AE.plus <- function(tol = 0.001) {
    h = rep(0, n)
    diff = 1
    while (diff > tol) {
        h.old <- h
        h <- Psi.plus(h)
        diff <- sum(abs(h.old - h))
    }
    return(h)
}

out <- Psi.plus.iterate(Tf)
times <- c(1, 5, 10, 15, 20) + 1
matplot(xs, t(out[times, ]), xlim = c(alpha, beta), ylim = c(0, 1), type = "l",
        lty = 1, xlab = "size x", ylab = "extinction probabilities", bty = "n",
        lwd = 2, col=rev(viridis(7)[2:6]))
h <- AE.plus()
lines(xs, h, type = "l", lwd = 5, xlim = c(alpha, beta), col = viridis(7)[1])
par(mar = c(0, 0, 0, 0))
legend("topright", c(as.character(times - 1), expression(Inf)),
       col = c(rev(viridis(7)[1:6])), lty = 1,
       lwd = c(rep(2, length(times)), 5), bty = "n")

## ----compare, fig.width=6, fig.height=4, echo=FALSE, fig.cap="Extinction probabilities as a function of time for populations with $100$ individuals of size $x=1$, and populations with $5$ or $8$ individuals of size $x=60$."----
Tf <- 40
N <- rep(0, n)
N[1] <- 100
out1 <- PVA(N, Tf)
N <- rep(0, n)
N[136] <- 8
out2 <- PVA(N, Tf)
N <- rep(0, n)
N[136] <- 5
out3 <- PVA(N, Tf)
out <- cbind(out1, out2, out3)

par(cex.lab = 1.25, cex.axis = 1.25)
matplot(0:Tf, out, typ = "l", lwd = 3, bty = "n", xlab = "time t",
        ylab = "probability of extinction by year t", col = viridis(4)[1:3])
legend("topleft", c("100 small", "8 large", "5 large"), lty = 1:3, col = viridis(4)[1:3],
       lwd = 3, bty = "n", seg.len = 6)

## ----sampling, cache=TRUE, fig.width=7, fig.height=4, echo=FALSE, fig.cap="Extinction probabilities for founding populations of different abundance. For each founding population abundance $N$, $500$ samples consisting of $N$ randomly choosen individuals from the data set were used to create a founding population of $N$ individuals. Extinction probablities by year $5$ and year $10$ were calculated for each of these sample populations. In A, extinction probability is plotted as a box plot for different $N$ values. In B, extinction probability is plotted against the mean size of an individual for a population abundance $N=25$."----
sizes <- cryptantha$size[which(cryptantha$size > 0 & cryptantha$size < beta)]  # empirical data
T1 <- 5
T2 <- 10
sample.sizes <- 5 * (1:5)
no.samples <- 500
stuff <- c()
mean.size <- c()
stuff2 <- c()
temp <- numeric(no.samples)
temp2 <- temp
temp3 <- temp
for (j in 1:length(sample.sizes)) {
    for (i in 1:no.samples) {
        sampling <- sample(sizes, size = sample.sizes[j])
        out2 <- hist(sampling, breaks = xs, plot = FALSE)
        N <- c(0, out2$counts)
        out <- PVA(N, T2)
        temp[i] <- out[T1 + 1]
        temp2[i] <- out[T2 + 1]
        temp3[i] <- mean(sampling)
    }
    stuff <- cbind(stuff, temp)
    stuff2 <- cbind(stuff2, temp2)
    mean.size <- cbind(mean.size, temp3)
}
par(cex.lab = 1.25, cex.axis = 1.25, mfrow = c(1, 2), mar = c(4.5, 4.5, 1, 1))
boxplot(log10(stuff), notch = FALSE, names = sample.sizes, col = "gray60", range = 0,
    xlab = "initial population abundance", ylab = expression(paste(log[10], " extinction probability")),
    ylim = c(min(log10(c(stuff, stuff2))), max(c(log10(stuff), log10(stuff2)))))
boxplot(log10(stuff2), notch = FALSE, names = sample.sizes, col = "gray90", range = 0,
    add = TRUE)
legend("bottomleft", c("in 10 years", "in 5 years"), pch = 22, pt.bg = c("gray90",
    "gray60"), cex = 1.25, bty = "n")
legend("topright", "A", bty = "n", cex = 1.25)
par(mar = c(4.5, 2.25, 1, 2.25))
i <- 5
plot(mean.size[, i], log10(stuff[, i]), ylim = c(min(log10(c(stuff[, i], stuff2[,
    i]))), max(log10(c(stuff2[, i], stuff[, i])))), xlab = "mean size of an individual",
    ylab = "", pch = 21, bg = "gray60", col="gray20")
points(mean.size[, i], log10(stuff2[, i]), pch = 21, bg = "gray90", col="gray20")
legend("bottomleft", legend=c("in 10 years", "in 5 years"), pch = c(21, 21), pt.bg = c("gray90",
    "gray60"), cex = 1.25, bty = "n")
legend("topright", "B", bty = "n", cex = 1.25)

