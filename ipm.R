# This script creates functions for the survival, growth, and fecundity kernels used to create the
# indiviudal-based-integro-projection model.

# Adapted from code of Jessica Metcalf, Sean McMahon, Rob Salguero-Gomez, and Eelke Jongenans from an ESA
# workshop in 2012.  Their code, in turn, is adapted from code by Johan Dahlgren from the publication
# Nicole et al. 2011. Interdependent effects of habitat quality and climate on population growth of an
# endangered plant. J Ecol 99 1211-1218.  The IPM is of the form N[y,t+1]=Integal of s[x]*G[y,x]+F[y,x]dx
# where s[x] is the probability of surviving for an individual of size x G[y,x] is the probability density
# for growing to size y given individual is of size x and has survived F[y,x] is the density of offsrping
# of size y produced by individuals of size x

library(IPMpack)
library(fitdistrplus)

# Read in the data from Salguero-Gomez. See help(dataIPMpackCryptantha) details.
data("dataIPMpackCryptantha")
cryptantha <- dataIPMpackCryptantha

# Establish a grid based on the range of sizes in the data
alpha <- min(cryptantha$size, cryptantha$sizeNext, na.rm = TRUE)
beta <- max(cryptantha$size, cryptantha$sizeNext, na.rm = TRUE)
xs <- seq(alpha, beta, length = 100)

# Estimate the survival kernel S(x)
model.survival <- glm(surv ~ size, data = cryptantha, family = binomial)
Survival <- function(x) {
    predict(model.survival, data.frame(size = x), type = "response")
}

# Estimate the growth kernel G(x,y)
model.growth <- lm(sizeNext ~ size, data = cryptantha)
Growth <- function(y, x) {
    dnorm(y,
          mean = predict(model.growth, data.frame(size = x), type = "response"),
          sd = sd(model.growth$residuals))
}

# Estimate the fecundity kernel, composed of binomial flowering probability and
# poisson fecundity among flowering indiviudals
model.flowering <- glm(fec0 ~ size, data = cryptantha, family = binomial)
Flowering <- function(x) {
    predict(model.flowering, data.frame(size = x), type = "response")
}

model.fecundity <- glm(fec1 ~ size, data = subset(cryptantha, fec0 == 1),
                       family = poisson)
Fecundity_flowered <- function(x) {
    exp(model.fecundity$coefficients[1] + model.fecundity$coefficients[2] * x)
}

# Estimate probability of seed establishment to generate Recruitment function
establishment.prob <- sum(is.na(cryptantha$size)) / sum(cryptantha$fec1, na.rm = TRUE)
Fecundity <- function(x) {
  as.numeric(Fecundity_flowered(x)) * establishment.prob
}


# Estimate mean recruit size and variation around that mean
recruit.size <- fitdist(cryptantha$sizeNext[is.na(cryptantha$size)],"gamma")
