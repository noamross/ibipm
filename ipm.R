# This script creates functions for the survival, growth, and fecundity kernels used to create the
# indiviudal-based-integro-projection model.

library(Matrix)
library(nlme)
library(MASS)
library(IPMpack)
library(fitdistrplus)

# Read in the data from Salguero-Gomez. See help(dataIPMpackCryptantha) details.
data("dataIPMpackCryptantha")
cryptantha <- dataIPMpackCryptantha


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
