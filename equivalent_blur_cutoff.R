## determine the blur level that would be required to equate normal performance to the LV group.

rm(list=ls())

library(rstan)
library(ggplot2)
library(plyr)

source('~/Dropbox/R_Functions/mcmc_helper_funcs.R')

load('data.Rdata')
load('bayesian_percentCorrect_3level_sample.RData')

# extract permuted mcmc chains:
params <- extract(fit1,permuted=TRUE)
gamma <- 0.1
nSamples <- length(params$mu[,1,1])
# the model is an intercept plus a slope for blur level, for each subject. 

#----------------
# what is the average performance in the unblurred case for LV group?
b1 <- params$mu[,2,1]
b2 <- params$mu[,2,2]
x <- log10(max(unique(dat$highCut)))
gamma <- 0.1
eta <- b1 + b2*x
p <- gamma + (1 - gamma) * pnorm(eta)

target_p <- mean(p)
target_eta <- qnorm((gamma - target_p) / (gamma - 1))

#---------------
# what blur level produces target_p for controls?

b1 <- params$mu[,1,1]
b2 <- params$mu[,1,2]

x <- (target_eta - b1) / b2
(blur_cut <- 10^(mean(x)))

