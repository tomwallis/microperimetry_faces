## Script to do Bayesian glm on data based on group4 variable. Called "2level" because I don't use a hierarchical layer to pool across groups (i.e. this fits a model to each subject's data with a pool term for group).
#
# TSAW. tsawallis@gmail.com

rm(list=ls())

library(rstan)
library(plyr)

# set up for parallel sampling:
library(doMC)
options(cores=4) #set this appropriate to your system
registerDoMC()

#------------------
load("data.Rdata")

# the basic model for each subject is an intercept plus a blur*block interaction:
X <- model.matrix(~ log10(dat$highCut))

gamma <- .1 # chance performance level.

# generate vector associating each subject with a group:
groupDat <- ddply(dat,.(subject),summarise,group=unique(group4))
groupVec <- as.numeric(groupDat$group)

#-----------------
# data formatted for Stan:
stanDat <- list(N = nrow(dat), 
                y = dat$correct,
                D = ncol(X),
                x = X,
                ss = as.numeric(dat$subject),
                S = length(unique(dat$subject)),
                gg = groupVec,
                G = length(unique(dat$group4)),
                gamma = gamma)

# due to the ordering of the trials, the group4 levels correspond to the numerical levels as:
# 1 = control, 2 = LV, 3 = LV:F, 4 = LV:PRL

# stan_rdump output for help from list:
N <- nrow(dat)
y <- dat$correct
D <- ncol(X)
x <- X
ss <- as.numeric(dat$subject)
S <- length(unique(dat$subject))
gg <- groupVec
G <- length(unique(dat$group4))
# stan_rdump(list('N', 'y', 'D', 'x', 'ss', 'S', 'gg', 'G', 'gamma'),file='myDump.Rdump')

#-------------------
# specify model for stan:
stanModel1 <- '
  data {
    int<lower=0> N; # number of data points, across all subjects.
    int<lower=0> S; # number of subjects
    int<lower=1> D; # number of predictor variables (design matrix cols).
    int<lower=1> G; # number of groups.
    int<lower=1,upper=S> ss[N]; # vector associating each trial with a subject.
    int<lower=1,upper=G> gg[S]; # vector associating each subject with a group.
    int<lower=0,upper=1> y[N]; # y value for each trial
    
    real<lower=0,upper=1> gamma; # lower bound of performance (specified in mAFC).    
    matrix[N,D] x; # design matrix for linear model.
  }

parameters {
  matrix<lower=-10,upper=10>[G,D] mu; # group level beta means for each group.
  matrix<lower=0,upper=10>[G,D] sigma; # group level standard deviation on the beta weights.
  matrix[S,D] e_beta; # parameterise as errors in beta using the "matt trick".
  matrix<lower=0,upper=10>[G,D] sigma_mean;
  matrix<lower=0,upper=10>[G,D] sigma_sd;
}

transformed parameters {
   matrix[S,D] beta;
   matrix[G,D] sigma_alpha;
   matrix[G,D] sigma_beta;

   for (d in 1:D){
     for (g in 1:G){
        sigma_alpha[g,d] <- pow(sigma_mean[g,d],2) / pow(sigma_sd[g,d],2);
        sigma_beta[g,d] <- sigma_mean[g,d] / pow(sigma_sd[g,d],2);
     }

     for (s in 1:S){
       beta[s,d] <- mu[gg[s],d] + e_beta[s,d]*sigma[gg[s],d]; # "matt trick" beta.
     }
   }
}

model {
  vector[N] p; # real-valued p for probability predicted by glm; will be replaced in loop.

  # priors on group level mu and sigma:
  for (g in 1:G){
    row(mu,g) ~ normal(0,3);
    row(sigma,g) ~ gamma(row(sigma_alpha,g),row(sigma_beta,g)); 
  }

  # unit normal prior on errors in beta:
  for (s in 1:S) {
    row(e_beta,s) ~ normal(0,1);
  }

  for (n in 1:N) {
    p[n] <- gamma + (1 - gamma) * Phi( dot_product(row(x,n), row(beta,ss[n])) );
  }
  y ~ bernoulli(p); # fitting as a cumulative normal (probit) model.
  
}'

# output model text file:
# write(stanModel1,file='model.txt')

#-------------------

# first fit bug test run:
initial_fit <- stan(model_code = stanModel1, data = stanDat,iter = 1, chains = 1, init=0)


# # now proper, with more iterations:
iter <- 20000
warmup <- iter/2
num_chains <- 4
thin <- 10

#prep for parallel case
fun <- function(){
  fit <- stan(fit = initial_fit, data = stanDat, iter = iter, warmup = warmup, chains = 1, verbose = F, init=0, thin = thin)
  return(fit)
}

#parallel case
start <- proc.time()[3]
parallel_fit <- foreach(i = 1:num_chains) %dopar% fun()
(proc.time()[3] - start)/60/60

# squish fit objects together:
fit1 <- sflist2stanfit(parallel_fit)

save(fit1,file='bayesian_percentCorrect_2level_sample.RData')

print(fit1)
plot(fit1)

