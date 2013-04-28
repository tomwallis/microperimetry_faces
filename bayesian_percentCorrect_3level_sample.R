## Script to do Bayesian glm on data based on group4 variable. Called "3level" because I use a hierarchical layer to pool across groups (i.e. this fits a model to each subject's data with a pool term for group and a higher pool term for all subjects).
#
# TSAW.

rm(list=ls())

library(ggplot2)
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
  vector<lower=-10,upper=10>[D] mu_mean; 
  vector<lower=0,upper=10>[D] mu_sd; 
  
  vector<lower=0,upper=10>[D] sigma_mean; 
  vector<lower=0,upper=10>[D] sigma_sd; 
  
  matrix[G,D] mu; # group level beta means.
  matrix<lower=0,upper=10>[G,D] sigma; # group level standard deviation on the beta weights.
  matrix[S,D] e_beta; # parameterise as errors in beta using the "matt trick".
}

transformed parameters {
  matrix[S,D] beta;
  vector[D] sigma_alpha;
  vector[D] sigma_beta;

  # convert mu, sigma and e_beta into betas:
  for (d in 1:D){
    sigma_alpha[d] <- pow(sigma_mean[d],2) / pow(sigma_sd[d],2);
    sigma_beta[d] <- sigma_mean[d] / pow(sigma_sd[d],2);

    for (s in 1:S){
      beta[s,d] <- mu[gg[s],d] + e_beta[s,d]*sigma[gg[s],d]; # "matt trick" beta.
    }
  }
}

model {
  vector[N] p; # real-valued p for probability predicted by glm; will be replaced in loop.
  
  mu_mean ~ normal(0,3); 

  # priors on group level mu and sigma:
  for (g in 1:G){
    row(mu,g) ~ normal(mu_mean,mu_sd);
    row(sigma,g) ~ gamma(sigma_alpha,sigma_beta); 
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

#-------------------

# first fit bug test run:
initial_fit <- stan(model_code = stanModel1, data = stanDat,iter = 1, chains = 1, init=0)


# # now proper, with more iterations:
iter <- 100000
warmup <- iter/2
num_chains <- 4
thin <- 25

#prep for parallel case
fun <- function(){
  fit <- stan(fit = initial_fit, data = stanDat, iter = iter, warmup = warmup, chains = 1, verbose = F, init=0,thin=thin)
  return(fit)
}

#parallel case
start <- proc.time()[3]
parallel_fit <- foreach(i = 1:num_chains) %dopar% fun()
proc.time()[3] - start

# squish fit objects together:
fit1 <- sflist2stanfit(parallel_fit)

save(fit1,file='bayesian_percentCorrect_3level_sample.RData')

print(fit1)
plot(fit1)

#---------------
# check sigma params:
params <- extract(fit1)
# #calculate descriptives of population sigma distribution:
# mu <- log(params$sigma_median[,2])
# sigma <- params$sigma_sigma[,2]
# hist(mode <- exp(mu - sigma^2))
# mean <- exp(mu - sigma^2 / 2)
# median <- exp(mu)
# sd <- sqrt( (exp(sigma^2)-1) * exp(2*mu + sigma^2))
# 
# x <- exp(seq(log(0.01),log(20),length=200))
# fun <- function(x,mu,sigma){
#   y <- dlnorm(x,mu,sigma)
# }
# curves <- sapply(x,fun,mu,sigma)
# 
# # create a data frame, sampled from only a few hundred curves:
# curveFrame <- expand.grid(x=x,curve=1:200)
# curveFrame$y <- NA
# ind <- sample(nrow(curves),size=200)
# for (i in 1 : nrow(curves)){
#   curveFrame$y[curveFrame$curve==i] <- curves[ind[i],]
# }
# 
# # plot curves:
# p <- ggplot(curveFrame,aes(x=x,y=y,colour=factor(curve))) + geom_path(alpha=.5)
# p <- p + scale_colour_discrete(guide=FALSE)
# p <- p + scale_x_log10()
# p
