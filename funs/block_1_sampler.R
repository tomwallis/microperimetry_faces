# script to sample only block 1 data. 

rm(list=ls())

library(rstan)
library(plyr)

#------------------
load("data.Rdata")

# Exclude all but block 1-----------------
dat <- subset(dat,block=='1')

# set up data for stan-----------------

# the basic model for each subject is an intercept plus a blur*block interaction:
X_raw <- model.matrix(~ scale(log10(dat$highCut),scale=FALSE))

gamma <- .1 # chance performance level.

# convert into binomial (ncorrect, ntrials) rather than bernoulli to speed up:
X <- data.frame(X_raw)
betaNames <- colnames(X)
X$y <- dat$correct 
X$subject <- dat$subject

X_binom <- ddply(X,(c(betaNames,'subject')),summarise,yBinom=sum(y),N=length(y))
X_binom_design <- subset(X_binom,select=c(-yBinom,-N,-subject)) # should have same number of columns as X_raw.

# generate vector associating each subject with a group:
groupDat <- ddply(dat,.(subject),summarise,group=unique(group4))
groupVec <- as.numeric(groupDat$group)

stanDat <- list(N = nrow(X_binom), 
                y = X_binom$yBinom,
                n_trials = X_binom$N,
                D = ncol(X_raw),
                x = X_binom_design,
                ss = as.numeric(X_binom$subject),
                S = length(unique(dat$subject)),
                gg = groupVec,
                G = length(unique(dat$group4)),
                gamma = gamma)

# due to the ordering of the trials, the group4 levels correspond to the numerical levels as:
# 1 = control, 2 = LV, 3 = LV:F, 4 = LV:PRL
#-------------------
# specify model for stan:
model <- '
data {
int<lower=0> N; # number of data points, across all subjects.
int<lower=0> S; # number of subjects
int<lower=1> D; # number of predictor variables (design matrix cols).
int<lower=1> G; # number of groups.
int<lower=1,upper=S> ss[N]; # vector associating each trial with a subject.
int<lower=1,upper=G> gg[S]; # vector associating each subject with a group.
int<lower=0> y[N]; # number correct for each cell.
int<lower=0> n_trials[N]; # number of trials for each cell.
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
p[n] <- gamma + (1 - gamma) * inv_logit( dot_product(row(x,n), row(beta,ss[n])) );
}

y ~ binomial(n_trials,p); 

}'

# Do sampling-------------------

source('~/Dropbox/R_Functions/mcmc_helper_funcs.R')
fit1 <- stan_sample(model_string = model, stan_data = stanDat, iter = 100000, n_saved_samples=2000)

save(fit1,file='block_1_samples.RData')

print(fit1)

#----------------
# plot diagnostics using ggmcmc:
library(ggmcmc)
sims <- ggs(fit1,description='blur + block model')
# ggmcmc(sims,file='final_model_convergence_checks_mu.pdf',family = c('mu'))
# ggmcmc(sims,file='final_model_convergence_checks_beta.pdf',family = c('beta'))

r_means <- ggs_running(sims,family='mu')
ggsave('block_1_convergence_running_means.pdf',width=10,height=15)
traces <- ggs_traceplot(sims,family='mu')
ggsave('block_1_convergence_traces.pdf',width=10,height=15)
