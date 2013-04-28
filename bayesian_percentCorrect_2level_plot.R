## plot results from 2 level model.

rm(list=ls())

library(rstan)
library(ggplot2)
library(plyr)

source('~/Dropbox/R_Functions/hdiOfMCMC.R')

load('data.Rdata')
load('bayesian_percentCorrect_2level_sample.RData')

# extract permuted mcmc chains:
params <- extract(fit1,permuted=TRUE)
gamma <- 0.1
nSamples <- length(params$mu[,1,1])
# the model is an intercept plus a slope for blur level, for each subject. 

# first generate mean and 95% credible limits of curve fits to each subject.
# create an x vector for plotting:
x <- seq(log10(10),log10(200),length=200) # blur in log10(highcut) units.
linPred <- function(x,b1,b2){
  eta <- b1 + b2*x
}

thisQuant <- function(p,x,probs=c(0.025,.975)){
  q <- matrix(rep(NA,times=length(x)*length(probs)),nrow=length(probs))
  for (i in 1:length(x)){
    q[,i] <- quantile(p[,i],probs=probs)
  }
  return(q)
}

thisHDI <- function(p,x,credMass=0.95){
  q <- matrix(rep(NA,times=length(x)*2),nrow=2)
  for (i in 1:length(x)){
    q[,i] <- hdiOfMCMC(p[,i],credMass=credMass)
  }
  return(q)
}

predFrame <- expand.grid(subject=unique(dat$subject),highCut=10^x)
predFrame$correct <- NA
predFrame$p.lower <- NA
predFrame$p.upper <- NA

for (subj in 1:length(unique(dat$subject))){
  thisSubj <- unique(dat$subject)[subj]
  b1 <- params$beta[,subj,1]
  b2 <- params$beta[,subj,2]
  eta <- sapply(x,linPred,b1,b2)
  p <- gamma + (1 - gamma) * pnorm(eta)
  # p is a matrix with rows = samples and columns = y(x). summarise into mean and 95% quantiles.
  #   p.mean <- colMeans(p)
  p.mean <- thisQuant(p,x,probs=.5)
  p.quants <- thisQuant(p,x)
  p.hdi <- thisHDI(p,x)
  
  predFrame$correct[predFrame$subject==thisSubj] <- p.mean
  predFrame$p.lower[predFrame$subject==thisSubj] <- p.quants[1,]
  predFrame$p.upper[predFrame$subject==thisSubj] <- p.quants[2,]
}

# group-level predictions:
predFrame.group <- expand.grid(group=unique(dat$group4),highCut=10^x)
predFrame.group$y <- NA
predFrame.group$p.lower <- NA
predFrame.group$p.upper <- NA

for (i in 1:length(unique(dat$group4))){
  thisGroup <- unique(dat$group4)[i]
  b1 <- params$mu[,i,1]
  b2 <- params$mu[,i,2]
  eta <- sapply(x,linPred,b1,b2)
  p <- gamma + (1 - gamma) * pnorm(eta)
  # p is a matrix with rows = samples and columns = y(x). summarise into mean and 95% quantiles.
  #   p.mean <- colMeans(eta)
  p.mean <- thisQuant(p,x,probs=.5)
  p.quants <- thisQuant(p,x)
  p.hdi <- thisHDI(p,x)
  
  predFrame.group$y[predFrame.group$group==thisGroup] <- p.mean
  predFrame.group$p.lower[predFrame.group$group==thisGroup] <- p.quants[1,]
  predFrame.group$p.upper[predFrame.group$group==thisGroup] <- p.quants[2,]
}


#---------------------
# plot individual subjects along with predictions:
sumDat_subject <- ddply(dat,.(subject,group4,highCut),summarise,y=mean(correct),ntrials=length(correct),ncorrect=sum(correct))

# beta error bars:
sumDat_subject$betaUpper <- qbeta(.975,sumDat_subject$ncorrect,(sumDat_subject$ntrials - sumDat_subject$ncorrect))
sumDat_subject$betaLower <- qbeta(.025,sumDat_subject$ncorrect,(sumDat_subject$ntrials - sumDat_subject$ncorrect))

p <- ggplot(dat,aes(x=highCut,y=correct)) + facet_wrap(~ subject,ncol=5)
p <- p + geom_ribbon(data=predFrame,aes(ymin=p.lower,ymax=p.upper),size=1,colour=NA,fill="black",alpha=.3)
p <- p + geom_line(data=predFrame,size=0.8,colour="black",size=0.8)
p <- p + geom_linerange(data=sumDat_subject,aes(y=y,ymax=betaUpper,ymin=betaLower))
p <- p + geom_point(data=sumDat_subject,aes(y=y),size=2.5,colour="black")
p <- p + scale_x_log10(name="High frequency cutoff (cpi)",limits=c(10,200)) 
p <- p + scale_y_continuous(breaks=seq(0,1,length=5),name="Proportion correct",limits=c(0,1))
p <- p + theme_grey(base_size=20)
# p

ggsave("bayes_2level_performance_by_subject.pdf",width=10,height=7.5)

#---------------------
# predictions on group level.

# group mean and sd performance averaged over subjects:
sumDat_group <- ddply(sumDat_subject,.(group4,highCut),summarise,ymean=mean(y),sem=sd(y)/sqrt(length(y)))

# rename all the "group4"s so they're all the same:
dat$group <- dat$group4
sumDat_group$group <- sumDat_group$group4
sumDat_subject$group <- sumDat_subject$group4

p3 <- ggplot(predFrame.group,aes(x=highCut,y=y)) + facet_wrap(~ group,ncol=4)
p3 <- p3 + geom_ribbon(aes(ymin=p.lower,ymax=p.upper),size=1,colour=NA,fill="black",alpha=.3)
p3 <- p3 + geom_line(size=0.8,colour="black")
p3 <- p3 + geom_point(data=sumDat_subject,size=1.5,alpha=.5,shape=3)
p3 <- p3 + geom_linerange(data=sumDat_group,aes(y=ymean,ymax=ymean+sem,ymin=ymean-sem))
p3 <- p3 + geom_point(data=sumDat_group,aes(y=ymean),size=2.5,colour="black",shape=16)
p3 <- p3 + scale_x_log10(name="High frequency cutoff (cpi)",limits=c(10,200)) 
p3 <- p3 + scale_y_continuous(breaks=seq(0,1,length=5),name="Proportion correct",limits=c(0,1))
p3 <- p3 + theme_grey(base_size=20)
p3

ggsave('bayes_2level_performance_by_group.pdf',width=6,height=4)


#---------------------
# arrange and plot group level beta means:

groupFrame <- expand.grid(group=unique(dat$group4),param=c('Intercept','Blur slope'),y=rep(NA,length=nSamples))

for (group in 1:length(unique(dat$group4))){
  thisGroup <- unique(dat$group4)[group]
  b1 <- params$mu[,group,1]
  b2 <- params$mu[,group,2]
  
  groupFrame$y[groupFrame$group==thisGroup & groupFrame$param=='Intercept'] <- b1
  groupFrame$y[groupFrame$group==thisGroup & groupFrame$param=='Blur slope'] <- b2
}

# violin plots by group:
p2 <- ggplot(groupFrame,aes(x=group,y=y)) + facet_wrap(~ param,ncol=2,scales="free_y")
p2 <- p2 + geom_violin(colour=NA,fill="black")
p2 <- p2 + theme_grey(base_size=20) + xlab('Group') + ylab('Credible parameter value')
p2

ggsave("bayes_2level_params.pdf",width=10,height=7.5)

#---------------------
# spit out confidence intervals on the difference scores between group parameters:

diffFrame <- data.frame()

for (group1 in 1: (length(unique(groupFrame$group)) - 1) ) {
  for (group2 in (group1+1): length(unique(groupFrame$group))) {
    
    thisgroup1 <- unique(groupFrame$group)[group1]
    thisgroup2 <- unique(groupFrame$group)[group2]
    
    for (param in 1:length(unique(groupFrame$param))){
      thisFrame <- data.frame(comparison=paste(group1,' - ',group2,sep=""))
      
      thisParam <- unique(groupFrame$param)[param]
      diffVec <- groupFrame$y[groupFrame$group==thisgroup1 & groupFrame$param==thisParam] - groupFrame$y[groupFrame$group==thisgroup2 & groupFrame$param==thisParam]
      diff_mean <- median(diffVec)
      diff_ci <- quantile(diffVec,probs=c(.025,.975))
      
      thisFrame$param <- thisParam
      thisFrame$diffCentral <- diff_mean
      thisFrame$diffL <- diff_ci[1]
      thisFrame$diffU <- diff_ci[2]
      
      diffFrame <- rbind(diffFrame,thisFrame)
    }
  }
}

print(diffFrame)

# examine group differences!
