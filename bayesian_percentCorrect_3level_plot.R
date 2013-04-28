## plot results from 3 level model.

rm(list=ls())

library(rstan)
library(ggplot2)
library(plyr)

source('~/Dropbox/R_Functions/hdiOfMCMC.R')

load('data.Rdata')
load('bayesian_percentCorrect_3level_sample.RData')

# extract permuted mcmc chains:
params <- extract(fit1,permuted=TRUE)
gamma <- 0.1
nSamples <- length(params$mu[,1,1])
# the model is an intercept plus a slope for blur level, for each subject. 

# extract unpermuted chains:
# traceplot(fit1)

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

# rename subject codes to match new labels (1--12 excluding dropped cases):
new_subjs <- c("C01", "C02", "C03", "C04", "C05", "C06", "C07", "C08", "C09", "C10", "C11", "C12", "P01", "P02", "P03", "P04", "P05", "P06", "P07", "P08", "P09", "P10", "P11", "P12")

sumDat_subject$subject <- as.character(sumDat_subject$subject)
predFrame$subject <- as.character(predFrame$subject)

for (i in 1:length(levels(dat$subject))){
  old_subj <- as.character(levels(dat$subject)[i])
  sumDat_subject$subject[sumDat_subject$subject==old_subj] <- new_subjs[i]
  predFrame$subject[predFrame$subject==old_subj] <- new_subjs[i]
}
sumDat_subject$subject <- factor(sumDat_subject$subject)
predFrame$subject <- factor(predFrame$subject)

# reorder subject levels to make it easier to arrange by group:
sumDat_subject$subject <- reorder(sumDat_subject$subject,as.numeric(sumDat_subject$group))

predFrame$y <- predFrame$correct

x_breaks <- c(16, 32, 186)

p <- ggplot(sumDat_subject,aes(x=highCut,y=y)) + facet_wrap(~ subject,ncol=4)
p <- p + geom_ribbon(data=predFrame,aes(ymin=p.lower,ymax=p.upper),size=1,colour=NA,fill="black",alpha=.3)
p <- p + geom_line(data=predFrame,size=0.8,colour="black",size=0.8)
p <- p + geom_linerange(aes(y=y,ymax=betaUpper,ymin=betaLower))
p <- p + geom_point(aes(y=y),size=2.5,colour="black")
p <- p + scale_x_log10(name="High frequency cutoff (cy / deg)",limits=c(10,200),breaks=x_breaks,labels=round(x_breaks/4)) 
p <- p + scale_y_continuous(breaks=seq(0,1,length=3),name="Proportion correct",limits=c(0,1))
p <- p + theme_grey(base_size=12)
p <- p + theme(panel.grid.minor=element_line(size=0))
# p

ggsave("bayes_3level_performance_by_subject.pdf",width=6,height=6)
ggsave("bayes_3level_performance_by_subject.svg",width=6,height=6)

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
p3 <- p3 + scale_x_log10(name="High frequency cutoff (cy / deg)",limits=c(10,200),breaks=x_breaks,labels=round(x_breaks/4)) 
p3 <- p3 + scale_y_continuous(breaks=seq(0,1,length=5),name="Proportion correct",limits=c(0,1))
p3 <- p3 + theme_grey(base_size=12)
p3 <- p3 + theme(panel.grid.minor=element_line(size=0))
# p3

ggsave('bayes_3level_performance_by_group.pdf',width=3.5,height=2.5)

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
# a summary frame for mean and median:
sumDistFrame <- ddply(groupFrame,.(group,param),summarise,mean=mean(y),median=median(y))

# exclude intercept; just plot blur slope:
groupFrame <- groupFrame[groupFrame$param=='Blur slope',]
sumDistFrame <- sumDistFrame[sumDistFrame$param=='Blur slope',]

p2 <- ggplot(groupFrame,aes(x=group,y=y))
p2 <- p2 + geom_violin(size=0.8)
# p2 <- p2 + geom_boxplot(notch=TRUE)
p2 <- p2 + geom_point(data=sumDistFrame,aes(y=mean),colour="black",shape=1,size=4)
p2 <- p2 + geom_point(data=sumDistFrame,aes(y=median),colour="black",shape=3,size=4)
p2 <- p2 + theme_grey(base_size=20) + xlab('Group') + ylab('Credible slope values')
p2

ggsave("bayes_3level_params.pdf",width=8,height=6)


#---------------------
# spit out confidence intervals on the difference scores between group parameters:

groupFrame <- expand.grid(group=unique(dat$group4),param=c('Intercept','Blur slope'),y=rep(NA,length=nSamples))

groupN <- ddply(sumDat_subject,.(group4),summarise,N=length(y) / 3)

for (group in 1:length(unique(dat$group4))){
  thisGroup <- unique(dat$group4)[group]
  # how many subjects in this group?:
  (thisN <- groupN$N[groupN$group4==thisGroup])
  
  b1 <- params$mu[,group,1]
  b2 <- params$mu[,group,2]
  
  b1.sd <- params$sigma[,group,1]
  b2.sd <- params$sigma[,group,2]
  
  groupFrame$mu[groupFrame$group==thisGroup & groupFrame$param=='Intercept'] <- b1
  groupFrame$mu[groupFrame$group==thisGroup & groupFrame$param=='Blur slope'] <- b2
  
  groupFrame$sigma[groupFrame$group==thisGroup & groupFrame$param=='Intercept'] <- b1.sd
  groupFrame$sigma[groupFrame$group==thisGroup & groupFrame$param=='Blur slope'] <- b2.sd
  
  groupFrame$N[groupFrame$group==thisGroup] <- rep(thisN,times=nSamples*2)
}

diffFrame <- data.frame()

for (group1 in 1: (length(unique(groupFrame$group)) - 1) ) {
  for (group2 in (group1+1): length(unique(groupFrame$group))) {
    
    thisgroup1 <- unique(groupFrame$group)[group1]
    thisgroup2 <- unique(groupFrame$group)[group2]
    
    for (param in 1:length(unique(groupFrame$param))){
      thisFrame <- data.frame(comparison=paste(group1,' - ',group2,sep=""))
      
      thisParam <- unique(groupFrame$param)[param]
      diffVec <- groupFrame$mu[groupFrame$group==thisgroup1 & groupFrame$param==thisParam] - groupFrame$mu[groupFrame$group==thisgroup2 & groupFrame$param==thisParam]
      diff_mean <- median(diffVec)
      diff_ci <- quantile(diffVec,probs=c(.025,.975))
      diff_hdi <- hdiOfMCMC(diffVec)
      
      # compute effect size as (mu1 - mu2 / sqrt( (sigma1^2 + sigma2^2) / 2))
      sigma1 <- groupFrame$sigma[groupFrame$group==thisgroup1 & groupFrame$param==thisParam]
      sigma2 <- groupFrame$sigma[groupFrame$group==thisgroup2 & groupFrame$param==thisParam]
      
      denomVec <- sqrt( (sigma1^2 + sigma2^2) / 2)
      
      effectSizeVec <- diffVec / denomVec
      effectSize_mean <- median(effectSizeVec)
      effectSize_ci <- quantile(effectSizeVec,probs=c(.025,.975))
      
      # compute weighted effect size (weighted by number of samples in data groups) as: (mu1 - mu2 / sqrt( (sigma1^2*(N1 - 1) + sigma2^2*(N2 - 1)) / (N1 + N2 -2))) 
      group1N <- groupFrame$N[groupFrame$group==thisgroup1][1]
      group2N <- groupFrame$N[groupFrame$group==thisgroup2][1]
      
      denomWeighted <- sqrt( (sigma1^2*(group1N - 1) + sigma2^2*(group2N - 1)) / (group1N + group2N -2))
      
      effectSizeWeighted_Vec <- diffVec / denomWeighted
      effectSizeWeighted_mean <- median(effectSizeWeighted_Vec)
      effectSizeWeighted_ci <- quantile(effectSizeWeighted_Vec,probs=c(.025,.975))
      
      thisFrame$param <- thisParam
      thisFrame$diffCentral <- diff_mean
      thisFrame$diffL <- diff_ci[1]
      thisFrame$diffU <- diff_ci[2]
      thisFrame$hdiL <- diff_hdi[1]
      thisFrame$hdiU <- diff_hdi[2]
      thisFrame$effectSize_Central <- effectSize_mean
      thisFrame$effectSize_L <- effectSize_ci[1]
      thisFrame$effectSize_U <- effectSize_ci[2]
      thisFrame$effectSizeWeighted_Central <- effectSizeWeighted_mean
      thisFrame$effectSizeWeighted_L <- effectSizeWeighted_ci[1]
      thisFrame$effectSizeWeighted_U <- effectSizeWeighted_ci[2]
      
      diffFrame <- rbind(diffFrame,thisFrame)
    }
  }
}

print(diffFrame)

# Examine group differences!

# is group 1 (CS) slope different to zero?
(cs.mid <- median(params$mu[,1,2]))
(cs.ci <- quantile(params$mu[,1,2],probs=c(.025,.975)))

# is group 2 (LV) slope different to zero?
(lv.mid <- median(params$mu[,2,2]))
(lv.ci <- quantile(params$mu[,2,2],probs=c(.025,.975)))

# is group 3 (LV:F) slope different to zero?
(f.mid <- median(params$mu[,3,2]))
(f.ci <- quantile(params$mu[,3,2],probs=c(.025,.975)))

