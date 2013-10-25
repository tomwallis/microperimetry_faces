## Analyse performance of ideal observer.
# TSAW

rm(list=ls())

library(R.matlab)
library(ggplot2)
library(scales)
library(mgcv)
library(MASS)
library(MPDiR)
library(psyphy)

wd <- getwd()

dat <- readMat(paste(wd,'/ideal_observer/idealPerformance.mat',sep=""))
dat <- dat$dataMat
dat <- data.frame(dat)
names(dat) <- c('blur','nominalNoiseContrast','noiseSD','faceIdent','correct')
dat$blur <- factor(dat$blur)

# remove unused variable levels:
dat$noiseSD <- NULL
dat$faceIdent <- NULL

##------------------
# correct measured ideal observer performance by rule of succession to make it comparable to observer data:
# add one correct and one incorrect trial at each combination of factor levels:

RoS_correct <- expand.grid(nominalNoiseContrast=unique(dat$nominalNoiseContrast),blur=factor(c(1,2,3)),correct=c(0,1))
dat <- rbind(dat,RoS_correct)


# express stimulus dimension as SNR to allow performance to approach chance as SNR -> 0. Signal contrast was .05.
dat$snr <- .05 / dat$nominalNoiseContrast

##------------------
# fit a psychometric function using psyphy toolbox functions for asymptotes.
# Previously tried setting lower asymptote using an offset in glm (see Knoblauch & Maloney, p. 143), but this only works when there is a single predictor (stimulus intensity). Adding a factor ruins the fits. 

glm.probit <- glm(correct ~ log10(snr)*blur, family=binomial(link=mafc.probit(.m=10)),data=dat)
glm.cauchy <- glm(correct ~ log10(snr)*blur, family=binomial(link=mafc.cauchit(.m=10)),data=dat)
glm.weib <- glm(correct ~ log10(snr)*blur, family=binomial(link=mafc.weib(.m=10)),data=dat)
glm.cloglog <- glm(correct ~ log10(snr)*blur, family=binomial(link=mafc.cloglog(.m=10)),data=dat) # this seems to be the same as weibull.

sapply(list(glm.probit,glm.cauchy,glm.weib,glm.cloglog),AIC)

anova(glm.probit,glm.cauchy,glm.weib,glm.cloglog,test="Chisq")

##------------------
# Step through the most complex model using the MASS package to drop less useful terms:
reducedModel <- stepAIC(glm.weib)

##------------------
# examine binomial diagnostics on reduced model. WARNING: this takes a while.
# diags <- binom.diagnostics(reducedModel,nsim=10000)
# plot(diags)

##------------------
# Prediction curves from model fit:
nd <- expand.grid(snr=10^seq(log10(min(dat$snr)),log10(max(dat$snr)),length=300),blur=factor(c(1,2,3)))
pred.pc <- predict(reducedModel,newdata=nd,type="response",se.fit=TRUE)
pred.lin <- predict(reducedModel,newdata=nd,se.fit=TRUE)

nd$pred.pc <- pred.pc$fit
nd$pred.pc.upper <- pred.pc$fit + pred.pc$se.fit * 1.96
nd$pred.pc.lower <- pred.pc$fit - pred.pc$se.fit * 1.96

nd$pred.lin <- pred.lin$fit
nd$pred.lin.upper <- pred.lin$fit + pred.lin$se.fit * 1.96
nd$pred.lin.lower <- pred.lin$fit - pred.lin$se.fit * 1.96


##-------------------
# plot against nominal noise contrast:
p <- ggplot(dat,aes(x=snr,y=correct,shape=blur,linetype=blur)) 
p <- p + geom_ribbon(data=nd,aes(ymin=pred.pc.upper,ymax=pred.pc.lower,y=NULL,colour=NULL),alpha=.25) 
p <- p + geom_line(data=nd,aes(y=pred.pc,x=snr),size=1,alpha=.75)
p <- p + stat_summary(fun.data = "mean_cl_boot",B=3000,size=1,linetype=1)
p <- p + xlab('SNR') + scale_y_continuous(name='Proportion Correct',breaks=seq(0,1,by=.25),limits=c(0,1))
axisbreaks <- trans_breaks("log10", function(x) 10^x)
axislabels <- trans_format("log10", math_format(10^.x))
p <- p + scale_x_log10(breaks = axisbreaks(dat$snr),labels = axislabels(axisbreaks(dat$snr)))
p <- p + geom_vline(aes(xintercept=10),linetype=2,colour="black",size=1)
p <- p + theme_gray(base_size=20)
p <- p + scale_shape_discrete(name="Blur",labels=c("Unblurred","Blur 1","Blur 2"))
p <- p + scale_linetype_discrete(name="Blur",labels=c("Unblurred","Blur 1","Blur 2"))
p <- p + theme(legend.position=c(0.7,0.2))
p

ggsave("performance_ideal.pdf",width=6,height=6)

##-------------------
# save best model fit into a file for future use.
save(reducedModel,dat,nd,file="ideal_observer_fits.RData")
