## compute correlations between performance measures and visual sensitivity measures.

rm(list=ls())

library(plyr)
library(ggplot2)
library(boot)
library(xtable)


load('data.Rdata')
source('~/Dropbox/R_Functions/multiplot.R')
source('~/Dropbox/R_Functions/pleasingColours.R')

#----------
# compute summaries for each subject:

sumDat_subject <- ddply(dat,.(subject,group4),summarise,age=mean(age),absolute_scotoma=mean(absolute_scotoma),fixation_stab2=mean(fixation_stab2),pr_better_eye=mean(pr_better_eye),ac_better_eye=mean(ac_better_eye),pc=mean(correct),ntrials=length(correct),ncorrect=sum(correct))

# beta error bars:
sumDat_subject$betaUpper <- qbeta(.975,sumDat_subject$ncorrect,(sumDat_subject$ntrials - sumDat_subject$ncorrect))
sumDat_subject$betaLower <- qbeta(.025,sumDat_subject$ncorrect,(sumDat_subject$ntrials - sumDat_subject$ncorrect))

#-----------
# scatterplots.

# Absolute scotoma against performance:
p1 <- ggplot(sumDat_subject,aes(x=absolute_scotoma,y=pc,fill=group4,shape=group4)) 
p1 <- p1 + geom_linerange(aes(ymax=betaUpper,ymin=betaLower),alpha=.5)
p1 <- p1 + geom_point(size=4,colour="black")
p1 <- p1 + scale_shape_manual(name="",values=c(21,22,23,24))
p1 <- p1 + scale_fill_manual(name="",values=cbPalette)
p1 <- p1 + scale_y_continuous(name="Proportion correct identification",limits=c(0,1),breaks=seq(0,1,by=.25))
p1 <- p1 + scale_x_continuous(name="Proportion of perimetry targets missed",limits=c(0,1),breaks=seq(0,1,by=.25))
p1 <- p1 + theme_grey(base_size=20)
# p1

# fixation stability against performance:
p2 <- ggplot(sumDat_subject,aes(x=fixation_stab2,y=pc,fill=group4,shape=group4)) 
p2 <- p2 + geom_linerange(aes(ymax=betaUpper,ymin=betaLower),alpha=.5)
p2 <- p2 + geom_point(size=4,colour="black")
p2 <- p2 + scale_shape_manual(name="",values=c(21,22,23,24))
p2 <- p2 + scale_fill_manual(name="",values=cbPalette)
p2 <- p2 + scale_y_continuous(name="Proportion correct identification",limits=c(0,1),breaks=seq(0,1,by=.25))
p2 <- p2 + scale_x_continuous(name="Fixation stability (2 degrees)")
p2 <- p2 + theme_grey(base_size=20)
# p2

# PR against performance:
p3 <- ggplot(sumDat_subject,aes(x=pr_better_eye,y=pc,fill=group4,shape=group4)) 
p3 <- p3 + geom_linerange(aes(ymax=betaUpper,ymin=betaLower),alpha=.5)
p3 <- p3 + geom_point(size=4,colour="black")
p3 <- p3 + scale_shape_manual(name="",values=c(21,22,23,24))
p3 <- p3 + scale_fill_manual(name="",values=cbPalette)
p3 <- p3 + scale_y_continuous(name="Proportion correct identification",limits=c(0,1),breaks=seq(0,1,by=.25))
p3 <- p3 + scale_x_continuous(name="Log contrast sensitivity")
p3 <- p3 + theme_grey(base_size=20)
# p3

# AC against performance:
p4 <- ggplot(sumDat_subject,aes(x=ac_better_eye,y=pc,fill=group4,shape=group4)) 
p4 <- p4 + geom_linerange(aes(ymax=betaUpper,ymin=betaLower),alpha=.5)
p4 <- p4 + geom_point(size=4,colour="black")
p4 <- p4 + scale_shape_manual(name="",values=c(21,22,23,24))
p4 <- p4 + scale_fill_manual(name="",values=cbPalette)
p4 <- p4 + scale_y_continuous(name="Proportion correct identification",limits=c(0,1),breaks=seq(0,1,by=.25))
p4 <- p4 + scale_x_continuous(name="LogMAR acuity")
p4 <- p4 + theme_grey(base_size=20)
# p4

svg('correlations_raw.svg',width=14,height=12)
multiplot(p1,p2,p3,p4,cols=2)
dev.off()
# ggsave('correlations.pdf',width=10,height=8)

#-----------
# bootstrapped correlation coefficients:

predictors <- c('targMiss','Fix_2','CS','VA','PC')
varNames <- c('absolute_scotoma','fixation_stab2','pr_better_eye','ac_better_eye','pc')

corFun <- function(dat,inds,variable1,variable2){
  # compute correlation coefficient:
  t <- paste('corr <- cor(dat$',variable1,'[inds], dat$',variable2,'[inds], use= "complete.obs")',sep="")
  eval(parse(text=t))
  return(corr)
}

ind <- 1
correlations <- list()

for (pred1 in 1:(length(predictors)-1)){
  for (pred2 in (pred1+1):length(predictors)){
    var1 <- paste('sumDat_subject$',varNames[pred1],sep='')
    var2 <- paste('sumDat_subject$',varNames[pred2],sep='')
    t <- paste('thisCor <- cor(',var1,', ',var2,', use= "complete.obs")',sep="")
    eval(parse(text=t))
    
    thisBoots <- boot(sumDat_subject,corFun,5000,variable1=varNames[pred1],variable2=varNames[pred2])
    thisCIs <- boot.ci(thisBoots)
    
    # correlation coefficient with test:
    t <- paste('corTest <- cor.test(',var1,', ',var2,', method="pearson")',sep="")
    print(eval(parse(text=t)))
    
    # correlation coefficient with test (spearman):
    t <- paste('corTest <- cor.test(',var1,', ',var2,', method="spearman")',sep="")
    print(eval(parse(text=t)))
    
    correlations[[ind]] <- list(var1=predictors[pred1],var2=predictors[pred2],cor=round(thisCor,digits=2),ci_L=round(thisCIs$bca[4],digits=2),ci_U=round(thisCIs$bca[5],digits=2))
    ind <- ind + 1
  }
}

# #-----------
# # arrange correlations and confidence intervals in a table:
# 

# since we want a pretty custom table, will try to do by hand using paste. Ugh.

# ??? no idea how to get this into xtable.

# Will just cut and paste manually into libreoffice. ugh.
