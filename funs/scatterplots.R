## compute correlations between performance measures and visual sensitivity measures.

source(paste0(getwd(),'/funs/setup_for_calc.R'))

#----------
# compute summaries for each subject:

sum_subject <- ddply(dat,.(subject,group4),summarise,age=mean(age),absolute_scotoma=mean(absolute_scotoma),fixation_stab2=mean(fixation_stab2),pr_better_eye=mean(pr_better_eye),ac_better_eye=mean(ac_better_eye),pc=mean(correct),ntrials=length(correct),ncorrect=sum(correct))

# beta error bars:
sum_subject$betaUpper <- qbeta(.975,sum_subject$ncorrect,(sum_subject$ntrials - sum_subject$ncorrect))
sum_subject$betaLower <- qbeta(.025,sum_subject$ncorrect,(sum_subject$ntrials - sum_subject$ncorrect))

#-----------
# scatterplots.

point_size <- 2.5

# Absolute scotoma against performance:
p1 <- ggplot(sum_subject,aes(x=absolute_scotoma,y=pc,fill=group4,shape=group4)) 
p1 <- p1 + geom_linerange(aes(ymax=betaUpper,ymin=betaLower))
p1 <- p1 + geom_point(size=point_size,colour="black")
p1 <- greyfig(p1,name='',guide=FALSE)
p1 <- p1 + scale_y_continuous(name="Proportion correct",limits=c(0,1),breaks=seq(0,1,by=.25))
p1 <- p1 + scale_x_continuous(name="Perimetry targets missed (prop)",limits=c(0,1),breaks=seq(0,1,by=.25))
p1 <- p1 + theme_grey(base_size=11)
# p1

# fixation stability against performance:
p2 <- ggplot(sum_subject,aes(x=fixation_stab2,y=pc,fill=group4,shape=group4)) 
p2 <- p2 + geom_linerange(aes(ymax=betaUpper,ymin=betaLower))
p2 <- p2 + geom_point(size=point_size,colour="black")
p2 <- greyfig(p2,name='',guide=FALSE)
p2 <- p2 + scale_y_continuous(name="Proportion correct",limits=c(0,1),breaks=seq(0,1,by=.25))
p2 <- p2 + scale_x_continuous(name="Fixation stability (2 degrees)")
p2 <- p2 + theme_grey(base_size=11)
# p2

# PR against performance:
p3 <- ggplot(sum_subject,aes(x=pr_better_eye,y=pc,fill=group4,shape=group4)) 
p3 <- p3 + geom_linerange(aes(ymax=betaUpper,ymin=betaLower))
p3 <- p3 + geom_point(size=point_size,colour="black")
p3 <- greyfig(p3,name='',guide=FALSE)
p3 <- p3 + scale_y_continuous(name="Proportion correct",limits=c(0,1),breaks=seq(0,1,by=.25))
p3 <- p3 + scale_x_continuous(name="Log contrast sensitivity")
p3 <- p3 + theme_grey(base_size=11)
# p3

# AC against performance:
p4 <- ggplot(sum_subject,aes(x=ac_better_eye,y=pc,fill=group4,shape=group4)) 
p4 <- p4 + geom_linerange(aes(ymax=betaUpper,ymin=betaLower))
p4 <- p4 + geom_point(size=point_size,colour="black")
p4 <- greyfig(p4,name='',guide=FALSE)
p4 <- p4 + scale_y_continuous(name="Proportion correct",limits=c(0,1),breaks=seq(0,1,by=.25))
p4 <- p4 + scale_x_continuous(name="LogMAR acuity")
p4 <- p4 + theme_grey(base_size=11)
# p4

pdf(paste0(getwd(),'/figs/correlations.pdf'),width=6,height=6)
multiplot(p1,p2,p3,p4,cols=2)
dev.off()
# ggsave('correlations.pdf',width=10,height=8)

#-----------
# reviewer comment: scatter plot of acuity versus fixation stability. Groups cluster?

fig <- ggplot(sum_subject,aes(x=ac_better_eye,y=fixation_stab2,fill=group4,shape=group4)) 
fig <- fig + geom_point(size=point_size,colour="black")
fig <- greyfig(fig,name='',guide=FALSE)
fig <- fig + scale_y_continuous(name="Fixation stability (2 degrees)")
fig <- fig + scale_x_continuous(name="LogMAR acuity")
fig <- fig + theme_grey(base_size=11)
fig <- fig + theme(legend.position='top')
# fig
ggsave(paste0(getwd(),'/figs/acuity_stability_scatter.pdf'),width=3.5,height=3.5)


