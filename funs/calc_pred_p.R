## compute model predictions at subject and group levels.

source(paste0(getwd(),'/funs/setup_for_calc.R'))

# Calculate predictions -----------------------------
plot_grain <- 200

x <- seq(log10(10),log10(200),length=plot_grain) # blur in log10(highcut) units.

# centre as in design matrix:
x_centred <- x - mean(log10(dat$highCut))

# block factors:
blocks <- c(-1, 0, 1)

subject_pred_p <- expand.grid(subject=unique(dat$subject),highCut=10^x,block=c('1','2','3'))
subject_pred_p$y <- NA
subject_pred_p$ymin <- NA
subject_pred_p$ymax <- NA

for (subj in 1:length(unique(dat$subject))){
  thisSubj <- unique(dat$subject)[subj]
  b1 <- params$beta[,subj,1]
  b2 <- params$beta[,subj,2]
  b3 <- params$beta[,subj,3]
  for (j in 1:length(blocks)){
    this_index <- subject_pred_p$subject==thisSubj & subject_pred_p$block==levels(subject_pred_p$block)[j]
    eta <- sapply(x_centred,linPred,blocks[j],b1,b2,b3)
    p <- gamma + (1 - gamma) * p_pred(eta)
    # p is a matrix with rows = samples and columns = y(x). summarise into mean and 95% quantiles.
    p_mean <- colMeans(p)
    p_hdi <- col_hdi(p)
    
    subject_pred_p$y[this_index] <- p_mean
    subject_pred_p$ymin[this_index] <- p_hdi[1,]
    subject_pred_p$ymax[this_index] <- p_hdi[2,]
  }
}

# group-level predictions:
group_pred_p <- expand.grid(group=unique(dat$group4),highCut=10^x,block=c('1','2','3'))
group_pred_p$y <- NA
group_pred_p$ymin <- NA
group_pred_p$ymax <- NA

for (i in 1:length(unique(dat$group4))){
  thisGroup <- unique(dat$group4)[i]
  b1 <- params$mu[,i,1]
  b2 <- params$mu[,i,2]
  b3 <- params$mu[,i,3]
  for (j in 1:length(blocks)){
    this_index <- group_pred_p$group==thisGroup & group_pred_p$block==levels(group_pred_p$block)[j]
    eta <- sapply(x_centred,linPred,blocks[j],b1,b2,b3)
    p <- gamma + (1 - gamma) * p_pred(eta)
    # p is a matrix with rows = samples and columns = y(x). summarise into mean and 95% quantiles.
    p_mean <- colMeans(p)
    p_hdi <- col_hdi(p)
    group_pred_p$y[this_index] <- p_mean
    group_pred_p$ymin[this_index] <- p_hdi[1,]
    group_pred_p$ymax[this_index] <- p_hdi[2,]
  }
}

# Do plots -------------------------------------------------------

#------------------------
# arrange some data:
sum_subject <- ddply(dat,.(subject,group4,highCut,block),summarise,y=mean(correct),ntrials=length(correct),ncorrect=sum(correct))

# beta error bars (with laplacian rule of succession correction):
sum_subject$betaUpper <- qbeta(.975,sum_subject$ncorrect+1,(sum_subject$ntrials - sum_subject$ncorrect)+1)
sum_subject$betaLower <- qbeta(.025,sum_subject$ncorrect+1,(sum_subject$ntrials - sum_subject$ncorrect)+1)

# rename subject codes to match new labels (1--12 excluding dropped cases):
new_subjs <- c("C01", "C02", "C03", "C04", "C05", "C06", "C07", "C08", "C09", "C10", "C11", "C12", "P01", "P02", "P03", "P04", "P05", "P06", "P07", "P08", "P09", "P10", "P11", "P12")

levels(sum_subject$subject) <- new_subjs
levels(subject_pred_p$subject) <- new_subjs

# reorder subject levels to make it easier to arrange by group:
sum_subject$subject <- reorder(sum_subject$subject,as.numeric(sum_subject$group))


#---------------------
# plot individual subjects along with predictions:

x_breaks <- c(16, 32, 186)
fig <- ggplot(sum_subject,aes(x=highCut,y=y,colour=block)) + facet_wrap(~ subject,ncol=4)
fig <- fig + geom_ribbon(data=subject_pred_p,aes(ymin=ymin,ymax=ymax,colour=NULL,fill=block),size=1,alpha=.3)
fig <- fig + geom_line(data=subject_pred_p,size=0.8,size=0.8)
fig <- fig + geom_linerange(aes(y=y,ymax=betaUpper,ymin=betaLower))
fig <- fig + geom_point(aes(y=y),size=2.5)
fig <- fig + scale_x_log10(name="High frequency cutoff (cy / deg)",limits=c(10,200),breaks=x_breaks,labels=round(x_breaks/4))
fig <- fig + scale_y_continuous(breaks=seq(0,1,length=3),name="Proportion correct",limits=c(0,1))
fig <- plcol(fig,name='Block')
fig <- fig + theme_grey(base_size=11)
fig <- fig + theme(legend.position='top')
# fig <- fig + theme(panel.grid.minor=element_line(size=0))
# fig

ggsave(paste0(getwd(),"/figs/performance_by_subject_raw.pdf"),width=6,height=6)
ggsave(paste0(getwd(),"/figs/performance_by_subject_raw.svg"),width=6,height=6)

#---------------------
# predictions on group level, averaging over block...
# group mean and sd performance averaged over subjects:
sum_group <- ddply(sum_subject,.(group4,highCut),summarise,ymean=mean(y),sem=sd(y)/sqrt(length(y)))

# rename all the "group4"s so they're all the same:
dat$group <- dat$group4
sum_group$group <- sum_group$group4
sum_subject$group <- sum_subject$group4

# average model predictions over blocks:
group_pred_p_block_av <- ddply(group_pred_p,.(group,highCut),summarise,y=mean(y),ymin=mean(ymin),ymax=mean(ymax))

# average subject data over blocks:
sum_subj_av_over_blocks <- ddply(sum_subject,.(subject,highCut,group),summarise,y=mean(y))

fig <- ggplot(group_pred_p_block_av,aes(x=highCut,y=y)) + facet_wrap(~ group,ncol=3)
fig <- fig + geom_ribbon(aes(ymin=ymin,ymax=ymax),size=1,colour=NA,fill="black",alpha=.3)
fig <- fig + geom_line(size=0.8,colour="black")
# fig <- fig + geom_point(data=sum_subject,size=1.5,alpha=.5,shape=3)
# plot individual lines joining each subject's points. A bit stupid that I have to loop it:
for (subj in 1:length(levels(sum_subj_av_over_blocks$subject))){
  s_dat <- subset(sum_subj_av_over_blocks,subject==levels(sum_subj_av_over_blocks$subject)[subj])
  fig <- fig + geom_line(data=s_dat,alpha=.3)
}
fig <- fig + geom_linerange(data=sum_group,aes(y=ymean,ymax=ymean+sem,ymin=ymean-sem))
fig <- fig + geom_point(data=sum_group,aes(y=ymean),size=2.5,colour="black",shape=16)
fig <- fig + scale_x_log10(name="High frequency cutoff (cy / deg)",limits=c(10,200),breaks=x_breaks,labels=round(x_breaks/4))
fig <- fig + scale_y_continuous(breaks=seq(0,1,length=5),name="Proportion correct",limits=c(0,1))
fig <- fig + theme_grey(base_size=11)
fig <- fig + theme(panel.grid.minor=element_line(size=0))
# fig

ggsave(paste0(getwd(),'/figs/performance_by_group.pdf'),width=3.5,height=2.5)

#--------------------
# scatter plot for intercept and slope samples?
# take a random sample from the parameter samples so plots aren't too big:
n_plot_points <- nSamples
sample_samples <- sample(1:nSamples,size=n_plot_points)

groupFrame <- expand.grid(group=unique(dat$group4),b1=rep(NA,length=n_plot_points))
groupFrame$b2 <- rep(NA,length=n_plot_points)
groupFrame$b3 <- rep(NA,length=n_plot_points)

param_summaries <- expand.grid(group=unique(dat$group4))

for (group in 1:length(unique(dat$group4))){
  thisGroup <- unique(dat$group4)[group]
  b1 <- params$mu[sample_samples,group,1]
  b2 <- params$mu[sample_samples,group,2]
  b3 <- params$mu[sample_samples,group,3]
  
  groupFrame$b1[groupFrame$group==thisGroup] <- b1
  groupFrame$b2[groupFrame$group==thisGroup] <- b2
  groupFrame$b3[groupFrame$group==thisGroup] <- b3
  
  #   # also summarise parameter means and credible intervals:
  #   param_summaries$x[param_summaries$group==thisGroup] <- quantile(params$mu[,group,1],probs=0.5)
  #   param_summaries$xmin[param_summaries$group==thisGroup] <- quantile(params$mu[,group,1],probs=0.025)
  #   param_summaries$xmax[param_summaries$group==thisGroup] <- quantile(params$mu[,group,1],probs=0.975)
  #   param_summaries$y[param_summaries$group==thisGroup] <- quantile(params$mu[,group,2],probs=0.5)
  #   param_summaries$ymin[param_summaries$group==thisGroup] <- quantile(params$mu[,group,2],probs=0.025)
  #   param_summaries$ymax[param_summaries$group==thisGroup] <- quantile(params$mu[,group,2],probs=0.975)
  
  # also summarise parameter means and credible intervals:
  param_summaries$x[param_summaries$group==thisGroup] <- mean(params$mu[,group,1])
  param_summaries$xmin[param_summaries$group==thisGroup] <- hdi(params$mu[,group,1])[1]
  param_summaries$xmax[param_summaries$group==thisGroup] <- hdi(params$mu[,group,1])[2]
  param_summaries$y[param_summaries$group==thisGroup] <- mean(params$mu[,group,2])
  param_summaries$ymin[param_summaries$group==thisGroup] <- hdi(params$mu[,group,2])[1]
  param_summaries$ymax[param_summaries$group==thisGroup] <- hdi(params$mu[,group,2])[2]
  
  param_summaries$y2[param_summaries$group==thisGroup] <- mean(params$mu[,group,3])
  param_summaries$ymin2[param_summaries$group==thisGroup] <- hdi(params$mu[,group,3])[1]
  param_summaries$ymax2[param_summaries$group==thisGroup] <- hdi(params$mu[,group,3])[2]
}


# scatterplots:
fig <- ggplot(param_summaries,aes(x=x,y=y,shape=group))
fig <- fig + geom_point(data=groupFrame,aes(x=b1,y=b2,colour=group),size=0.8,alpha=0.75)
fig <- fig + geom_point(size=2,aes(fill=group))
fig <- fig + geom_errorbar(aes(ymin=ymin,ymax=ymax))
fig <- fig + geom_errorbarh(aes(xmin=xmin,xmax=xmax))
# fix legend:
fig <- plcol(fig,name='')
fig <- fig + scale_shape_manual(name='',values=c(21,22,23))
fig <- fig + xlab('Intercept') + ylab('Blur slope')
fig <- fig + theme_grey(base_size=11)
fig <- fig + theme(legend.position='top')
# fig

p1 <- fig
# ggsave("coefficient_scatterplot.pdf",width=3.5,height=3.5)

fig <- ggplot(param_summaries,aes(x=x,y=y2,shape=group))
fig <- fig + geom_point(data=groupFrame,aes(x=b1,y=b3,colour=group),size=0.8,alpha=0.75)
fig <- fig + geom_point(size=2,aes(fill=group))
fig <- fig + geom_errorbar(aes(ymin=ymin2,ymax=ymax2))
fig <- fig + geom_errorbarh(aes(xmin=xmin,xmax=xmax))
# fix legend:
fig <- plcol(fig,name='')
fig <- fig + scale_shape_manual(name='',values=c(21,22,23))
fig <- fig + xlab('Intercept') + ylab('Learning slope')
fig <- fig + theme_grey(base_size=11)
fig <- fig + theme(legend.position='top')

p2 <- fig

pdf(paste0(getwd(),'/figs/coefficient_scatterplot.pdf'),width=5,height=3)
multiplot(p1,p2,cols=2)
dev.off()
