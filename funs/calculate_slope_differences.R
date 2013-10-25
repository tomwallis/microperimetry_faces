# calculate differences in slope.

source(paste0(getwd(),'/funs/setup_for_calc.R'))

# Compute ----------------------------
#---------------------
# b2:
b2 <- params$mu[,,2]
colnames(b2) <- c('CS','LV','LV:F')

b2_diffs <- array(data=NA,dim=c(nSamples,length(unique(dat$group4))),
                  dimnames=list(sample=NULL,
                                group=c('CS - LV','CS - LV:F','LV - LV:F')))

b2_diffs[,1] <- b2[,1] - b2[,2]
b2_diffs[,2] <- b2[,1] - b2[,3]
b2_diffs[,3] <- b2[,2] - b2[,3]

#---------------------
# b3:
b3 <- params$mu[,,3]
colnames(b3) <- c('CS','LV','LV:F')

b3_diffs <- array(data=NA,dim=c(nSamples,length(unique(dat$group4))),
                  dimnames=list(sample=NULL,
                                group=c('CS - LV','CS - LV:F','LV - LV:F')))

b3_diffs[,1] <- b3[,1] - b3[,2]
b3_diffs[,2] <- b3[,1] - b3[,3]
b3_diffs[,3] <- b3[,2] - b3[,3]


#---------------------
# stick together into one frame:
slope_frame <- expand.grid(comparison=c('CS - LV','CS - LV:F','LV - LV:F'),param=c('Blur slope','Learning slope'))
slope_frame$y <- NA
slope_frame$ymin <- NA
slope_frame$ymax <- NA

for (i in 1: NROW(slope_frame)){
  if(slope_frame$param[i]=='Blur slope'){
    slope_frame$y[i] <- mean(b2_diffs[,i])
    slope_frame$ymin[i] <- hdi(b2_diffs[,i])[1]
    slope_frame$ymax[i] <- hdi(b2_diffs[,i])[2]
  } else {
    slope_frame$y[i] <- mean(b3_diffs[,i-3])
    slope_frame$ymin[i] <- hdi(b3_diffs[,i-3])[1]
    slope_frame$ymax[i] <- hdi(b3_diffs[,i-3])[2]
  }
}

# Plot ------------------------------------------------
fig <- ggplot(slope_frame,aes(x=comparison,y=y)) + facet_wrap(~ param,ncol=2,scales='free_y')
fig <- fig + geom_hline(aes(y_intercept=0))
fig <- fig + geom_errorbar(aes(ymin=ymin,ymax=ymax))
fig <- fig + geom_point()
fig <- plcol(fig,name='')
fig <- fig + scale_shape_manual(name='',values=c(21,22,23))
fig <- fig + xlab('Comparison') + ylab('Slope difference')
fig <- fig + theme_grey(base_size=11)
fig <- fig + theme(legend.position='top')
fig <- fig + theme(axis.text.x = element_text(angle = 45, hjust = 1))
fig

ggsave(paste0(getwd(),'/figs/parameter_difference_plot.pdf'),width=3.5,height=3)

