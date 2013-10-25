## calculate mean performance differences at each blur and block.

source(paste0(getwd(),'/funs/setup_for_calc.R'))

# Compute mean performance differences -----------------
x <- c(log10(16), log10(32), log10(186))

# centre as in design matrix:
x <- x - mean(log10(dat$highCut))

# block factors:
blocks <- c(-1, 0, 1)

group_samples <- array(data=NA,dim=c(nSamples,length(x),length(unique(dat$group4)),length(blocks)),
                       dimnames=list(sample=NULL,
                                     blur=c('4','8','46'),
                                     group=c('CS','LV','LV:F'),
                                     block=c('1','2','3')))

# at the group level, calculate the predicted proportion correct for each mcmc sample at each blur level and block.
for (i in 1:length(unique(dat$group4))){
  thisGroup <- unique(dat$group4)[i]
  b1 <- params$mu[,i,1]
  b2 <- params$mu[,i,2]
  b3 <- params$mu[,i,3]
  #   b4 <- params$mu[,i,4]
  for (j in 1:length(blocks)){
    eta <- sapply(x,linPred,blocks[j],b1,b2,b3)
    p <- gamma + (1 - gamma) * p_pred(eta)
    # p is a matrix with rows = samples and columns = y(x).
    
    # save groups to array:
    group_samples[,,i,j] <- p
  }
  
}

colMeans(group_samples)

# now take difference score matrices between the three groups:
diff_array <- array(data=NA,dim=c(nSamples,length(x),length(unique(dat$group4)),length(blocks)),
                    dimnames=list(sample=NULL,
                                  blur=c('4','8','46'),
                                  group=c('CS - LV','CS - LV:F','LV - LV:F'),
                                  block=c('1','2','3')))
for (j in 1:length(blocks)){
  diff_array[,,1,j] <- group_samples[,,1,j] - group_samples[,,2,j]
  diff_array[,,2,j] <- group_samples[,,1,j] - group_samples[,,3,j]
  diff_array[,,3,j] <- group_samples[,,2,j] - group_samples[,,3,j]
}


colMeans(diff_array)

# now summarise quantiles of the difference scores:
diff_frame <- expand.grid(comparison=c('CS - LV','CS - LV:F','LV - LV:F'),blur_level=c('4','8','46'),block=c('1','2','3'))
diff_frame$y <- NA
diff_frame$ymin <- NA
diff_frame$ymax <- NA

for (i in 1:length(unique(dat$group4))){
  for (j in 1:length(unique(diff_frame$blur_level))){
    this_blur <- levels(diff_frame$blur_level)[j]
    for (k in 1:length(blocks)){
      this_block <- levels(diff_frame$block)[k]
      this_index <- diff_frame$blur_level==this_blur & diff_frame$comparison==levels(diff_frame$comparison)[i] & diff_frame$block==this_block
      # take median and 95% quantile range:
      #       this_quantiles <- col_quantile(diff_array[,,i,k],probs=c(0.025,0.5,0.975)) # probabilities in rows, blur level in columns.
      #       diff_frame$ymin[this_index] <- this_quantiles[1,j]
      #       diff_frame$y[this_index] <- this_quantiles[2,j]
      #       diff_frame$ymax[this_index] <- this_quantiles[3,j]
      
      # take mean and HDI:
      this_hdi <- col_hdi(diff_array[,,i,k])
      this_mean <- colMeans(diff_array[,,i,k])
      diff_frame$ymin[this_index] <- this_hdi[1,j]
      diff_frame$y[this_index] <- this_mean[j]
      diff_frame$ymax[this_index] <- this_hdi[2,j]
    }
    
  }
}

# rename columns:
colnames(diff_frame) <- c('comparison','blur','block','mean','minHDI','maxHDI')

# Plot mean performance differences ------------------------------

levels(diff_frame$block) <- c('Block 1','Block 2','Block 3')

fig <- ggplot(diff_frame,aes(x=comparison,y=mean,shape=blur,colour=blur,fill=blur)) + facet_wrap(~ block, ncol=3)
fig <- fig + geom_hline(aes(y_intercept=0))
fig <- fig + geom_linerange(aes(ymin=minHDI,ymax=maxHDI),position=position_dodge(width=0.7))
fig <- fig + geom_point(position=position_dodge(width=0.7))
fig <- greyfig(fig,name='Blur')
fig <- fig + xlab('Comparison') + ylab('Performance difference')
fig <- fig + theme_grey(base_size=11)
fig <- fig + theme(legend.position='top')
fig <- fig + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# fig

ggsave(paste0(getwd(),'/figs/mean_performance_difference_plot.pdf'),width=5,height=4)
