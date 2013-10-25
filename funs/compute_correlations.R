# compute bayesian correlation coefficients.

source(paste0(getwd(),'/funs/setup_for_calc.R'))

# compute summaries for each subject ----------

sum_subject <- ddply(dat,.(subject,group4),summarise,age=mean(age),absolute_scotoma=mean(absolute_scotoma),fixation_stab2=mean(fixation_stab2),pr_better_eye=mean(pr_better_eye),ac_better_eye=mean(ac_better_eye),pc=mean(correct),ntrials=length(correct),ncorrect=sum(correct))

# beta error bars:
sum_subject$betaUpper <- qbeta(.975,sum_subject$ncorrect,(sum_subject$ntrials - sum_subject$ncorrect))
sum_subject$betaLower <- qbeta(.025,sum_subject$ncorrect,(sum_subject$ntrials - sum_subject$ncorrect))

input_dat <- subset(sum_subject,subset=(group4!='CS'),select=c(absolute_scotoma,fixation_stab2,pr_better_eye,ac_better_eye,pc))


# Compute correlations ----------------------------
fit <- bcor_mcmc(input_dat,prior='lkj_corr(2.0)',iter=100000,n_saved_samples=2000,parallel=TRUE)
save(fit,file=paste0(getwd(),'/output/correlations_samples.RData'))

bcor_plot(input_dat,fit,out_file=paste0(getwd(),'/figs/correlation_histograms.pdf'))
table_list <- bcor_table(input_dat,fit,variable_names = c('TM','FS','CS','VA','PC'))