## function to import data files into R.
#
# TSAW

library(R.matlab)
library(stringr)
library(plyr)

subjects <- c('C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12','C13',
              'P01','P02','P03','P04','P05','P06','P07','P09','P10','P11','P12','P13')

# I've excluded P08 as they only did 7 trials in the first block before quitting...

ages <- c(29,20,46,68,40,60,61,34,45,23,34,29,
          52,57,34,26,30,57,70,32,70,69,58,31)

controlAges <- c(29,20,46,68,40,60,61,34,45,23,34,29)
patientAges <- c(52,57,34,26,30,57,70,32,70,69,58,31)

# demographics: age t-test:
t.test(controlAges,patientAges)

#---------------
# visual sensitivity data:
absolute_scotoma <- c(rep(NA,times=length(controlAges)),
                      0.02,0.46,0.67,0.13,0.44,0.00,0.42,0.02,0.83,0.06,0.04,0.02)

fixation_stab2 <- c(rep(NA,times=length(controlAges)),
                    57,9,100,41,42,100,18,41,76,100,72,100)

fixation_stab4 <- c(rep(NA,times=length(controlAges)),
                    94,37,100,99,100,100,100,80,100,100,100,100)

pr_better_eye <- c(rep(NA,times=length(controlAges)),
                   1.50,1.50,1.35,1.65,1.05,1.20,1.05,1.65,0.45,1.35,1.05,1.65)

ac_better_eye <- c(rep(NA,times=length(controlAges)),
                   0.88,0.88,0.18,0.88,1.08,0.18,0.70,0.88,1.00,0.54,0.54,0.00)

#--------------
# patients who received instructions after the first block:
instruct1 <- c('P02','P04','P05')

dataPath <- paste(getwd(),'/data',sep="")

# create an empty data frame to aggregate data:
dat <- data.frame()

for (i in 1 : length(subjects)){
  # use regular expression to match subject data files:
  letter1 <- substr(subjects[i],1,1)
  letter2 <- substr(subjects[i],2,2)
  letter3 <- substr(subjects[i],3,3)
  regexp <- paste("^[",letter1,"][",letter2,"][",letter3,"]",sep="")
  fullPaths <- list.files(path=dataPath,pattern=regexp,full.names=TRUE)
  
  # is this data from a patient or control?
  if(letter1=="C"){
    group1 <- "control"
    group2 <- "control"
  }
  
  if(letter1=="P"){
    group1 <- "patient"
    ifelse(any(subjects[i]==instruct1),{
      group2 <- "instruction"
    },{group2 <- "no instruction"})
  }  
  

  for (j in 1 : length(fullPaths)){
    
    ifelse (subjects[i] == 'C01',{
      # If an RDAta file:
      load(fullPaths[j])
      thisDat <- thisBlock
    },{
      # If a .mat file:
      thisMat <- readMat(fullPaths[j])
      
      thisDat <- thisMat[2]
      thisDat <- data.frame(thisDat$dataMat)
      
    })
    
    names(thisDat) <- c("trialNum","correct","correctFaceGridPos","stimFace","responseFace","contrast","blurLevel")    
    thisDat$subject <- subjects[i]
    thisDat$block <- j
    thisDat$group1 <- group1
    thisDat$group2 <- group2
    thisDat$age <- ages[i]
    thisDat$absolute_scotoma <- absolute_scotoma[i]
    thisDat$fixation_stab2 <- fixation_stab2[i]
    thisDat$fixation_stab4 <- fixation_stab4[i]
    thisDat$pr_better_eye <- pr_better_eye[i]
    thisDat$ac_better_eye <- ac_better_eye[i]
    
    dat <- rbind(dat,thisDat)
  }
  
}


dat$correctFaceGridPos <- factor(dat$correctFaceGridPos)
dat$stimFace <- factor(dat$stimFace)
dat$responseFace <- factor(dat$responseFace)
dat$blurLevel <- factor(dat$blurLevel)
dat$subject <- factor(dat$subject)
dat$block <- factor(dat$block)
dat$group1 <- factor(dat$group1)
dat$group2 <- factor(dat$group2)

# create a continuous variable for blur based on ideal filter cutoffs:
dat$highCut <- NA
#  "unblurred" will use Nyquist limit of 372/2:
dat$highCut[dat$blurLevel==1] <- 186
dat$highCut[dat$blurLevel==2] <- 32
dat$highCut[dat$blurLevel==3] <- 16

# add a continuous covariate for block:
dat$block.continuous <- as.numeric(dat$block)

# #-------------------
# # create a new grouping factor based on perimetry assessments:
# # 'CS' = control subject; no low vision diagnosis.
# # 'LV' = low vision, uses diseased fovea, no PRL.
# # 'LV:PRL' = no functional fovea, uses PRL
# # 'LV:F' = has a functional fovea despite clinical diagnosis.
# dat$group3 <- as.character(dat$group1)
# dat$group3[dat$group3=='control'] <- 'CS'
# dat$group3[dat$subject=='P01'] <- 'LV'
# dat$group3[dat$subject=='P02'] <- 'LV'
# dat$group3[dat$subject=='P03'] <- 'LV:F'
# dat$group3[dat$subject=='P04'] <- 'LV:PRL'
# dat$group3[dat$subject=='P05'] <- 'LV'
# dat$group3[dat$subject=='P06'] <- 'LV:F'
# dat$group3[dat$subject=='P07'] <- 'LV:PRL'
# dat$group3[dat$subject=='P09'] <- 'LV'
# dat$group3[dat$subject=='P10'] <- 'LV'
# dat$group3[dat$subject=='P11'] <- 'LV:F'
# dat$group3[dat$subject=='P12'] <- 'LV'
# dat$group3[dat$subject=='P13'] <- 'LV:F'
# 
# dat$group3 <- factor(dat$group3)

#---- group 4: collapse prl people into LV group.
dat$group4 <- as.character(dat$group1)
dat$group4[dat$group4=='control'] <- 'CS'
dat$group4[dat$subject=='P01'] <- 'LV'
dat$group4[dat$subject=='P02'] <- 'LV'
dat$group4[dat$subject=='P03'] <- 'LV:F'
dat$group4[dat$subject=='P04'] <- 'LV'
dat$group4[dat$subject=='P05'] <- 'LV'
dat$group4[dat$subject=='P06'] <- 'LV:F'
dat$group4[dat$subject=='P07'] <- 'LV'
dat$group4[dat$subject=='P09'] <- 'LV'
dat$group4[dat$subject=='P10'] <- 'LV'
dat$group4[dat$subject=='P11'] <- 'LV:F'
dat$group4[dat$subject=='P12'] <- 'LV'
dat$group4[dat$subject=='P13'] <- 'LV:F'

dat$group4 <- factor(dat$group4)


sumEachGroup <- ddply(dat,.(group4),summarise,nSubjs=length(unique(subject)))

#----------------  
# generate some summary data to check trial / block numbers:
(sumDat <- ddply(dat,.(subject,block,age),summarise,nTrials=length(correct)))

summary(dat)

# save as R data frame:
save(dat,file=paste0(getwd(),"/output/data.Rdata"))

# export as .csv for use in other programs:
write.table(dat,file=paste0(getwd(),"/output/data.txt"),row.names=FALSE)

