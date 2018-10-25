##Analysis Script for the HypnoBlink Experiment
##By subcomarc, with the resources of J?r?me Sackur
## Needs file "blinkData.R"

############### ################################
#Preloading settings, directories and libraries#
############### ################################

#setwd("C:/Users/lscpuser/Google Drive/R/")
#source("C:/Users/lscpuser/Google Drive/R/R_Lib_Sackur.r")
#setwd("C:/Users/subcomarc/Google Drive/R")
#source("C:/Users/subcomarc/Google Drive/R/R_Lib_Sackur.r")
setwd("/home/subcomarc/Google_Drive/R/")
source("/home/subcomarc/Google_Drive/R/R_Lib_Sackur.r")

library("nlme")
library('xtable')
library('knitr')
library('lme4')
library('car')
library('multcomp')
library('lsmeans')
library("Hmisc")
library("HH")


options(xtable.floating = FALSE)
options(xtable.timestamp = "")


##################### ###################################################
#Load preprocessed data (if this is the first time, DO PREPPING INSTEAD)#
#################### ####################################################

#load("blinkData.r")
load("/home/subcomarc/Google_Drive/R/blinkData.r")
#load("/media/subcomarc/201ED7AF1ED77BEA/Documents and Settings/subcomarc/Google Drive/R/blinkData.r")
#data$hypnoC[data$hypnoC==4]<-0 # If you are going to plot, so that no hypno is next to low

####################### ########################################
#Prepping the data for stats and adding some convenient columns#
###################### ######################################### ##################

#data <- read.table(file="C:/Users/lscpuser/Google Drive/R/BlinkTable.txt", header=T)
data <- read.table(file="C:/Users/subcomarc/Google Drive/R/BlinkTable.txt", header=T)

data <- droplevels(data[data$SubName!=101,])
data$suj <- data$SubName
sujs <- unique(data$suj); length(sujs)
ntrial <- data$WhichTrial
nblock <- data$WhichBlock

data$target1 <- 'neutral'
data$target1[data$FaceCombination== 3 | data$FaceCombination == 4] <- "angry"
data$target2 <- 'neutral'
data$target2[data$FaceCombination== 1 | data$FaceCombination == 4 | data$FaceCombination == 5] <- "angry"


data$nface <- 2
data$nface[data$FaceCombination>4] <- 1
data$accn <- data$FacesCorrect
data$accOr <- data$AllAnswers
data$rtOr <- data$RTimeWS*1000
hist(data$rtOr[data$rtOr < 10000])
data$rtn <- data$RTimeHM*1000
hist(data$rtn[data$rtn<20000])

## 3: high: 2: medium; 1: low
data$hypno <- data$HLRAWcol
data$hypnoC <- data$HLcol
#data$hypnoC[data$hypnoC==rm0]<-4
data$lag <- data$AllLags
data$soa <- data$AllSOAs
data$lag <- ordered(data$lag)

## new hypno scores with 8 as highs

data$NewH[data$hypnoC==1]<-1
data$NewH[data$HLRAWcol %in% c(5,6,7)]<-2
data$NewH[data$HLRAWcol %in% c(8,9,10,11,12)]<-3


## outliers

e <- evalq(aggregate(list(accn=accn, accOr=accOr, hypnoC=hypnoC), list(suj=suj), mean), data)
evalq(plot(accn, accOr, col=hypnoC, type='p', xlim=c(0, 1), ylim=c(0,1)), e)
l <- lm(accOr~accn, e)
abline(l) 

e$acc <- (e$accOr+e$accn)/2
hist(e$acc); mean(e$acc); sd(e$acc)

sum(e$acc < .6)
bad <- as.vector(e$suj[e$acc < .6])

e <- e[!e$suj %in% bad,]
min(e$acc); hist(e$acc)

data <- droplevels(data[!data$suj %in% bad,]) ## getting rid of accuracy outliers (3 subs)

sum(data$rtOr > 4000)

hist(data$rtOr[data$rtOr<4000])

e <- evalq(aggregate(list(medianRtOr=rtOr, hypno=hypno, hypnoC=hypnoC), list(suj=suj), median),
           data[data$rtOr < 10000,])
evalq(plot(hypno, medianRtOr, col=hypnoC,type='p'), e)
boxplot(e$medianRtOr)
sd(e$medianRtOr)
limRt <- mean(e$medianRtOr)+3*sd(e$medianRtOr)

bad <- as.vector(e$suj[e$medianRtOr > limRt])
data <- droplevels(data[!data$suj %in% bad,]) ## getting rid of RT outliers (1 sub)

e <- evalq(aggregate(list(medianRtOr=rtOr, hypno=hypno, hypnoC=hypnoC), list(suj=suj), median),
           data)
e$sd <- evalq(aggregate(rtOr, list(suj=suj), sd), data)$x

data <- merge(data, e)

data$outlier <- data$rtOr > data$medianRtOr + 3 * data$sd

e <- evalq(aggregate(list(outlier=outlier), list(suj=suj), mean), data)
any(e$outlier>.1); max(e$outlier)
data <- droplevels(data[!data$outlier & data$rtOr>150,])
hist(data$rtOr[data$rtOr < 300])
data$hypnoC <- ordered(data$hypnoC)
hist(data$rtOr)
data$lagC <- "early"
data$lagC[data$lag>=5] <- "late"

e <- evalq(aggregate(list(medianRtOr=rtOr), list(suj=suj), median),
           data)
f <- evalq(aggregate(list(hypnoC=hypnoC), list(suj=suj), unique),
           data)
e <- merge(e,f)

boxplot(e$medianRtOr)
bad <- as.vector(e$suj[e$medianRtOr > 1200])
data <- droplevels(data[!data$suj %in% bad,])

## Fake lag 0 for no T1

data$lagComplete <- as.numeric(data$lag)
data$lagComplete[data$nface==1] <- 0
data$lagComplete <- ordered(data$lagComplete)

data$lagCatComplete <- "early"
data$lagCatComplete[data$lag>=5] <- "late"
data$lagCatComplete[data$nface==1] <- "a_noT1" # a for alphebetical order for plotting purposes


data$HypnoYN<-1 # 
data$HypnoYN[data$suj>=300]<-0 #

## 3 hyps vs group without hypnosis of all hyps
data$OldComp<-1
data$OldComp[data$suj>=400]<-0

## 3 hyps for everyone, hypnosis vs no hypnosis
data$NewComp<-1
data$NewComp[data$suj>=300]<-0
data$NewComp[data$suj>=400]<-1


#data$newhypnoC<-data$hypnoC # 
#data$newhypnoC[data$hypnoC==4]<-2 # 
#data$newhypnoC[data$suj %in% c(301,302)]<-1 # 
#data$newhypnoC[data$suj %in% c(303,304)]<-3 # 



save(data, file='blinkData.R')




###################### #########################
#All contrasts and dataframes for linear models#
##################### ##########################
data$target1 <- factor(data$target1)
data$target2 <- factor(data$target2)
data$hypnoC <- factor(data$hypnoC)
#data$hypnoC <- ordered(data$hypnoC)
#data$suj <- factor(data$suj)
data$lag <- factor(data$lag)
data$lagC <- factor(data$lagC)
data$lagComplete <- factor(data$lagComplete)
data$lagCatComplete <- factor(data$lagCatComplete)
data$HypnoYN <- factor(data$HypnoYN)
data$NewH <- factor(data$NewH)
data$nface <-factor(data$nface)
#data$newhypnoC <- factor(data$newhypnoC)
#data$HLRAWcol<-ordered(data$HLRAWcol)

contrasts(data$target1)=contr.sum(2)
contrasts(data$target2)=contr.sum(2)
contrasts(data$hypnoC)=contr.sum(4)
contrasts(data$lag)=contr.sum(6)
contrasts(data$lagC)=contr.sum(2)
contrasts(data$lagComplete)=contr.sum(7)
contrasts(data$lagCatComplete)=contr.sum(3)
contrasts(data$HypnoYN)=contr.sum(2)
contrasts(data$NewH)=contr.sum(3)
contrasts(data$nface)=contr.sum(2)
#contrasts(data$newhypnoC)=contr.sum(3)
#contrasts(data$HLRAWcol)=contr.sum(12)

data$HLRAWcolSCALED<-scale(data$HLRAWcol, center = TRUE, scale = TRUE)

###################### ###################### ######################
#STATISTICAL ANALYSES#
###################### ###################### ######################


#Remember: accOr= Accuracy Orientation Task. accn= Accuracy Number of Faces Task.
#rtOr= RT orientation task. rtn= RT number of faces task. hypnoC= hypnosis category (no hypno 4, or 0 
#if you plotting)

#Attentional Blink????????
## #######################

##Full models for acc and RT, plus check for best model

##Acc Anova(AccInc3)

AccInc <- glmer(accOr~lagC+(1|suj), family='binomial', 
                droplevels(data[data$NewH %in% c(1,2,3) & data$NewComp %in% c(1) & 
                                  data$nface==2 & data$HypnoYN %in% c(1),]))

AccInc2 <- glmer(accOr~lagC*target1+(1|suj), family='binomial', 
                droplevels(data[data$NewH %in% c(1,2,3) & data$NewComp %in% c(1) & 
                                  data$nface==2 & data$HypnoYN %in% c(1),]))

AccInc3 <- glmer(accOr~lagC*target1*target2+(1|suj), family='binomial', 
                droplevels(data[data$NewH %in% c(1,2,3) & data$NewComp %in% c(1) & 
                                  data$nface==2 & data$HypnoYN %in% c(1),]))

AccInc4 <- glmer(accOr~target2*NewH+(1|suj), family='binomial', 
                droplevels(data[data$NewH %in% c(1,2,3) & data$NewComp %in% c(1) & 
                                  data$nface==2 & data$HypnoYN %in% c(0),]))

AccInc45 <- glmer(accOr~lagC*target1+(1|suj), family='binomial', 
                 droplevels(data[data$NewH %in% c(1,2,3) & data$NewComp %in% c(1) & 
                                   data$nface==2 & data$HypnoYN %in% c(1),]))

AccInc5 <- glmer(accOr~lagC*target1*target2*HypnoYN+(1|suj), family='binomial', 
                droplevels(data[data$NewH %in% c(1,2,3) & data$NewComp %in% c(1) & 
                                  data$nface==2,]))

AccInc6 <- glmer(accOr~lagC*target1*target2*HypnoYN+(1|suj), family='binomial', 
                 droplevels(data[data$NewH %in% c(1,2,3) & data$NewComp %in% c(1) & 
                                   data$nface==2,]))

AccInc7 <- glmer(accOr~lagC*target1*target2+HypnoYN*NewH+(1|suj), family='binomial', 
                 droplevels(data[data$NewH %in% c(1,2,3) & data$NewComp %in% c(1) & 
                                   data$nface==2,]))

AccInc8 <- glmer(accOr~lagC*target1*target2+HypnoYN+NewH+(1|suj), family='binomial', 
                 droplevels(data[data$NewH %in% c(1,2,3) & data$NewComp %in% c(1) & 
                                   data$nface==2,]))


AllModelsAcc<-anova(AccInc,AccInc2,AccInc3,AccInc4,AccInc5,AccInc6)
ModelsInterHYN<-anova(AccInc6,AccInc5)
ModelsInterNwH<-anova(AccInc4,AccInc5)
ModelsInteractionHypno<-anova(AccInc8,AccInc5)
ModelsInteractionHypno2<-anova(AccInc7,AccInc5)

Anova(AccInc4)

Anova(AccInc4)
lsm.options(pbkrtest.limit = 9557)
lsmeans(AccInc4, pairwise ~ target2*NewH, adjust="tukey")

#Test plot like scatter

e <- evalq(aggregate(list(accOr=accOr), list(suj=suj, lagC=lagC, target1=target1), mean),droplevels(data[data$NewH %in% c(1,2,3) & data$NewComp %in% c(1) & 
                                                                                  data$nface==2 & data$HypnoYN %in% (0),]))

ggplot() +
  #scale_shape_discrete(solid=T, legend=F)
  geom_point(data=e, mapping=aes(x=lagC, y=accOr, color=target1, size="0.5"))


##RT

RTInc <- lmer(rtOr~lagC+(1|suj), 
                droplevels(data[data$NewH %in% c(1,2,3) & data$NewComp %in% c(1) & 
                                  data$nface==2,]))

RTInc2 <- lmer(rtOr~lagC*target1+(1|suj), 
                 droplevels(data[data$NewH %in% c(1,2,3) & data$NewComp %in% c(1) & 
                                   data$nface==2,]))

RTInc3 <- lmer(rtOr~lagC*target1*target2*NewH+(1|suj), 
                 droplevels(data[data$NewH %in% c(1,2,3) & data$NewComp %in% c(1) & 
                                   data$nface==2 & data$HypnoYN %in% c(0),]))

RTInc4 <- lmer(rtOr~lagC*target1*target2*NewH+(1|suj), 
                 droplevels(data[data$NewH %in% c(1,2,3) & data$NewComp %in% c(1) & 
                                   data$nface==2,]))

RTInc5 <- lmer(rtOr~lagC*target1*target2*NewH*HypnoYN+(1|suj), 
                 droplevels(data[data$NewH %in% c(1,2,3) & data$NewComp %in% c(1) & 
                                   data$nface==2,]))

RTInc6 <- lmer(rtOr~lagC*target1*target2*HypnoYN+(1|suj), 
               droplevels(data[data$NewH %in% c(1,2,3) & data$NewComp %in% c(1) & 
                                 data$nface==2,]))

RTInc7 <- lmer(rtOr~lagC*target1*target2+HypnoYN*NewH+(1|suj), 
                 droplevels(data[data$NewH %in% c(1,2,3) & data$NewComp %in% c(1) & 
                                   data$nface==2,]))

RTInc8 <- lmer(rtOr~lagC*target1*target2+HypnoYN+NewH+(1|suj), 
                 droplevels(data[data$NewH %in% c(1,2,3) & data$NewComp %in% c(1) & 
                                   data$nface==2,]))


AllModelsRT<-anova(RTInc,RTInc2,RTInc3,RTInc4,RTInc5,RTInc6)
ModelsInterHYN<-anova(RTInc6,RTInc5)
ModelsInterNwH<-anova(RTInc4,RTInc5)
ModelsInteractionHypno<-anova(RTInc8,RTInc5)
ModelsInteractionHypno2<-anova(RTInc7,RTInc5)


##Whitout Hypnosis
#Accuracy

AccInc <- glmer(accOr~lagC*target1*target2+(1|suj), family='binomial', 
                data=droplevels(data[data$hypnoC %in% c(4),]))




#data$HLRAWcol<-scale(data$HLRAWcol, center = TRUE, scale = TRUE)
AccInc <- glmer(accOr~lagC*target1*target2*HLRAWcol*HypnoYN+(1|suj), family='binomial', 
                data=data)

AccInc <- glmer(accOr~lagC*HypnoYN+lagC:HypnoYN+(1|suj), family='binomial', 
                data=data)

AccInc <- glmer(accOr~lagC*target1*target2*HLRAWcol*HypnoYN+(1|suj), family='binomial', 
                data=data)


Anova(AccInc)
lsm.options(pbkrtest.limit = 9557)
lsmeans(AccInc, pairwise ~ target1*lagC | target1*target2, adjust="tukey")
lsmeans(AccInc, pairwise ~ target1*lagC, adjust="tukey")
lsmeans(AccInc, pairwise ~ target1*target2, adjust="tukey")


e <- evalq(aggregate(list(accOr=accOr), list(suj=suj, target1=target1, lagC=lagC, target2=target2, NewH=NewH), 
                     mean),droplevels(data[data$HypnoYN %in% c(0),]))
evalq(interaction.plot(lagC, target1, accOr), e)
evalq(interaction.plot(lagC, NewH, accOr), e)
evalq(interaction.plot(lagC, HypnoYN, accOr), e)
evalq(interaction.plot(HLRAWcol, HypnoYN, accOr), e)

## Draw SE bars without participant means

grandMean <- mean(e$accOr)
partMeans <-evalq(aggregate(list(partMean=accOr), list(suj=suj), mean),
                 droplevels(data[data$hypnoC %in% c(4),]))
e <- merge(e, partMeans)
e$naccOr <- e$accOr-e$partMean+grandMean
f <- evalq(aggregate(list(seaccOr=naccOr), list(target1=target1,target2=target2), se), e)
f$seaccOr <- f$seaccOr*sqrt(6/3) # compute se and apply Morey's 2008 correction (For a 2X2 design M=4).
f$maccOr <- evalq(aggregate(list(maccOr=accOr), list(target1=target1,target2=target2), mean), e)$maccOr

evalq(interaction.plot(target1, target2, accOr, ylim=c(0.6,0.79)) ,e)
j <- 0.001
xs <- c(1-j,2-j,1+j,2+j)
y1=f$maccOr+f$seaccOr
y0=f$maccOr-f$seaccOr

arrows(x0=xs, y0=y0, x1=xs, y1=y1, length=0,
       angle=90,code=3, lty=1, lwd=2)


#Response Time

RTInc <- lmer(rtOr~lagC*target1*target2+(1|suj), data=droplevels(data[data$hypnoC %in% c(4),]))
qqnorm(residuals(RTInc))
plot(fitted(RTInc),residuals(RTInc))
Anova(RTInc)
lsm.options(pbkrtest.limit = 9557)
lsmeans(RTInc, pairwise ~ target1*target2*lagC, adjust="tukey")

e <- evalq(aggregate(list(accOr=accOr), list(suj=suj, target1=target1, target2=target2, lagC=lagC, NewH=NewH), 
                     mean), droplevels(data[data$HypnoYN %in% c(0) & data$NewComp %in% c(1),]))
evalq(interaction.plot(NewH, lagC, accOr), e)

## Draw SE bars without participant means

grandMean <- mean(e$rtOr)
partMeans <-evalq(aggregate(list(partMean=rtOr), list(suj=suj), mean),
                  droplevels(data[data$hypnoC %in% c(4),]))
e <- merge(e, partMeans)
e$nrtOr <- e$rtOr-e$partMean+grandMean
f <- evalq(aggregate(list(sertOr=nrtOr), list(lagC=lagC,target1=target1), se), e)
f$sertOr <- f$sertOr*sqrt(6/3) # compute se and apply Morey's 2008 correction (For a 2X2 design M=4).
f$mrtOr <- evalq(aggregate(list(mrtOr=rtOr), list(lagC=lagC,target1=target1), mean), e)$mrtOr

evalq(interaction.plot(lagC, target1, rtOr) ,e)
j <- 0.001
xs <- c(1-j,2-j,3-j,4-j,1+j,2+j,3+j,4+j)
y1=f$mrtOr+f$sertOr
y0=f$mrtOr-f$sertOr

arrows(x0=xs, y0=y0, x1=xs, y1=y1, length=0,
       angle=90,code=3, lty=1, lwd=2)



##Compare no-Hypnosis with all hypnosis participants pooled together
#Accuracy
AccInc <- glmer(accOr~lagC*target1*target2*hypnoC*HypnoYN+(1|suj), family='binomial', 
                droplevels(data[data$hypnoC %in% c(1,2,3) & data$NewComp %in% c(1),]))




AccInc <- glmer(accOr~lagC*target1*target2*HLRAWcol*HypnoYN+(1|suj), family='binomial', 
                droplevels(data[data$hypnoC %in% c(1,2,3) & data$NewComp %in% c(1),]))

AccInc <- glmer(accOr~lagC*target1*target2*NewH+(1|suj), family='binomial', 
                droplevels(data[data$NewH %in% c(1,2,3) & data$NewComp %in% c(1) & 
                                data$HypnoYN==1 & data$nface==2,]))


AccInc2 <- glmer(accOr~lagC*target1*target2*NewH+(1|suj), family='binomial', 
                droplevels(data[data$NewH %in% c(1,2,3) & data$NewComp %in% c(1) & data$nface==2
                                & data$HypnoYN==0,]))

AccInc2 <- glmer(accOr~lagC*target1*target2*HypnoYN+(1|suj), family='binomial', 
                 droplevels(data[data$NewH %in% c(1,2,3) & data$NewComp %in% c(1) & data$nface==2
                                 ,]))

RtOr3 <- lmer(rtOr~lagC*target1*target2*HypnoYN*NewH+(1|suj), 
              droplevels(data[data$NewH %in% c(1,2,3) & data$NewComp %in% c(1) & data$nface==2
                              ,]))

RtOr2 <- lmer(rtOr~lagC*target1*target2*HypnoYN+(1|suj), 
                 droplevels(data[data$NewH %in% c(1,2,3) & data$NewComp %in% c(1) & data$nface==2
                                 ,]))
Anova(RtOr2)

Anova(AccInc2)

AccInc <- glmer(accOr~lagC*target2*hypnoC+(1|suj), family='binomial', 
                droplevels(data[data$hypnoC %in% c(1,2,3,4) & data$OldComp %in% c(1) & data$nface==1,]))

AccInc <- glmer(accOr~lagC*target1*target2+(1|suj), family='binomial', 
                droplevels(data[data$hypnoC %in% c(4) & data$OldComp %in% c(1),]))

AccInc2 <- glmer(accOr~lagC*hypnoC+(1|suj), family='binomial', 
                droplevels(data[data$hypnoC %in% c(1,2,3,4) & data$OldComp %in% c(1),]))

AccInc2 <- glmer(accOr~lagC*target2+(1|suj), family='binomial', 
                 droplevels(data[data$hypnoC %in% c(3,4) & data$OldComp %in% c(1) 
                                 & data$lagC %in% c("early"),]))




anova(AccInc,AccInc2)
Anova(AccInc2, type=2)
Anova(AccInc, type=2)

lsm.options(pbkrtest.limit = 9557)
lsmeans(AccInc, pairwise ~ target2*hypnoC*HypnoYN, adjust="tukey")
lsmeans(AccInc2, pairwise ~ NewH*target2, adjust="tukey")

AccInc2 <- glmer(accOr~hypnoC*target2+(1|suj), family='binomial', 
                 droplevels(data[data$hypnoC %in% c(1,2,3) & data$OldComp %in% c(1) 
                                 ,]))

AccInc3 <- glmer(accOr~lagC*hypnoC*target2+(1|suj), family='binomial', 
                 droplevels(data[data$hypnoC %in% c(1,2,3,4) & data$OldComp %in% c(1) 
                                 ,]))

#& data$lagC %in% c("late")

anova(AccInc2,AccInc3)
Anova(AccInc)
e <- evalq(aggregate(list(accOr=accOr), list(suj=suj, hypnoC=hypnoC, 
                    target1=target1, lagC=lagC, target2=target2), mean),
                    droplevels(data[data$hypnoC %in% c(1,2,3) & data$OldComp %in% c(1) ,]))

e <- evalq(aggregate(list(accOr=accOr), list(suj=suj, lagC=lagC, target1=target1, 
                                             target2=target2, NewH=NewH), mean),
           droplevels(data[data$NewH %in% c(1,2,3) & data$NewComp %in% c(1) 
                           & data$HypnoYN==0 & data$nface==2,]))

e <- evalq(aggregate(list(rtOr=rtOr), list(suj=suj, lagC=lagC, target1=target1, NewH=NewH, 
                                             target2=target2), mean),
           droplevels(data[data$NewH %in% c(1,2,3) & data$NewComp %in% c(1) 
                           & data$HypnoYN %in% c(0) & data$nface==2,]))


e <- evalq(aggregate(list(accOr=accOr), list(suj=suj, hypnoC=hypnoC, 
                                             target1=target1, lagC=lagC), mean),
           droplevels(data[data$hypnoC %in% c(3,4) & data$OldComp %in% c(1) ,]))


evalq(plot(target2, accOr), e)

evalq(interaction.plot(target1, target2, accOr), e[e$lagC %in% c("early"),])
evalq(interaction.plot(NewH,lagC,accOr), e)
evalq(interaction.plot(NewH, lagC, rtOr), e[e$HypnoYN==1 & e$target2=="neutral",])

evalq(interaction.plot(target1, lagC, rtOr), e)


evalq(interaction.plot(hypnoC, target2, accOr), e[e$HypnoYN==1,])
evalq(interaction.plot(lagC, hypnoC, accOr, ylim=c(0.6,0.9)), e[e$HypnoYN==0,])
evalq(interaction.plot(lagC, hypnoC, accOr, ylim=c(0.6,0.9)), e[e$HypnoYN==1,])
evalq(interaction.plot(target1, hypnoC, accOr, ylim=c(0.6,0.9)), e[e$HypnoYN==0,])
evalq(interaction.plot(target1, hypnoC, accOr, ylim=c(0.6,0.9)), e[e$HypnoYN==1,])
evalq(interaction.plot(target2, hypnoC, accOr, ylim=c(0.7,0.85)), e[e$HypnoYN==0,])
evalq(interaction.plot(target2, hypnoC, accOr, ylim=c(0.7,0.85)), e[e$HypnoYN==1,])

evalq(interaction.plot(target1, target2, accOr), e)
evalq(interaction.plot(target1, target2, accOr), e)
evalq(interaction.plot(lagC, hypnoC, accOr), e)


AccInc <- glmer(accOr~lagC*target1*target2*NewH+(1|suj), family='binomial', 
                droplevels(data[data$NewH %in% c(1,2,3) & data$NewComp %in% c(1) & 
                                  data$HypnoYN==1 & data$nface==2,]))


Anova(AccInc)

e <- evalq(aggregate(list(accOr=accOr), list(suj=suj, lagCatComplete=lagCatComplete, target1=target1, 
                                             target2=target2, NewH=NewH), mean),
           droplevels(data[data$NewH %in% c(1,2,3) & data$NewComp %in% c(1) 
                           & data$HypnoYN==0 & data$nface==2,]))



## Draw SE bars without participant means

grandMean <- mean(e$accOr)
partMeans <-evalq(aggregate(list(partMean=accOr), list(suj=suj), mean),
                  data)
e <- merge(e, partMeans)
e$naccOr <- e$accOr-e$partMean+grandMean
f <- evalq(aggregate(list(seaccOr=naccOr), list(NewH=NewH,lagC=lagC), se), e)
f$seaccOr <- f$seaccOr*sqrt(6/3) # compute se and apply Morey's 2008 correction (For a 2X2 design M=4).
f$maccOr <- evalq(aggregate(list(maccOr=accOr), list(NewH=NewH,lagC=lagC), mean), e)$maccOr

evalq(interaction.plot(NewH, lagC, accOr, ylim=c(0.65,0.9)) ,e)

#means<-aggregate(e[e$target2=="angry",], by=list(e$NewH[e$target2=="angry"]),mean)
#barplot(means$rtOr, ylim=c(700,1000))
#means<-aggregate(e[e$target2=="neutral",], by=list(e$NewH[e$target2=="neutral"]),mean)
#barplot(means$rtOr, ylim=c(700,1000))


j <- 0.001
xs <- c(1-j,2-j,3-j,1+j,2+j,3+j)
#xs <- c(1-j,2-j,1+j,2+j)
y1=f$maccOr+f$seaccOr
y0=f$maccOr-f$seaccOr

arrows(x0=xs, y0=y0, x1=xs, y1=y1, length=0,
       angle=90,code=3, lty=1, lwd=2)

#Response Time

RTInc <- lmer(rtOr~lagC*target1*target2+(1|suj), data=droplevels(data[data$hypnoC %in% c(4),]))
qqnorm(residuals(RTInc))
plot(fitted(RTInc),residuals(RTInc))
Anova(RTInc)
lsm.options(pbkrtest.limit = 9557)
lsmeans(RTInc, pairwise ~ target1*target2*lagC, adjust="tukey")

e <- evalq(aggregate(list(rtOr=rtOr), list(suj=suj, NewH=NewH, HypnoYN=HypnoYN,target1=target1, target2=target2, lagC=lagC), 
                     mean),droplevels(data[data$HypnoYN %in% c(0) & data$NewComp %in% c(1) & data$nface==2,]))
evalq(interaction.plot(NewH, lagC, rtOr), e)

## Draw SE bars without participant means

grandMean <- mean(e$rtOr)
partMeans <-evalq(aggregate(list(partMean=rtOr), list(suj=suj), mean),
                  droplevels(data[data$NewH %in% c(1,2,3),]))
e <- merge(e, partMeans)
e$nrtOr <- e$rtOr-e$partMean+grandMean
f <- evalq(aggregate(list(sertOr=nrtOr), list(NewH=NewH,lagC=lagC), se), e)
f$sertOr <- f$sertOr*sqrt(6/3) # compute se and apply Morey's 2008 correction (For a 2X2 design M=4).
f$mrtOr <- evalq(aggregate(list(mrtOr=rtOr), list(NewH=NewH,lagC=lagC), mean), e)$mrtOr

evalq(interaction.plot(NewH, lagC, rtOr, ylim=c(700,1030)) ,e)
j <- 0.001
xs <- c(1-j,2-j,3-j,1+j,2+j,3+j)
y1=f$mrtOr+f$sertOr
y0=f$mrtOr-f$sertOr

arrows(x0=xs, y0=y0, x1=xs, y1=y1, length=0,
       angle=90,code=3, lty=1, lwd=2)


#### ####
##                   ##
## HYPNOSIS ANALYSES ##
##                   ##

#Accuracy Inclination task
##########################

AccInc <- glmer(accOr~lagC*target1*target2*HLRAWcol*HypnoYN+(1|suj), family='binomial', data=data)
RtInc <- lmer(rtOr~lagC*target1*target2*HLRAWcol*HypnoYN+(1|suj), data=data)
Anova(RtInc)
lsmeans(RtInc, pairwise ~ lagC*target1*target2*HLRAWcol*HypnoYN, adjust="tukey")



AccInc <- glmer(accOr~lagC*target1*target2*hypnoC*HypnoYN+(1|suj), family='binomial', data=data)
AccInc <- glmer(accOr~lagC*target1*target2*newhypnoC*HypnoYN+(1|suj), family='binomial', data=data)
AccInc <- glmer(accOr~hypnoC+HypnoYN+(1|suj), family='binomial', data=data)
AccInc <- glmer(accOr~lagC*target1*target2*hypnoC*HypnoYN+(1|suj), family='binomial', data=data)
#qqnorm(residuals(AccInc))

### #### 
fix.formula <- accOr~lagC*target1*target2*hypnoC*HypnoYN+(1|suj)
x <- model.matrix(fix.formula, data)
fix.fit <- lm (fix.formula, data, method = "qr", singular.ok = TRUE)
oo <- lm (fix.formula, data, method = "qr", singular.ok = FALSE)
coef <- fix.fit$coef
p <- length(coef)
no.NA <- sum(is.na(coef))
rank <- fix.fit$rank #Test the rank deficiency... hypnoC is nested inside HypnoYN!!
### ####
summary(AccInc)
#plot(fitted(AccInc),residuals(AccInc))
Anova(AccInc)
lsm.options(pbkrtest.limit = 9557)
lsmeans(AccInc, pairwise ~ target2|hypnoC*target1|lagC, adjust="tukey")
lsmeans(AccInc, pairwise ~ hypnoC|target1*target2|lagC, adjust="tukey")
lsmeans(AccInc, pairwise ~ hypnoC*lagC|target1|target2, adjust="tukey")
lsmeans(AccInc, pairwise ~ hypnoC*target2, adjust="tukey")

e <- evalq(aggregate(list(accOr=accOr), list(suj=suj, hypnoC=hypnoC, HLRAWcol=HLRAWcol,
                                             target1=target1, target2=target2, hypno=hypno, 
                                             lagC=lagC, HypnoYN=HypnoYN),mean),data)
evalq(interaction.plot(lagC, hypnoC, accOr), e[e$HypnoYN=="yes",])
evalq(interaction.plot(hypnoC, target2, accOr), e[e$HypnoYN=="no",])
evalq(interaction.plot(hypno, HypnoYN, accOr), e[e$target2=="neutral",])
evalq(interaction.plot(hypnoC, target2, accOr), e[e$lagC=="late",])
evalq(interaction.plot(hypnoC, target2, accOr), e[e$target1=="neutral",])

## Draw SE bars without participant means

grandMean <- mean(e$accOr)
partMeans <-evalq(aggregate(list(partMean=accOr), list(suj=suj), mean),
                  data)
e <- merge(e, partMeans)
e$naccOr <- e$accOr-e$partMean+grandMean
f <- evalq(aggregate(list(seaccOr=naccOr), list(hypnoC=hypnoC,target2=target2), se), e)
f$seaccOr <- f$seaccOr*sqrt(6/3) # compute se and apply Morey's 2008 correction (For a 2X2 design M=4).
f$maccOr <- evalq(aggregate(list(maccOr=accOr), list(hypnoC=hypnoC,target2=target2), mean), e)$maccOr

evalq(interaction.plot(hypnoC, target2, accOr) ,e)


j <- 0.001
xs <- c(1-j,2-j,3-j,4-j,1+j,2+j,3+j,4+j)
y1=f$maccOr+f$seaccOr
y0=f$maccOr-f$seaccOr

arrows(x0=xs, y0=y0, x1=xs, y1=y1, length=0,
       angle=90,code=3, lty=1, lwd=2)

#Accuracy Face Number task
##########################

AccNF <- glmer(accn~lagC*target1*target2*hypnoC+(1|suj), family='binomial', data=data)
#qqnorm(residuals(AccInc))
#plot(fitted(AccInc),residuals(AccInc))
Anova(AccNF)
lsm.options(pbkrtest.limit = 9557)
lsmeans(AccNF, pairwise ~ hypnoC | hypnoC*target2, adjust="tukey")

e <- evalq(aggregate(list(accn=accn), list(suj=suj, hypnoC=hypnoC, target1=target1, target2=target2, lagC=lagC), 
                     mean),data)
evalq(interaction.plot(target2, hypnoC, accn), e)
evalq(interaction.plot(hypnoC, target2, accn), e[e$target1=="angry",])
evalq(interaction.plot(hypnoC, target2, accn), e[e$target1=="neutral",])

## Draw SE bars without participant means

grandMean <- mean(e$accn)
partMeans <-evalq(aggregate(list(partMean=accn), list(suj=suj), mean),
                  data)
e <- merge(e, partMeans)
e$naccn <- e$accn-e$partMean+grandMean
f <- evalq(aggregate(list(seaccn=naccn), list(hypnoC=hypnoC,target2=target2), se), e)
f$seaccn <- f$seaccn*sqrt(6/3) # compute se and apply Morey's 2008 correction (For a 2X2 design M=4).
f$maccn <- evalq(aggregate(list(maccn=accn), list(hypnoC=hypnoC,target2=target2), mean), e)$maccn

evalq(interaction.plot(hypnoC, target2, accn) ,e)
j <- 0.001
xs <- c(1-j,2-j,3-j,4-j,1+j,2+j,3+j,4+j)
y1=f$maccn+f$seaccn
y0=f$maccn-f$seaccn

arrows(x0=xs, y0=y0, x1=xs, y1=y1, length=0,
       angle=90,code=3, lty=1, lwd=2)



#Response Time Inclination task
###############################

RTInc <- lmer(rtOr~lagC*target1*target2*hypnoC+(1|suj), data=data)
qqnorm(residuals(RTInc))
plot(fitted(RTInc),residuals(RTInc))
Anova(RTInc)
lsm.options(pbkrtest.limit = 9557)
lsmeans(RTInc, pairwise ~ hypnoC | target1*target2*lagC, adjust="tukey")

e <- evalq(aggregate(list(rtOr=rtOr), list(suj=suj, hypnoC=hypnoC, target1=target1, target2=target2, lagC=lagC), 
                     mean),data)
evalq(interaction.plot(hypnoC, target2, rtOr), e)

## Draw SE bars without participant means

grandMean <- mean(e$rtOr)
partMeans <-evalq(aggregate(list(partMean=rtOr), list(suj=suj), mean),
                  data)
e <- merge(e, partMeans)
e$nrtOr <- e$rtOr-e$partMean+grandMean
f <- evalq(aggregate(list(sertOr=nrtOr), list(hypnoC=hypnoC,target2=target2), se), e)
f$sertOr <- f$sertOr*sqrt(6/3) # compute se and apply Morey's 2008 correction (For a 2X2 design M=4).
f$mrtOr <- evalq(aggregate(list(mrtOr=rtOr), list(hypnoC=hypnoC,target2=target2), mean), e)$mrtOr

evalq(interaction.plot(hypnoC, target2, rtOr) ,e)
means<-aggregate(e[e$target2=="angry",], by=list(e$hypnoC[e$target2=="angry"]),mean)
barplot(means$rtOr, ylim=c(700,1000))
means<-aggregate(e[e$target2=="neutral",], by=list(e$hypnoC[e$target2=="neutral"]),mean)
barplot(means$rtOr, ylim=c(700,1000))

j <- 0.001
xs <- c(1-j,2-j,3-j,4-j,1+j,2+j,3+j,4+j)
y1=f$mrtOr+f$sertOr
y0=f$mrtOr-f$sertOr

arrows(x0=xs, y0=y0, x1=xs, y1=y1, length=0,
       angle=90,code=3, lty=1, lwd=2)


#Response Time Face Number task
###############################

RTNF <- lmer(rtn~lagC*target1*target2*hypnoC+(1|suj), data=data)
qqnorm(residuals(RTNF))
plot(fitted(RTInc),residuals(RTInc))
Anova(RTNF)
lsm.options(pbkrtest.limit = 9557)
lsmeans(RTInc, pairwise ~ hypnoC | target1*target2*lagC, adjust="tukey")

e <- evalq(aggregate(list(rtn=rtn), list(suj=suj, hypnoC=hypnoC, target1=target1, target2=target2, lagC=lagC), 
                     mean),data)
evalq(interaction.plot(target2, hypnoC, rtn), e)

## Draw SE bars without participant means

grandMean <- mean(e$rtn)
partMeans <-evalq(aggregate(list(partMean=rtn), list(suj=suj), mean),
                  data)
e <- merge(e, partMeans)
e$nrtn <- e$rtn-e$partMean+grandMean
f <- evalq(aggregate(list(sertn=nrtn), list(hypnoC=hypnoC,target2=target2), se), e)
f$sertn <- f$sertn*sqrt(6/3) # compute se and apply Morey's 2008 correction (For a 2X2 design M=4).
f$mrtn <- evalq(aggregate(list(mrtn=rtn), list(hypnoC=hypnoC,target2=target2), mean), e)$mrtn

evalq(interaction.plot(hypnoC, target2, rtn) ,e)
j <- 0.001
xs <- c(1-j,2-j,3-j,4-j,1+j,2+j,3+j,4+j)
y1=f$mrtn+f$sertn
y0=f$mrtn-f$sertn

arrows(x0=xs, y0=y0, x1=xs, y1=y1, length=0,
       angle=90,code=3, lty=1, lwd=2)



#Statistics for the diff between T2Angry - T2Neutral
### ################################################

e <- evalq(aggregate(list(accOr=accOr), list(suj=suj, hypnoC=hypnoC, target1=target1, lagC=lagC, target2=target2), 
                     mean),droplevels(data[data$OldComp==1,]))
evalq(interaction.plot(lagC, target2, accOr), e[e$hypnoC %in% c(1,2,3,4),])
evalq(interaction.plot(lagC, hypnoC, accOr), e)

e$accdif<-e[e$target2=="angry",]$accOr-e[e$target2=="neutral",]$accOr

l <- aov(accdif~hypnoC+Error(suj/(target1*target2)), data=e)
l <- aov(accdif~target1*target2*hypnoC*lagC, data=e)
#l <- aov(accdif~target1*target2*hypnoC+Error(suj/(target1*target2)), data=e[e$lagC=="late",])
#l <- aov(accdif~target1*target2*hypnoC, data=e[e$lagC=="late",])
summary(l)

mod<-lm(accdif ~ hypnoC, data=e[e$hypnoC %in% c(1,2,3,4),])
Anova(mod)
lsmeans(mod, pairwise ~ hypnoC, adjust="tukey")

#l <- lme(accdif ~ lagC*target1*target2*hypnoC, data=e, random = ~1|suj)
#anova(l)
#require(multcomp)
#summary(glht(l, linfct=mcp(hypnoC="Tukey")))
#lsmeans(l, pairwise ~ hypnoC | target1*target2*lagC, adjust="tukey")

evalq(interaction.plot(hypnoC, lagC,accdif),e)
evalq(interaction.plot(hypnoC, lagC,accdif),e[e$target1=="angry",])
evalq(interaction.plot(hypnoC, lagC,accdif),e[e$target1=="neutral",])
#TukeyHSD(x=l, 'hypnoC', conf.level=0.95)

## ####
#DIFFUSION MODEL
################

## First we plot densities
##### #####

#By Category of H
densLCorr <- density(data$rtOr[data$accOr==1 & data$hypnoC==1])
densMCorr <- density(data$rtOr[data$accOr==1 & data$hypnoC==2])
densHCorr <- density(data$rtOr[data$accOr==1 & data$hypnoC==3])
#densNHCorr <- density(data$rtOr[data$accOr==1 & data$hypnoC==4])
densLInc <- density(data$rtOr[data$accOr==0 & data$hypnoC==1])
densMInc <- density(data$rtOr[data$accOr==0 & data$hypnoC==2])
densHInc <- density(data$rtOr[data$accOr==0 & data$hypnoC==3])
#densNHInc <- density(data$rtOr[data$accOr==0 & data$hypnoC==4])

#Hypno
densHCorr <- density(data$rtOr[data$accOr==1 & data$HypnoYN=="yes"])
densHInc <- density(data$rtOr[data$accOr==0 & data$HypnoYN=="yes"])
densNHCorr <- density(data$rtOr[data$accOr==1 & data$HypnoYN=="no"])
densNHInc <- density(data$rtOr[data$accOr==0 & data$HypnoYN=="no"])


M <- max(c(densLCorr$y,densMCorr$y,densHCorr$y,densNHCorr$y, densLInc$y,densMInc$y,densHInc$y,
           densNHInc$y))
#pdf("Densities.pdf", width = 12, height = 6)
plot(densLCorr, col=2, lwd=3, ylim=c(0, M), xlim=c(0, 3500))
lines(densMCorr, col=1, lwd=3); lines(densHCorr, col=3, lwd=3); lines(densNHCorr, col=4, lwd=3) 
lines(densLInc, col=2, lwd=1, lty=2); lines(densMInc, col=1, lwd=1, lty=2); lines(densHInc, col=3, lwd=1, lty=2)
lines(densNHInc, col=4, lwd=1, lty=2)
#dev.off()

M <- max(c(densHCorr$y,densNHCorr$y,densHInc$y,densNHInc$y))
#pdf("Densities.pdf", width = 12, height = 6)
plot(densHCorr, col=2, lwd=3, ylim=c(0, M), xlim=c(0, 3500))
lines(densNHCorr, col=1, lwd=3); lines(densHInc, col=2, lwd=3, lty=2); lines(densNHInc, col=1, lwd=3, lty=2)
#dev.off()


#### ####
#For running fast-dm (if it's the first time, otherwise load Fast-dm_diffusion_parameters.R)
###################### #####################################################################
sujs <- unique(data$suj[data$NewComp==1 & data$HypnoYN==0]);
i <- 1
#df <- data ## All data
df <- droplevels(data[data$NewH %in% c(1,2,3) & data$NewComp %in% c(1) 
                      & data$HypnoYN==0,]) ## All hypno ##  & data$nface==2 in case I wanna kill absent t1s
#df <- droplevels(data[data$NewH %in% c(1,2,3) & data$NewComp %in% c(1) 
#                      & data$HypnoYN==0,]) ## All no hypno ##  & data$nface==2 in case I wanna kill absent t1s
## ## We recode t1 et t2 levels, because it makes things easier
## ## aftewards.
df$T1 <- paste(df$target1, "t1", sep='')
df[df$nface==1,]$T1 <- 'absent'
df$T2 <- paste(df$target2, "t2", sep='')

#makepath = function(path){dir.create(dirname(path),recursive=TRUE);path} #because R in windows will not go make the file otherwise
## ## We write the data for each subject in different files, in different
## ## folder so as to parallelize by hand. DONT FORGET TO CHANGE THE F**CKING ROOT DIR
 for (s in sujs){
     #folder <- paste("C:/Users/subcomarc/Google Drive/R/DM/d", ceiling(i/12), "/", sep='') #put no more than 13 per folder
     temp <- evalq(data.frame(acc=accOr, rt=rtOr/1000, target1=T1,
                              target2=T2, lag=lagC),
                   df[df$suj==s,])
     #file=as.character(paste(folder, s, ".txt", sep=''))
     file=as.character(paste(s, ".txt", sep=''))
     #write.table(temp, makepath(file), row.names=FALSE, col.names=FALSE, quote=FALSE)
     write.table(temp, file, row.names=FALSE, col.names=FALSE, quote=FALSE)
     i <- i+1
 }

##RUN FASTDM AT THIS POINT

#paste("C:/Users/subcomarc/Google Drive/R/DM/d", n,

 #fdm <- do.call(rbind, lapply(1:5, function(n){ #fdmImp
              fdm <- readFastDMComplete(file='blinkFinal.txt',
                                           dvs=c('v','a','t0'),
                                           factors=c('target1', 'target2', 'lag'),
                                           flevels=list(c('absent', 'angryt1', 'neutralt1'),
                                           #flevels=list(c('angryt1', 'neutralt1'), #for dataframes without absent T1 trials
                                               c('angryt2', 'neutralt2'),
                                               c('early', 'late')))
  #        }))


 f <- evalq(aggregate(list(HLRAWcol=HLRAWcol, NewH=NewH, HypnoYN=HypnoYN), 
                      list(suj=suj), unique), data)
fdm <- merge(fdm, f)
 #save(fdm, file="Fast-dm_diffusion_parameters_hypno.R")
 #save(fdm, file="Fast-dm_diffusion_parameters_Nohypno.R")
 save(fdm, file="Fast-dm_diffusion_parameters_ALL.R")
###################### ###################### 

load("Fast-dm_diffusion_parameters_ALL.R")
load("Fast-dm_diffusion_parameters_hypno.R")
load("Fast-dm_diffusion_parameters_Nohypno.R")
 
 
fdm$target1<-factor(fdm$target1)  
fdm$target2<-factor(fdm$target2)  
fdm$lag<-factor(fdm$lag)
fdm$NewH<-factor(fdm$NewH)
fdm$HypnoYN<-factor(fdm$HypnoYN)


contrasts(fdm$target1)=contr.sum(3)
contrasts(fdm$target2)=contr.sum(2)
contrasts(fdm$lag)=contr.sum(2)
contrasts(fdm$NewH)=contr.sum(3)
contrasts(fdm$HypnoYN)=contr.sum(2)

#Test plot like scatter

fdm$NewH<-as.numeric(as.character(fdm$NewH))
fdm$NewH<-as.factor(fdm$NewH)
e <- evalq(aggregate(list(a=a), list(suj=suj, NewH=NewH), mean),droplevels(fdm[fdm$target1 %in% c("angryt1","neutralt1")
                                          & fdm$HypnoYN %in% c(1),]))
fit1<-lm(a~NewH,data=e)
ggplot(e, aes(x=NewH, y=a, color=NewH)) +
  #scale_shape_discrete(solid=T, legend=F)
  geom_point(shape=1)+
  stat_smooth(method='lm',  col="purple")

#scale_shape_discrete(solid=T, legend=F)
geom_point(data=e, mapping=aes(x=NewH, y=a, color=NewH, size="0.5"))

#Separate H from no H with stats over everyone
load("Fast-dm_diffusion_parameters_hypno.R")
load("Fast-dm_diffusion_parameters_Nohypno.R")
load("Fast-dm_diffusion_parameters_ALL.R")
DMa <- lmer(a~lag*target1*target2*NewH*HypnoYN+(1|suj), data=droplevels(fdm[fdm$target1 %in% c("angryt1","neutralt1"),]))
DMv <- lmer(v~lag*target1*target2*NewH*HypnoYN+(1|suj), data=droplevels(fdm[fdm$target1 %in% c("angryt1","neutralt1"),]))
DMt0 <- lmer(t0~lag*target1*target2*NewH*HypnoYN+(1|suj), data=droplevels(fdm[fdm$target1 %in% c("angryt1","neutralt1"),]))

DMa <- lmer(a~lag*target1*target2*NewH+(1|suj), data=droplevels(fdm[fdm$target1 %in% c("angryt1","neutralt1"),]))
DMv <- lmer(v~lag*target1*target2*NewH+(1|suj), data=droplevels(fdm[fdm$target1 %in% c("angryt1","neutralt1"),]))
DMt0 <- lmer(t0~lag*target1*target2*NewH+(1|suj), data=droplevels(fdm[fdm$target1 %in% c("angryt1","neutralt1"),]))


Anova(DMa)
Anova(DMv)
Anova(DMt0)

#For a (threshold separation, liberal vs conservative)

load("Fast-dm_diffusion_parameters_ALL.R")
load("Fast-dm_diffusion_parameters_hypno.R")
load("Fast-dm_diffusion_parameters_Nohypno.R")


DMa <- lmer(a~lag*target1*target2*NewH+(1|suj), data=fdm)


 DMa <- lmer(a~lag*target1*target2*NewH+(1|suj), data=droplevels(fdm[fdm$target1 %in% c("angryt1","neutralt1") 
                                                                             & fdm$HypnoYN %in% c(0),]))
 DMa <- lmer(a~lag*target1*target2+(1|suj), data=droplevels(fdm[fdm$target1 %in% c("angryt1","neutralt1"),]))
 DMa <- lmer(a~lag*target1*target2*NewH+(1|suj), fdm)

 qqnorm(residuals(DMa))
 plot(fitted(DMa),residuals(DMa))
 Anova(DMa)
 lsm.options(pbkrtest.limit = 9557)
 lsmeans(DMa, pairwise ~ HypnoYN | lag*target2 |target1, adjust="tukey")
 
 e <- evalq(aggregate(list(a=a), list(suj=suj, 
                                      NewH=NewH,
                                      HypnoYN=HypnoYN), 
                      mean),droplevels(fdm[fdm$target1 %in% c("angryt1","neutralt1"),]))
 e <- evalq(aggregate(list(a=a), list(suj=suj, 
                                      NewH=NewH), 
                      mean),droplevels(fdm[fdm$target1 %in% c("angryt1","neutralt1")
                                           & fdm$HypnoYN %in% c(0),]))
 e <- evalq(aggregate(list(a=a), list(suj=suj, hypnoC=hypnoC,target1=target1, target2=target2, lag=lag), 
                      mean),fdm)
 evalq(interaction.plot(hypno, HypnoYN, a, ylim=c(-2,9)), e[e$lag=="late",])
 evalq(interaction.plot(lag, target1, a), e[e$target1 %in% c("angryt1","neutralt1"),])
 evalq(interaction.plot(NewH, target2, a), e)
 
 Pt<- evalq(aggregate(list(a=a), list(NewH=NewH), mean), e)
 evalq(qplot(NewH, a, ylim=c(1,2.5), geom="jitter"), e)
 evalq(plot(NewH, a, border="white", ylim=c(1.40,1.80)), e)
 evalq(lines(Pt))
 

 
 require(ggplot2)
 qplot(NewH, a , data=e)
 lines(e$NewH, e$a)
 
 ## Draw SE bars without participant means
 
 grandMean <- mean(e$a)
 partMeans <-evalq(aggregate(list(partMean=a), list(suj=suj), mean),
                   fdm)
 e <- merge(e, partMeans)
 e$na <- e$a-e$partMean+grandMean
 f <- evalq(aggregate(list(sea=na), list(NewH=NewH), se), e)
 f$sea <- f$sea*sqrt(6/3) # compute se and apply Morey's 2008 correction (For a 2X2 design M=4).
 f$ma <- evalq(aggregate(list(ma=a), list(NewH=NewH), mean), e)$ma
 

 #evalq(lines(Pt))
 #evalq(interaction.plot(hypnoC, target2, a) ,e)
 
 evalq(plot(ma, type="l", ylim=c(1.40,1.80)), f) #ylim=c(1.40,1.80)),
 evalq(interaction.plot(NewH, HypnoYN, a), e) #ylim=c(1.40,1.80)),
 
 
 j <- 0.001
 xs <- c(1-j,2-j,3-j,1+j,2+j,3+j)
 y1=f$ma+f$sea
 y0=f$ma-f$sea
 
 arrows(x0=xs, y0=y0, x1=xs, y1=y1, length=0,
        angle=90,code=3, lty=1, lwd=2)
 
  
#For v (drift rate. Average slope of the information accumulation process)

load("Fast-dm_diffusion_parameters_ALL.R")
load("Fast-dm_diffusion_parameters_hypno.R")
load("Fast-dm_diffusion_parameters_Nohypno.R")

 
 DMv <- lmer(v~lag*target1*target2*hypno*HypnoYN+(1|suj), data=fdm)
 
 DMv <- lmer(v~lag*target1*target2*NewH+(1|suj), data=droplevels(fdm[fdm$target1 %in% c("angryt1","neutralt1"),]))
 DMv <- lmer(v~lag*target1*target2+(1|suj), data=droplevels(fdm[fdm$target1 %in% c("angryt1","neutralt1"),]))
 qqnorm(residuals(DMv))
 plot(fitted(DMv),residuals(DMv))
 Anova(DMv)
 lsm.options(pbkrtest.limit = 9557)
 lsmeans(DMv, pairwise ~ lag*target1, adjust="tukey")
 
 e <- evalq(aggregate(list(v=v), list(suj=suj, target1=target1, lag=lag), 
                      mean),droplevels(fdm[fdm$target1 %in% c("angryt1","neutralt1"),]))
 evalq(interaction.plot(lag, target2, v), e[e$target1=="absent",])
 evalq(interaction.plot(lag, target1, v), e)
 
 ## Draw SE bars without participant means
 
 grandMean <- mean(e$v)
 partMeans <-evalq(aggregate(list(partMean=v), list(suj=suj), mean),
                   fdm)
 e <- merge(e, partMeans)
 e$nv <- e$v-e$partMean+grandMean
 f <- evalq(aggregate(list(sev=nv), list(lag=lag,target1=target1), se), e)
 f$sev <- f$sea*sqrt(6/3) # compute se and apply Morey's 2008 correction (For a 2X2 design M=4).
 f$mv <- evalq(aggregate(list(mv=v), list(lag=lag,target1=target1), mean), e)$mv
 
 evalq(interaction.plot(lag, target1, ylim=c(0.29,1.1), v) ,e)
 j <- 0.001
 xs <- c(1-j,2-j,1+j,2+j)
 y1=f$mv+f$sev
 y0=f$mv-f$sev
 
 arrows(x0=xs, y0=y0, x1=xs, y1=y1, length=0,
        angle=90,code=3, lty=1, lwd=2)
 
#For t0 (non-decision time or response time constant (in seconds). 
#Lower bound for the duration of all non-decisional processes (encoding and response execution))

load("Fast-dm_diffusion_parameters_ALL.R")
load("Fast-dm_diffusion_parameters_hypno.R")
load("Fast-dm_diffusion_parameters_Nohypno.R")
 
 DMt0 <- lmer(t0~lag*target1*target2*hypno*HypnoYN+(1|suj), data=fdm)
 
 DMt0 <- lmer(t0~lag*target1*target2*NewH+(1|suj), data=droplevels(fdm[fdm$target1 %in% c("angryt1","neutralt1"),]))
 DMt0 <- lmer(t0~lag*target1*target2+(1|suj), data=droplevels(fdm[fdm$target1 %in% c("angryt1","neutralt1"),]))
 qqnorm(residuals(DMt0))
 plot(fitted(DMt0),residuals(DMt0))
 Anova(DMt0)
 lsm.options(pbkrtest.limit = 9557)
 lsmeans(DMt0, pairwise ~ NewH*lag, adjust="tukey")
 
 e <- evalq(aggregate(list(t0=t0), list(suj=suj, NewH=NewH, lag=lag), 
                      mean),droplevels(fdm[fdm$target1 %in% c("angryt1","neutralt1"),]))
 evalq(interaction.plot(lag, NewH, t0), e)
 
 ## Draw SE bars without participant means
 
 grandMean <- mean(e$t0)
 partMeans <-evalq(aggregate(list(partMean=t0), list(suj=suj), mean),
                   fdm)
 e <- merge(e, partMeans)
 e$nt0 <- e$t0-e$partMean+grandMean
 f <- evalq(aggregate(list(set0=nt0), list(NewH=NewH,lag=lag), se), e)
 f$set0 <- f$set0*sqrt(6/3) # compute se and apply Morey's 2008 correction (For a 2X2 design M=4).
 f$mt0 <- evalq(aggregate(list(mt0=t0), list(NewH=NewH,lag=lag), mean), e)$mt0
 
 evalq(interaction.plot(NewH, lag, t0, ylim=c(0.21,0.29)) ,e)
 j <- 0.001
 xs <- c(1-j,2-j,3-j,1+j,2+j,3+j)
 y1=f$mt0+f$set0
 y0=f$mt0-f$set0
 
 arrows(x0=xs, y0=y0, x1=xs, y1=y1, length=0,
        angle=90,code=3, lty=1, lwd=2)
 
