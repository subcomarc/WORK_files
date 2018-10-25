##Analysis Script for the HypnoTunnel Experiment
##By subcomarc, with the resources of J?r?me Sackur
## Needs file "TunnelData.r"

################################################
#Preloading settings, directories and libraries#
################################################

#setwd("C:/Users/lscpuser/Google Drive/R/")
#data <- read.table(file="C:/Users/lscpuser/Google Drive/R/DataTable RAW.txt", header=T)
#source("C:/Users/lscpuser/Google Drive/R/R_Lib_Sackur.r")
#setwd("C:/Users/subcomarc/Google Drive/R/")
setwd("/home/subcomarc/Google_Drive/R/")
#data <- read.table(file="C:/Users/subcomarc/Google Drive/R/DataTable RAW.txt", header=T)
#source("C:/Users/subcomarc/Google Drive/R/R_Lib_Sackur.r")
source("/home/subcomarc/Google_Drive/R/R_Lib_Sackur.r")

library('xtable')
library('knitr')
library('lme4')
library('car')
library('multcomp')
library('lsmeans')
library("Hmisc")
library("HH")
library("xlsx")

options(xtable.floating = FALSE)
options(xtable.timestamp = "")

##################################################################
#Load preprocessed data (for preprocessing see analysesTUNNEL2.r)#
##################################################################

load("TunnelData.r")

###################### ###################### ######################
#STATISTICAL ANALYSES#
###################### ###################### ######################

#Prepping the data for stats and adding some convenient columns

data<-data
#data$HarSc[data$suj==1121]<-10 #Correct one wrong hypnotizability raw score
data$NewRaw<-data$HarSc
#data<-droplevels(data[!data$suj %in% c(1091,1101,1121,2111),])#Eliminate subs with dubious hypno reactions during test (beware of ,2111)
#data<-droplevels(data[!data$suj %in% c(1091,1101,2111),])#Eliminate subs with dubious hypno reactions during test (beware of ,2111)
#data$HarSc[data$suj==1091]<-9
#data$HarSc[data$suj==1101]<-9
data$rep <- 1 #for repetition priming
data$rep[data$Target!=data$Prime] <- 0
data$perf <- data$rt #for having rt for all blocks on the same column
data$perf[data$block==2] <- data$rtVis[data$block==2]
data$RTimeVis<-round(data$RTimeVis*1000)
data$RTimeVis[data$block %in% c(1,3)]<-5000
data <- droplevels(data[data$RTimeVis<5001,])
data<- droplevels(data[!(data$block!=2 & data$rt>3500),])
#data<- droplevels(data[!(data$block!=2 & data$rt>1500),])
#data<- droplevels(data[!(data$block!=2 & data$rt<200),])
data$IsSeen<-1 #for treating visibility as a binary
data$IsSeen[data$Visibility<2]<-0 #for treating visibility as a binary
data$suj<-as.numeric(as.character(data$suj))
data$NewH<-2 #Medium
data$NewH[data$HarSc>7]<-3 #High
data$NewH[data$HarSc<4]<-1 #Low
data$Induction<-1
data$Induction[data$hyp==3]<-0
data$NewH[data$suj %in% c(1091,1101,1121)]<-3 #Actually highs
data$NewRaw[data$suj %in% c(1091,1101,1121)]<-9 #Actually highs
data$NewH[data$suj %in% c(2111)]<-1 #Actually Lows
data$NewRaw[data$suj %in% c(2111)]<-3 #Actually lows
data$NewRaw[data$suj %in% c(3061)]<-9 #Correct Raw hypno fault
data$NewRaw[data$suj %in% c(3181)]<-10 #
data$Subli<-0
data$Subli[data$dur %in% c(17)]<-1
data$AllH <- data$NewH
data$AllH[data$hyp %in% c(3) & data$NewH %in% c(3)] <- 4
data$AllH[data$Induction %in% c(1) & data$NewH %in% c(1)] <- 5

#Add Age and Gender from original dataset#

Subs <- unique(data$SubjName)
Gender <- c(1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1)#(Female=1, Male=0)
Age <- c(34, 24, 21, 24, 21, 24, 23, 28, 24, 22, 22, 26, 28, 24, 28, 32, 24, 24, 22, 19, 25, 24, 27, 25, 19, 22, 26, 26, 23, 22, 21, 24, 29, 35, 21, 27, 22, 24, 21, 23, 26, 32, 23, 20)
e <- data.frame(cbind(Subs, Gender, Age))
data$Gender <- 2
data$Age <- 2
for(i in Subs){
data[data$SubjName == i,]$Gender <- e[e$Subs == i,]$Gender
data[data$SubjName == i,]$Age <- e[e$Subs == i,]$Age
}
e$Hypnotiz <- 2
for(i in Subs){
  e[e$Subs == i,]$Hypnotiz <- data[data$SubjName == i,]$NewRaw[1]
}

#data<-droplevels(data[!(data$hyp==3 & data$HarSc<4),])#for a control group with only mediums
#data<-droplevels(data[!(data$hyp==3 & data$HarSc>8),])#for a control group with only mediums

#Save a page of the data for the online open repository: dataiosf

dataiosf <- evalq(data.frame(SubjName, BlockNum, TrialN, rt, rtVis, acc, 
                             IsPrimed, Visibility, dur, Target, Prime, Quadrant, 
                             PrimeNoPrime, Luminance, centralTask, periphTask, periphVis, 
                             NewRaw, NewH, Induction, Subli, Gender, Age, stringsAsFactors = FALSE),
                             data[(data$NewH %in% c(3)) | (data$NewH %in% c(1) & data$Induction %in% c(1)),])

dataiosf[dataiosf$BlockNum %in% c(1),]$BlockNum <- "OC"
dataiosf[dataiosf$BlockNum %in% c(2),]$BlockNum <- "SP"
dataiosf[dataiosf$BlockNum %in% c(3),]$BlockNum <- "OP"
dataiosf[dataiosf$BlockNum %in% c(4),]$BlockNum <- "OCSP"

names(dataiosf) <-  c("Participant","BlockName","TrialNumber","ResponseTime","ResponseTimeVisibilityTask", 
"Accuracy","Congruency","VisibilityofPeriphTarget","DurationofPeriphTarget","CentralTarget","PeriphTarget",
"Quadrant","PresenceofPeriphTargt", "Luminance", "CentralTask","PeriphTask", "PeriphVisibility",
"RawHarvadScore", "HypnotizabilityScore", "HypnoticInduction","SubliminalTarget","Gender",
"Age")

save(dataiosf, file="/home/subcomarc/Google_Drive/R/DataHT.R")
write.xlsx(dataiosf, "/home/subcomarc/Google_Drive/R/DataHT.xlsx")

#Prepare factors and data without baseline control

data2 <- data[data$dur>0 & data$PrimeNoPrime==1 & data$IsPrimed!=2,] #to have a dataset without the trials without prime


data$block <- factor(data$block)
data$suj <- factor(data$suj)
data$hyp <- factor(data$hyp)
data$NewH<-factor(data$NewH)
data$dur <- factor(data$dur)
data$Quadrant <- factor(data$Quadrant)
data$IsPrimed <- factor(data$IsPrimed)
data$rep<- factor(data$rep)
data$Induction <- factor(data$Induction)
data$Subli<- factor(data$Subli)
data$AllH <- factor(data$AllH)

data2$block <- factor(data2$block)
data2$hyp <- factor(data2$hyp)
data2$NewH<-factor(data2$NewH, ordered=TRUE)
data2$dur <- factor(data2$dur)
data2$Quadrant <- factor(data2$Quadrant)
data2$IsPrimed <- factor(data2$IsPrimed)
data2$rep<- factor(data2$rep)
data2$HypnoYN<-factor(data2$HypnoYN)
data2$Induction <- factor(data2$Induction)
data2$Subli<- factor(data2$Subli)
data2$AllH <- factor(data2$AllH)

data2$NewH<-ordered(data2$NewH)
#data2$NewRaw<-ordered(data2$NewRaw)



contrasts(data$Quadrant)=contr.sum(4)
contrasts(data$block)=contr.sum(4)
contrasts(data$hyp)=contr.sum(3)
contrasts(data$NewH)=contr.sum(3)
contrasts(data$dur)=contr.sum(5)
contrasts(data$rep)=contr.sum(2)
contrasts(data$IsPrimed)=contr.sum(3)
contrasts(data$Induction)=contr.sum(2)
contrasts(data$Subli)=contr.sum(2)
contrasts(data$AllH)=contr.sum(5)

contrasts(data2$Quadrant)=contr.sum(4)
contrasts(data2$block)=contr.sum(4)
contrasts(data2$hyp)=contr.sum(3)
contrasts(data2$NewH)=contr.sum(3)
contrasts(data2$dur)=contr.sum(4)
contrasts(data2$rep)=contr.sum(2)
contrasts(data2$IsPrimed)=contr.sum(2)
contrasts(data2$HypnoYN)=contr.sum(2)
contrasts(data2$Induction)=contr.sum(2)
contrasts(data2$Subli)=contr.sum(2)
contrasts(data2$AllH)=contr.sum(5)

#contrasts(data2$NewRaw)=contr.sum(13)

control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))

#control=glmerControl(optimizer="Nelder_Mead", optCtrl = list(maxfun = 100000))


######################################################
#SD and mean for all blocks for all hyp for ACC and RT
######################################################

Block1 <- evalq(aggregate(list(rt=rt, acc=acc), list(suj=suj,IsPrimed=IsPrimed, dur=dur, hyp=hyp)
                          , mean), data2[data2$block==1,])

Block2 <- evalq(aggregate(list(RTimeVis=RTimeVis, Visibility=Visibility), 
                          list(suj=suj, IsPrimed=IsPrimed, dur=dur, hyp=hyp), 
                          mean), data2[data2$block==2,])

Block3 <- evalq(aggregate(list(rt=rt, acc=acc), list(suj=suj,IsPrimed=IsPrimed, dur=dur, hyp=hyp), 
                          mean),data2[data2$block==3,])

Block4 <- evalq(aggregate(list(rt=rt, RTimeVis=RTimeVis, acc=acc, Visibility=Visibility), 
                          list(suj=suj, IsPrimed=IsPrimed, dur=dur, hyp=hyp),mean),
                data2[data2$block==4,])

Mean1=mean(Block1$acc);SD1=sd(Block1$acc);Rat1=Mean1/SD1
Mean2=mean(Block2$Visibility);SD2=sd(Block2$Visibility);Rat2=Mean2/SD2
Mean3=mean(Block3$acc);SD3=sd(Block3$acc);Rat3=Mean3/SD3
Mean4=mean(Block4$acc);SD4=sd(Block4$acc);Rat4=Mean4/SD4


#############
#Performances
#############

#Accuracy


#for highs and lows i+s

CentrTAcc<-glmer(acc~NewH*Induction*dur*IsPrimed+(1|suj), family="binomial", droplevels(data2[data2$block %in% c(1,4) 
                                                                    & data2$NewH %in% c(1,3),]))



##LATEST VERSION HERE VV

#With Raw

CentrTAcc<-glmer(acc~NewRaw*Induction*IsPrimed*dur+(1|suj), family="binomial", droplevels(data2[data2$block %in% c(1,4) 
                                                                                              & data2$NewRaw %in% c(0:4,8:12),]))

CentrTAcc2<-glmer(acc~Induction*IsPrimed+(1|suj), family="binomial", droplevels(data2[data2$block %in% c(1,4) 
                                                                                              & data2$NewRaw %in% c(8:12),]))
                                                                                             #& data2$NewRaw %in% c(0:12),]))

CentrTAcc3<-glmer(acc~NewRaw*IsPrimed+(1|suj), family="binomial", droplevels(data2[data2$block %in% c(1,4)
                                                                                   & data2$Induction %in% c(1)
                                                                                   & data2$NewRaw %in% c(0:4,8:12),]))


Anova(CentrTAcc3)



#With H

CentrTAcc<-glmer(acc~NewH*Induction*IsPrimed*dur+(1|suj), family="binomial", droplevels(data2[data2$block %in% c(1,4) 
                                                                                                & data2$NewRaw %in% c(0:4,8:12),]))

CentrTAcc2<-glmer(acc~NewH*IsPrimed+(1|suj), family="binomial", droplevels(data2[data2$block %in% c(1,4)
                                                                                 & data2$Induction %in% c(1)
                                                                                 #& data2$dur %in% c(17)
                                                                                             & data2$NewH %in% c(1,3),]))
                                                                                          #& data2$NewRaw %in% c(0:12),]))
CentrTAcc2<-glmer(acc~NewH*IsPrimed+(1|suj), family="binomial", datatest)

#& data2$NewRaw %in% c(0:12),]))



CentrTAcc2<-glmer(acc~Induction*IsPrimed+(1|suj), family="binomial", droplevels(data2[data2$block %in% c(1,4)
                                                                                     #& data2$Induction %in% c(1)
                                                                                     #& data2$dur %in% c(17)
                                                                                     & data2$NewH %in% c(3),]))
#& data2$NewRaw %in% c(0:12),]))
#for the tukey in CentrACC, if you want to get 01 just do the contrast subsetting to lows alone!!!!!!#

CentrTAcc2<-glmer(acc~NewH+IsPrimed+NewH:IsPrimed+(1|suj), family="binomial", droplevels(data2[data2$block %in% c(1,4)
                                                                                               & data2$Induction %in% c(1)
                                                                                               #& data2$dur %in% c(17,33,67,84)
                                                                                               & data2$NewH %in% c(1,3),]))
                                                                                               #& data2$NewRaw %in% c(0:12),]))
Anova(CentrTAcc2)


lsm.options(pbkrtest.limit = 9557)
lsmeans(CentrTAcc2, pairwise ~ IsPrimed, adjust="tukey")


##BIC AND BF TOO! exp((BICfull - BICnull)/2) BIC > 1 favors BICnull / BIC < 1 favors BICfull

CentrTAcc2<-glmer(acc~NewH+IsPrimed+(1|suj), family="binomial", droplevels(data2[data2$block %in% c(1,4)
                                                                                 & data2$Induction %in% c(1)
                                                                                 #& data2$dur %in% c(17,33,67,84)
                                                                                 & data2$NewH %in% c(1,3),]))
                                                                              #& data2$NewRaw %in% c(0:12),]))

CentrTAcc2null<-glmer(acc~IsPrimed+(1|suj), family="binomial", droplevels(data2[data2$block %in% c(1,4)
                                                                                 & data2$Induction %in% c(1)
                                                                                 #& data2$dur %in% c(17,33,67,84)
                                                                                 & data2$NewH %in% c(1,3),]))
#& data2$NewRaw %in% c(0:12),]))

anova(CentrTAcc2, CentrTAcc2null)
exp((1411.8 - 1404.4)/2) #exp((BICfull - BICnull)/2)

CentrTAcc3<-glmer(acc~NewH*Induction*IsPrimed*dur*block+(1|suj), family="binomial", droplevels(data2[data2$block %in% c(1,4) 
                                                                                                       & data2$NewRaw %in% c(0:4,8:12),]))




AllModelsRawCentrAcc<-anova(CentrTAcc, CentrTAcc2,CentrTAcc3 )

Anova(CentrTAcc2)

e <- evalq(aggregate(list(acc=acc), list(suj=suj, NewRaw=NewRaw,
            IsPrimed=IsPrimed, Induction=Induction, NewH=NewH), mean),
           droplevels(data2[data2$block %in% c(1,4) 
                                                                                                                             & data2$NewRaw %in% c(0:4,8:12),]))

evalq(interaction.plot(NewH,IsPrimed,acc), e)



lsm.options(pbkrtest.limit = 9557)
lsmeans(CentrTAcc2, pairwise ~ IsPrimed*NewRaw, adjust="tukey")
lsmeans(CentrTAcc3, pairwise ~IsPrimed*NewRaw)

#With Induction

CentrTAcc<-glmer(acc~Induction*IsPrimed*dur+(1|suj), family="binomial", droplevels(data2[data2$block %in% c(1,4) 
                                                                                              & data2$NewH %in% c(3),]))

CentrTAcc2<-glmer(acc~Induction*IsPrimed+(1|suj), family="binomial", droplevels(data2[data2$block %in% c(1,4) 
                                                                                 & data2$NewH %in% c(3),]))
#& data2$NewRaw %in% c(0:12),]))

##BIC AND BF TOO! exp((BICfull - BICnull)/2) BIC > 1 favors BICnull / BIC < 1 favors BICfull

CentrTAcc2<-glmer(acc~Induction*IsPrimed+(1|suj), family="binomial", droplevels(data2[data2$block %in% c(1,4) 
                                                                                      & data2$NewH %in% c(3)
                                                                                      #& data2$Induction %in% c(0,1)
                                                                                      ,]))
CentrTAcc2null<-glmer(acc~IsPrimed+(1|suj), family="binomial", droplevels(data2[data2$block %in% c(1,4) 
                                                                                      & data2$NewH %in% c(3),]))

anova(CentrTAcc2, CentrTAcc2null)
exp((984.55 - 967.98)/2)

CentrTAcc3<-glmer(acc~Induction*IsPrimed*dur*block+(1|suj), family="binomial", droplevels(data2[data2$block %in% c(1,4) 
                                                                                                     & data2$NewH %in% c(3),]))




AllModelsRawCentrAcc<-anova(CentrTAcc, CentrTAcc2,CentrTAcc3 )

Anova(CentrTAcc2)


##LATEST VERSION HERE ^^




CentrTAcc<-glmer(acc~NewH*Induction*IsPrimed+(1|suj), family="binomial", droplevels(data2[data2$block %in% c(1,4) 
                                                                                              & data2$NewH %in% c(1,3),]))

CentrTAcc2<-glmer(acc~hyp*dur*IsPrimed+(1|suj), family="binomial", droplevels(data2[data2$block %in% c(1,4) 
                                                                              & data2$NewH %in% c(1,3),]))
CentrTAcc3<-glmer(acc~hyp*IsPrimed+(1|suj), family="binomial", droplevels(data2[data2$block %in% c(1,4) 
                                                                          & data2$hyp %in% c(1,2),]))

#for suggestion-only condition (mediums)

CentrTAcc<-glmer(acc~dur*IsPrimed*block+(1|suj), family="binomial", droplevels(data2[data2$block %in% c(1,4) 
                                                                                         & data2$hyp %in% c(3),]))
CentrTAcc2<-glmer(acc~dur*IsPrimed+(1|suj), family="binomial", droplevels(data2[data2$block %in% c(1,4) 
                                                                                    & data2$hyp %in% c(3),]))
CentrTAcc3<-glmer(acc~IsPrimed+(1|suj), family="binomial", droplevels(data2[data2$block %in% c(1,4) 
                                                                                & data2$hyp %in% c(3),]))

#With no induction vs induction

CentrTAcc<-glmer(acc~Induction*dur*IsPrimed*block+(1|suj), family="binomial", droplevels(data2[data2$block %in% c(1,4) 
                                                                                         & data2$Induction %in% c(0,1),]))
CentrTAcc2<-glmer(acc~Induction*dur*IsPrimed+(1|suj), family="binomial", droplevels(data2[data2$block %in% c(1,4) 
                                                                                    & data2$Induction %in% c(0,1),]))
CentrTAcc3<-glmer(acc~Induction*IsPrimed+(1|suj), family="binomial", droplevels(data2[data2$block %in% c(1,4) 
                                                                                & data2$Induction %in% c(0,1),]))



AllModelsCentrACC<-anova(CentrTAcc,CentrTAcc2,CentrTAcc3)

#CentrTAcc3<-glmer(acc~HarSc*HypnoYN*IsPrimed+(1|suj), family="binomial", data=data2
#                  [data2$block %in% c(1,4),])

CentrTAcc3 <- glmer(acc~IsPrimed+(1|suj), family='binomial',
                    droplevels(data2[data2$block %in% c(1,4)& data2$Induction %in% c(0,1),]))

AllModelsCentrACC<-anova(CentrTAcc,CentrTAcc2,CentrTAcc3)

Anova(CentrTAcc3)
summary(CentrTAcc3)


lsm.options(pbkrtest.limit = 9557)
lsmeans(CentrTAcc2, pairwise ~ IsPrimed*hyp, adjust="tukey")
lsmeans(CentrTAcc3, pairwise ~IsPrimed*hyp)

#Plots

#e <- evalq(aggregate(list(acc=acc), list(suj=suj, dur=dur, IsPrimed=IsPrimed, hyp=hyp), mean),
#                       data2[data2$block %in% c(1,4),])

#With no induction vs induction
e <- evalq(aggregate(list(acc=acc), list(suj=suj, dur=dur, IsPrimed=IsPrimed, Induction=Induction, NewH=NewH, NewRaw=NewRaw), mean),
           droplevels(data2[data2$block %in% c(1,4) & data2$NewH %in% c(1,3),]))

e <- evalq(aggregate(list(acc=acc), list(suj=suj, IsPrimed=IsPrimed, Induction=Induction, NewH=NewH), mean),droplevels(data2[data2$block %in% c(1,4) 
                 & data2$NewRaw %in% c(0:4,8:12),]))

evalq(interaction.plot(NewH,IsPrimed,acc), e[e$Induction==0,])


e <- evalq(aggregate(list(acc=acc), list(suj=suj, IsPrimed=IsPrimed, hyp=hyp), mean),
           droplevels(data2[data2$block %in% c(1,4) & data2$hyp %in% c(1,2),]))

e <- evalq(aggregate(list(acc=acc), list(suj=suj, IsPrimed=IsPrimed, HarSc=HarSc), mean),
           data2[data2$block %in% c(1,4),])



## Plot with error bars

e <- evalq(aggregate(list(acc=acc), list(suj=suj, NewH=NewH, IsPrimed=IsPrimed), mean),
           droplevels(data2[data2$block %in% c(1,4) & data2$Induction %in% c(1) & data2$NewH %in% c(1,3),]))

grandMean <- mean(e$acc)
partMeans <-evalq(aggregate(list(partMean=acc), list(suj=suj), mean),
                  droplevels(data2[data2$block %in% c(1,4) & data2$Induction %in% c(0,1),]))

## partMeans are participants means

e <- merge(e, partMeans)
e$nacc <- e$acc-e$partMean+grandMean
## we merge them and remove partMeans from the cells

f <- evalq(aggregate(list(seACC=nacc), list(NewH=NewH, IsPrimed=IsPrimed), se), e)
f$seACC <- f$seACC*sqrt(4/3) # compute se and apply Morey's 2008 correction (I had a 2X2 design so M=4).
f$macc <- evalq(aggregate(list(macc=acc), list(NewH=NewH, IsPrimed=IsPrimed), mean), e)$macc

evalq(interaction.plot(NewH, IsPrimed, acc, ylim=c(0.9,1)) ,e)
j <- 0.0001
xs <- c(1-j,3-j,1+j,3+j)
y1=f$macc+f$seACC
y0=f$macc-f$seACC

arrows(x0=xs, y0=y0, x1=xs, y1=y1, length=0,
       angle=90,code=3, lty=1, lwd=2)

evalq(interaction.plot(Induction,IsPrimed,acc), e)

evalq(interaction.plot(hyp,IsPrimed,acc), e)
intxplot(acc ~ hyp, groups=IsPrimed, data=e, se=c(1,1,1,0.3,0.3,0.3),
         ylim=c(0,3))
intxplot(acc ~ hyp, groups=IsPrimed, data=e, se=c(3.27,4.47,3.21,4.37,3.41,4.34))

evalq(interaction.plot(dur, IsPrimed, acc),e[e$hyp==1,])
evalq(interaction.plot(dur, IsPrimed, acc),e[e$hyp==2,])
evalq(interaction.plot(dur, IsPrimed, acc),e[e$hyp==3,])

hyp1<-subset(e, hyp=="1",select=acc)
hyp3<-subset(e, hyp=="3",select=acc)
hyp2<-subset(e, hyp=="2",select=acc)
boxplot(cbind(hyp1,hyp2,hyp3), ylim=c(0.9,1))

#Block 4

CentrTAcc<-glmer(acc~hyp*dur*IsPrimed+(1|suj), family="binomial", data=data2[data2$block %in% c(4),])

qqnorm(residuals(CentrTAcc))
plot(fitted(CentrTAcc),residuals(CentrTAcc))

Anova(CentrTAcc)
lsm.options(pbkrtest.limit = 9557)
lsmeans(CentrTAcc, pairwise ~ hyp | dur*IsPrimed, adjust="tukey")

#Plots

e <- evalq(aggregate(list(acc=acc), list(dur=dur, IsPrimed=IsPrimed, hyp=hyp, block=block), mean),
           data2[data2$block %in% c(1, 4),])
hyp1<-subset(e, hyp=="1",select=acc)
hyp3<-subset(e, hyp=="3",select=acc)
hyp2<-subset(e, hyp=="2",select=acc)
boxplot(cbind(hyp1,hyp2,hyp3))

#Response Time

#Block 1

CentrTRT<-lmer(rt~hyp*dur*block*IsPrimed+(1|suj), data=droplevels(data2[data2$block %in% c(1,4) & data2$hyp %in% c(1,2),]))
CentrTRT2<-lmer(rt~hyp*dur*IsPrimed+(1|suj), data=droplevels(data2[data2$block %in% c(1,4) & data2$hyp %in% c(1,2),]))
CentrTRT3<-lmer(rt~hyp*IsPrimed+(1|suj), data=droplevels(data2[data2$block %in% c(1,4) & data2$hyp %in% c(1,2),]))

AllModelsCentrTRT<-anova(CentrTRT,CentrTRT2,CentrTRT3)

CentrTRT<-lmer(rt~hyp*dur*IsPrimed+(1|suj), data=droplevels(data2[data2$block %in% c(1) & data2$hyp %in% c(1,2),]))
CentrTRT2<-lmer(rt~hyp*IsPrimed+(1|suj), data=data2[data2$block %in% c(1),])
CentrTRT<-lmer(rt~hyp*dur*IsPrimed+(1|suj), data=data2[data2$block %in% c(4),])
CentrTRT2<-lmer(rt~hyp*IsPrimed+(1|suj), data=data2[data2$block %in% c(4),])

anova(CentrTRT,CentrTRT2)

qqnorm(residuals(CentrTRT))
plot(fitted(CentrTRT),residuals(CentrTRT))

Anova(CentrTRT2)
lsm.options(pbkrtest.limit = 9557)
lsmeans(CentrTRT, pairwise ~ hyp | dur*rep, adjust="tukey")

#Plots
e <- evalq(aggregate(list(rt=rt), list(dur=dur, IsPrimed=IsPrimed, hyp=hyp, block=block), mean),
           data2[data2$block %in% c(1),])
hyp1<-subset(e, hyp=="1",select=rt)
hyp3<-subset(e, hyp=="3",select=rt)
hyp2<-subset(e, hyp=="2",select=rt)
boxplot(cbind(hyp1,hyp2,hyp3))

#Block 4

CentrTRT<-lmer(rt~hyp*dur*IsPrimed+(1|suj), data=droplevels(data2[data2$block %in% c(4) & data2$hyp %in% c(1,2),]))

qqnorm(residuals(CentrTRT))
plot(fitted(CentrTRT),residuals(CentrTRT))

Anova(CentrTRT)
lsm.options(pbkrtest.limit = 9557)
lsmeans(CentrTRT, pairwise ~ dur | hyp*IsPrimed, adjust="tukey")

#Plots

e <- evalq(aggregate(list(rt=rt), list(dur=dur, IsPrimed=IsPrimed, hyp=hyp, block=block), mean),
           data2[data2$block %in% c(4),])
evalq(interaction.plot(dur,IsPrimed,rt), e)
evalq(plot(dur,rt, main="Central task Response Time Across Peripheral Target Duration", xlab="Duration (ms)", ylab="Response Times (ms)"), e)
hyp1<-subset(e, hyp=="1",select=rt)
hyp3<-subset(e, hyp=="3",select=rt)
hyp2<-subset(e, hyp=="2",select=rt)
boxplot(cbind(hyp1,hyp2,hyp3))

######################
#Subjective Visibility
######################

#Accuracy

##LATEST VERSION HERE VV

##With Raw

SubjVisAcc <- lmer(Visibility~NewRaw*Induction*IsPrimed*dur*block+(1|suj), 
      droplevels(data2[data2$block %in% c(2,4) 
           #& data2$NewRaw %in% c(0:4,8:12),]))
           & data2$NewRaw %in% c(0:12),]))

           
SubjVisAcc2 <- lmer(Visibility~Induction*IsPrimed*dur+(1|suj), 
                   droplevels(data2[data2$block %in% c(2,4) 
                                    & data2$NewRaw %in% c(9:12),]))

SubjVisAcc3 <- lmer(Visibility~NewRaw*Induction*IsPrimed+(1|suj), 
                   droplevels(data2[data2$block %in% c(2,4) 
                                    & data2$NewRaw %in% c(0:4,8:12),]))

##With H

SubjVisAcc <- lmer(Visibility~NewH*Induction*IsPrimed*dur*block+(1|suj), 
                   droplevels(data2[data2$block %in% c(2,4) 
                                    & data2$NewRaw %in% c(0:4,8:12),]))

SubjVisAcc2 <- lmer(Visibility~NewH*IsPrimed*dur+(1|suj), 
                    droplevels(data2[data2$block %in% c(2,4)
                                     #& data2$dur %in% c(17,33,67,84)
                                     & data2$Induction %in% c(1)
                                     & data2$NewH %in% c(1,3) ,]))
                                     #& data2$NewRaw %in% c(0:12),]))

#reviewer question: will it work just for 33 vs 67#
SubjVisAcc2 <- lmer(Visibility~NewH*dur+(1|suj), 
                    droplevels(data[data$block %in% c(2,4)
                                     & data$dur %in% c(0,17)
                                     & data$Induction %in% c(1)
                                     & data$NewH %in% c(1,3) ,]))

SubjVisAcc2bis <- lmer(Visibility~NewH+(1|suj), 
                    droplevels(data[data$block %in% c(2,4)
                                    & data$dur %in% c(0,17)
                                    & data$Induction %in% c(1)
                                    & data$NewH %in% c(1,3) ,]))

anova(SubjVisAcc2, SubjVisAcc2bis)
#Compare Low Vs. High without Induction#
SubjVisAcc2 <- lmer(Visibility~AllH*dur+(1|suj), 
                    droplevels(data[data$block %in% c(2,4)
                                    & data$dur %in% c(17,33,67,84)
                                    #& data$Induction %in% c(0,1)
                                    & data$AllH %in% c(3,4) 
                                    ,]))

SubjVisAcc2 <- lmer(Visibility~AllH*dur*IsPrimed+(1|suj), 
                    droplevels(data2[data2$block %in% c(2,4)
                                    & data2$dur %in% c(17,33,67,84)
                                    #& data$Induction %in% c(0,1)
                                    & data2$AllH %in% c(4,5) 
                                    ,]))
SubjVisAcc2null <- lmer(Visibility~dur+(1|suj), 
                    droplevels(data2[data2$block %in% c(2,4)
                                    & data2$dur %in% c(17,33,67,84)
                                    #& data$Induction %in% c(0,1)
                                    & data2$AllH %in% c(4,5) 
                                    ,]))

anova(SubjVisAcc2, SubjVisAcc2null)

#& data2$NewRaw %in% c(0:12),]))


                                ##0 and 17 are the same##
SubjVisAcc2 <- lmer(Visibility~dur*NewH+(1|suj), 
                    droplevels(data[data$block %in% c(2,4)
                                     & data$dur %in% c(0,33)
                                     & data$Induction %in% c(1)
                                     & data$NewH %in% c(1,3) ,]))
                                ##0 and 17 are the same##



Anova(SubjVisAcc2)



SubjVisAcc3 <- lmer(Visibility~NewH*Induction*IsPrimed+(1|suj), 
                    droplevels(data2[data2$block %in% c(2,4) 
                                     & data2$NewRaw %in% c(0:4,8:12),]))


AllModelsRawSubjVis<-anova(SubjVisAcc, SubjVisAcc2, SubjVisAcc3)

Anova(SubjVisAcc2)

##With Induction

SubjVisAcc <- lmer(Visibility~Induction*IsPrimed*dur+(1|suj), 
                   droplevels(data2[data2$block %in% c(2,4) 
                                    & data2$NewH %in% c(3),]))

SubjVisAcc2 <- lmer(Visibility~NewH*Induction*dur+(1|suj), 
                    droplevels(data[data$block %in% c(2,4)
                                     & data$NewH %in% c(1,3)
                                     & data$Induction %in% c(0,1)
                                     & data$dur %in% c(17,33)
                                     ,]))

SubjVisAcc2 <- lmer(Visibility~dur*Induction+(1|suj), 
                    droplevels(data2[data2$block %in% c(2,4)
                                    & data2$NewH %in% c(3)
                                    & data2$NewRaw %in% c(9,10)
                                    #& data2$Induction %in% c(0)
                                    #& data2$dur %in% c(17,33)
                                    ,]))

##Low With Induction vs High without Induction



SubjVisAcc2 <- lmer(Visibility~AllH*dur*IsPrimed+(1|suj), 
                    droplevels(data2[data2$block %in% c(2,4)
                                     & data2$AllH %in% c(4,5)
                                     #& data2$NewRaw %in% c(9,10)
                                     #& data2$Induction %in% c(0)
                                     #& data2$dur %in% c(17,33)
                                     ,]))



SubjVisAcc3 <- lmer(Visibility~Induction*IsPrimed+(1|suj), 
                    droplevels(data2[data2$block %in% c(2,4) 
                                     & data2$NewH %in% c(3),]))


AllModelsRawSubjVis<-anova(SubjVisAcc, SubjVisAcc2, SubjVisAcc3)

Anova(SubjVisAcc)

e <- evalq(aggregate(list(Visibility=Visibility), 
                     list(suj=suj, IsPrimed=IsPrimed, NewH=NewH,
                            NewRaw=NewRaw,dur=dur),
                            mean),droplevels(data2[data2$block %in% c(2,4) 
                            & data2$NewRaw %in% c(0:4,8:12),]))

evalq(interaction.plot(dur, Induction, Visibility), e)
evalq(interaction.plot(dur, IsPrimed, Visibility), e)
evalq(interaction.plot(dur, NewH, Visibility), e[e$Induction==0,])
evalq(plot(IsPrimed, Visibility), e)





lsm.options(pbkrtest.limit = 9557)
lsmeans(SubjVisAcc2, pairwise ~ NewH*dur, adjust="tukey")



##LATEST VERSION HERE ^^




SubjVisAcc <- lmer(Visibility~NewH*Induction*dur*IsPrimed+(1|suj),  droplevels(data2[data2$block %in% c(2,4) & data2$NewH %in% c(1,3),]))


#SubjVisAcc <- glmer(Visibility~hyp*dur*block*IsPrimed+(1|suj), family="poisson", data2[data2$block %in% c(2,4),])
#SubjVisAcc2 <- glmer(Visibility~hyp*dur*IsPrimed+(1|suj), family="poisson", data2[data2$block %in% c(2,4),])
#SubjVisAcc3 <- glmer(Visibility~hyp*IsPrimed+(1|suj), family="poisson", data2[data2$block %in% c(2,4),])

#for highs and lows i+s
SubjVisAcc <- lmer(Visibility~hyp*dur*block*IsPrimed+(1|suj),  droplevels(data2[data2$block %in% c(2,4) & data2$hyp %in% c(1,2),]))
SubjVisAcc2 <- lmer(Visibility~hyp*dur*IsPrimed+(1|suj), droplevels(data2[data2$block %in% c(2,4) & data2$hyp %in% c(1,2),]))
SubjVisAcc3 <- lmer(Visibility~hyp*IsPrimed+(1|suj), droplevels(data2[data2$block %in% c(2,4) & data2$hyp %in% c(1,2),]))

#for suggestion-only condition (mediums)
SubjVisAcc <- lmer(Visibility~dur*block*IsPrimed+(1|suj),  droplevels(data2[data2$block %in% c(2,4) & data2$hyp %in% c(3),]))
SubjVisAcc2 <- lmer(Visibility~dur*IsPrimed+(1|suj), droplevels(data2[data2$block %in% c(2,4) & data2$hyp %in% c(3),]))
SubjVisAcc3 <- lmer(Visibility~IsPrimed+(1|suj), droplevels(data2[data2$block %in% c(2,4) & data2$hyp %in% c(3),]))


#With no induction vs induction

SubjVisAcc <- lmer(Visibility~Induction*dur*block*IsPrimed+(1|suj),  droplevels(data2[data2$block %in% c(2,4) & data2$Induction %in% c(0,1),]))
SubjVisAcc2 <- lmer(Visibility~Induction*dur*IsPrimed+(1|suj), droplevels(data2[data2$block %in% c(2,4) & data2$Induction %in% c(0,1),]))
SubjVisAcc3 <- lmer(Visibility~Induction*IsPrimed+(1|suj), droplevels(data2[data2$block %in% c(2,4) & data2$Induction %in% c(0,1),]))


AllSubjVisAcc<-anova(SubjVisAcc, SubjVisAcc2, SubjVisAcc3)


SubjVisAcc2 <- lmer(Visibility~hyp*dur*IsPrimed+(1|suj), droplevels(data2[data2$block %in% c(2,4) & data2$hyp %in% c(1,2),]))
Anova(SubjVisAcc2)


summary(SubjVisAcc2)
summary(SubjVisAcc3)
print(summary(SubjVisAcc3), correlation=TRUE)
Anova(SubjVisAcc2)

lsm.options(pbkrtest.limit = 9557)
#lsmeans(SubjVisAcc2, pairwise ~ hypTOT|dur, adjust="tukey")
lsmeans(SubjVisAcc3, pairwise ~ hyp*dur, adjust="tukey")
lsmeans(SubjVisAcc3, pairwise ~ hyp:dur)


#Visibility plots and tables

e <- evalq(aggregate(list(Visibility=Visibility), list(suj=suj,dur=dur, IsPrimed=IsPrimed, 
                              hyp=hyp, block=block), mean), data2[data2$block %in% c(2,4),])

e <- evalq(aggregate(list(Visibility=Visibility), list(suj=suj,Subli=Subli, NewH=NewH
                                                       ), mean), data2[data2$block %in% c(2,4),])

e <- evalq(aggregate(list(Visibility=Visibility), list(suj=suj,dur=dur, IsPrimed=IsPrimed, 
                            HypnoYN=HypnoYN, HarSc=HarSc, block=block), mean), 
           data[data$block %in% c(2,4),])

#With no induction vs induction
e <- evalq(aggregate(list(Visibility=Visibility), list(suj=suj,dur=dur, IsPrimed=IsPrimed, 
                             Induction=Induction), mean), 
           droplevels(data2[data2$block %in% c(2,4) & data2$Induction %in% c(0,1),]))


## Plot with error bars

e <- evalq(aggregate(list(Visibility=Visibility), list(suj=suj, NewH=NewH, dur=dur), mean),
           droplevels(data[data$block %in% c(2,4) & data$NewH %in% c(1,3) & data$Induction %in% c(1),]))


e <- evalq(aggregate(list(Visibility=Visibility), list(suj=suj, AllH=AllH, dur=dur), mean),
           droplevels(data2[data2$block %in% c(2,4) & data2$AllH %in% c(4,5),]))

grandMean <- mean(e$Visibility)
#grandMean <- evalq(aggregate(list(Visibility=Visibility), list(dur=dur, suj=suj), mean),e)
partMeans <-evalq(aggregate(list(partMean=Visibility), list(suj=suj), mean),
                  droplevels(data2[data2$block %in% c(2,4) & data2$AllH %in% c(4,5),]))

## partMeans are participants means

e <- merge(e, partMeans)
e$nVisibility <- e$Visibility-e$partMean+grandMean
## we merge them and remove partMeans from the cells

f <- evalq(aggregate(list(seVisibility=nVisibility), list(AllH=AllH, dur=dur), se), e)
f$seVisibility <- f$seVisibility*sqrt(6/4) # compute se and apply Morey's 2008 correction (I had a 2X2 design so M=4).
f$mVisibility <- evalq(aggregate(list(mVisibility=Visibility), list(AllH=AllH, dur=dur), mean), e)$mVisibility

evalq(interaction.plot(dur, AllH, Visibility), e)
j <- 0.0001
xs <- c(17-j,17+j, 33-j,33+j, 67-j,67+j, 84-j,84+j)
y1=f$mVisibility+f$seVisibility
y0=f$mVisibility-f$seVisibility

arrows(x0=xs, y0=y0, x1=xs, y1=y1, length=0,
       angle=90,code=3, lty=1, lwd=2)

evalq(interaction.plot(Induction,IsPrimed,Visibility), e)


evalq(interaction.plot(HarSc,dur,Visibility), e[e$HypnoYN=="yes" & e$HarSc<=9 & e$HarSc>=2,])
evalq(interaction.plot(HarSc,dur,Visibility), e[e$HypnoYN=="no",])

evalq(interaction.plot(dur, hyp, Visibility), e)

evalq(interaction.plot(dur,HypnoYN,Visibility), e)

e <-droplevels(e[e$block %in% c(2,4),])

e <- evalq(aggregate(list(Visibility=Visibility), 
                     list(suj=suj,dur=dur, IsPrimed=IsPrimed, 
                    hyp=hyp), mean), data2[data2$block %in% c(4),])

evalq(interaction.plot(hyp,IsPrimed,Visibility), e)

intxplot(Visibility ~ dur, groups=IsPrimed, data=e, se=TRUE,
         ylim=c(1,3.5))

evalq(interaction.plot(hyp,block, Visibility), e)

e <- evalq(aggregate(list(IsSeen=IsSeen), list(dur=dur, IsPrimed=IsPrimed, hyp=hyp, block=block), mean),
           data2[data2$block %in% c(2,4),])

#Response Time

SubjVisRT <- lmer(RTimeVis~hyp*dur*block*IsPrimed+(1|suj), data2[data2$block %in% c(2,4),])

SubjVisRT2 <- lmer(RTimeVis~hyp*dur*IsPrimed+(1|suj), data2[data2$block %in% c(2,4),])

SubjVisRT <- lmer(RTimeVis~hyp*dur*IsPrimed+(1|suj), data2[data2$block %in% c(2),])
SubjVisRT2 <- lmer(RTimeVis~hyp*IsPrimed+(1|suj), data2[data2$block %in% c(2),])
SubjVisRT <- lmer(RTimeVis~hyp*dur*IsPrimed+(1|suj), data2[data2$block %in% c(4),])
SubjVisRT2 <- lmer(RTimeVis~hyp*IsPrimed+(1|suj), data2[data2$block %in% c(4),])


Anova(SubjVisRT)

anova(SubjVisRT,SubjVisRT2)

e <- evalq(aggregate(list(RTimeVis=RTimeVis), list(suj=suj,dur=dur, IsPrimed=IsPrimed, 
                                                       hyp=hyp, block=block), mean), 
           data[data$block %in% c(2,4),])

e <-droplevels(e[e$block %in% c(2,4),])

evalq(interaction.plot(hyp,block, RTimeVis), e)

Anova(SubjVisAcc2)
lsm.options(pbkrtest.limit = 9557)
lsmeans(SubjVisAcc2, pairwise ~ hyp| dur*IsPrimed, adjust="tukey")
lsmeans(SubjVisAcc2, pairwise ~ hyp*block | dur*IsPrimed, adjust="tukey")


#Response Time plots and tables

e <- evalq(aggregate(list(RTimeVis=RTimeVis), list(suj=suj, dur=dur, IsPrimed=IsPrimed, hyp=hyp), mean),
           data2[data2$block %in% c(2,4),])
evalq(interaction.plot(dur, hyp, RTimeVis), e)

e <- evalq(aggregate(list(RTimeVis=RTimeVis), list(suj=suj, dur=dur, IsPrimed=IsPrimed, hyp=hyp, block=block), mean),
           data2[data2$block %in% c(2,4),])
evalq(interaction.plot(hyp, block, RTimeVis), e)

e <- evalq(aggregate(list(RTimeVis=RTimeVis), list(dur=dur, IsSeen=IsSeen, IsPrimed=IsPrimed, hyp=hyp, block=block), mean),
           data2[data2$block %in% c(2,4) & data2$hyp==3,])
evalq(interaction.plot(dur, IsSeen, RTimeVis), e)


#######################
#Objective Visibility
######################

##LATEST VERSION HERE VV

##With Raw

ObjVisAcc <- glmer(acc~Induction*IsPrimed*dur+(1|suj), family="binomial", 
                   droplevels(data2[data2$block %in% c(3) 
                                    & data2$NewRaw %in% c(8:12),]))

ObjVisAcc2 <- glmer(acc~Induction*dur+(1|suj), family="binomial", 
                    droplevels(data2[data2$block %in% c(3) 
                                     & data2$NewRaw %in% c(8:12),]))


Anova(ObjVisAcc2)

##With H

ObjVisAcc <- glmer(acc~NewH*dur*IsPrimed+(1+IsPrimed|suj), family="binomial", 
                   droplevels(data2[data2$block %in% c(3) 
                                    & data2$NewH %in% c(1,3)
                                    #& data2$dur %in% c(17,33,67,84)
                                    & data2$Induction %in% c(1),]))



ObjVisAcc <- glmer(acc~Induction*IsPrimed+(1|suj), family="binomial", 
                   droplevels(data2[data2$block %in% c(3) 
                                    & data2$NewH %in% c(3)
                                    & data2$dur %in% c(17)
                                    & data2$Induction %in% c(0,1)
                                   ,]))

ObjVisAcc2 <- glmer(acc~NewH*IsPrimed+(1|suj), family="binomial", 
                    droplevels(data2[data2$block %in% c(3) 
                                     & data2$NewRaw %in% c(1,3),]))




AllModelsRawObjVis<-anova(ObjVisAcc, ObjVisAcc2)

Anova(ObjVisAcc)

##With Induction

ObjVisAcc <- glmer(acc~Induction*IsPrimed*dur+(1|suj), family="binomial", 
                   droplevels(data2[data2$block %in% c(3) 
                                    & data2$NewH %in% c(3)
                                    & data2$dur %in% c(17,33,67,84)
                                    ,]))

ObjVisAcc <- glmer(acc~Induction*IsPrimed*Subli+(1|suj), family="binomial", 
                   droplevels(data2[data2$block %in% c(3) 
                                    & data2$NewH %in% c(3)
                                    & data2$Induction %in% c(0,1)
                                    ,]))


ObjVisAcc2 <- glmer(acc~NewH*IsPrimed*dur+(1|suj), family="binomial", 
                    droplevels(data2[data2$block %in% c(3) 
                                     & data2$NewH %in% c(1,3)
                                     & data2$Induction %in% c(1),]))

ObjVisAcc2null <- glmer(acc~IsPrimed*Subli+(1|suj), family="binomial", 
                    droplevels(data2[data2$block %in% c(3) 
                                     & data2$NewH %in% c(3),]))

##BIC AND BF TOO! exp((BICfull - BICnull)/2) BIC > 1 favors BICnull / BIC < 1 favors BICfull

anova(ObjVisAcc2, ObjVisAcc2null)


AllModelsRawObjVis<-anova(ObjVisAcc, ObjVisAcc2)

Anova(ObjVisAcc)
exp((2345.5 - 2332.9)/2)

e <- evalq(aggregate(list(acc=acc), 
                     list(suj=suj, IsPrimed=IsPrimed,
                          NewH=NewH, dur=dur),
                     mean),droplevels(data2[data2$block %in% c(3) 
                                            & data2$NewH %in% c(1,3)
                                            & data2$Induction %in% c(1)
                                            ,]))

evalq(interaction.plot(dur, IsPrimed, acc), e)
evalq(interaction.plot(IsPrimed, Induction, acc), e)
evalq(interaction.plot(dur, NewH, acc), e[e$IsPrimed==1,])


evalq(interaction.plot(dur, NewH, acc), e[e$Induction==1,])
evalq(interaction.plot(dur, NewH, acc), e[e$Induction==0,])

evalq(interaction.plot(NewH, IsPrimed, acc), e[e$Induction==1,])
evalq(interaction.plot(NewH, IsPrimed, acc), e[e$Induction==0,])


lsm.options(pbkrtest.limit = 9557)
lsmeans(ObjVisAcc2, pairwise ~ NewH*IsPrimed, adjust="tukey")
lsmeans(ObjVisAcc, pairwise ~ NewRaw*IsPrimed, adjust="tukey")
lsmeans(ObjVisAcc, pairwise ~ Induction*IsPrimed, adjust="tukey")


##LATEST VERSION HERE ^^

#Accuracy

#for highs and lows i+s
ObjVisAcc <- glmer(acc~hyp*dur*IsPrimed+(1|suj), family='binomial', droplevels(data2[data2$block %in% c(3) & data2$hyp %in% c(1,2),]))
ObjVisAcc2 <- glmer(acc~hyp*IsPrimed+(1|suj), family='binomial', droplevels(data2[data2$block %in% c(3) & data2$hyp %in% c(1,2),]))

#for suggestion-only condition (mediums)
ObjVisAcc <- glmer(acc~dur*IsPrimed+(1|suj), family='binomial', droplevels(data2[data2$block %in% c(3) & data2$hyp %in% c(3),]))
ObjVisAcc2 <- glmer(acc~IsPrimed+(1|suj), family='binomial', droplevels(data2[data2$block %in% c(3) & data2$hyp %in% c(3),]))

#With no induction vs induction

ObjVisAcc <- glmer(acc~Induction*dur*IsPrimed+(1|suj), family='binomial', droplevels(data2[data2$block %in% c(3) & data2$Induction %in% c(0,1),]))
ObjVisAcc2 <- glmer(acc~Induction*IsPrimed+(1|suj), family='binomial', droplevels(data2[data2$block %in% c(3) & data2$Induction %in% c(0,1),]))


AllModelsObjVisAcc<-anova(ObjVisAcc,ObjVisAcc2)

ObjVisAcc <- glmer(acc~dur*IsPrimed+(1|suj), family='binomial', droplevels(data2[data2$block %in% c(3) & data2$hyp %in% c(1),]))

ObjVisAcc <- glmer(acc~NewH*Induction*dur*IsPrimed+(1|suj), family='binomial', droplevels(data2[data2$block %in% c(3) & data2$NewH %in% c(1,3),]))
Anova(ObjVisAcc)
Anova(ObjVisAcc2)

lsm.options(pbkrtest.limit = 9557)
lsmeans(ObjVisAcc, pairwise ~ Induction*IsPrimed, adjust="tukey")
lsmeans(ObjVisAcc, pairwise ~ hyp*block | dur*IsPrimed, adjust="tukey")



#Accuracy Plots and Tables

e <- evalq(aggregate(list(acc=acc), list(suj=suj, NewH=NewH, dur=dur, IsPrimed=IsPrimed, Induction=Induction), mean),
           droplevels(data2[data2$block %in% c(3) & data2$Induction %in% c(0,1) & data2$NewH %in% c(3),]))

#With no induction vs induction
e <- evalq(aggregate(list(acc=acc), list(suj=suj, Induction=Induction, dur=dur, IsPrimed=IsPrimed), mean),
           droplevels(data2[data2$block %in% c(3) & data2$Induction %in% c(0,1),]))






e <- evalq(aggregate(list(acc=acc), list(suj=suj, dur=dur, IsPrimed=IsPrimed, hyp=hyp), mean),
           droplevels(data2[data2$block %in% c(3) & data2$hyp %in% c(2,3),]))
e <- evalq(aggregate(list(acc=acc), list(suj=suj, dur=dur, IsPrimed=IsPrimed, hyp=hyp), mean),
           data2[data2$block %in% c(3)& data2$hyp %in% c(1,2),])
e <- evalq(aggregate(list(acc=acc), list(suj=suj, IsPrimed=IsPrimed, NewH=NewH, Induction=Induction), mean),
           droplevels(data2[data2$block %in% c(3) 
                 & data2$dur %in% c(17,33,67,84)
                 & data2$NewH %in% c(1,3)
                 & data2$Induction %in% c(1)
                 ,]))
evalq(interaction.plot(hyp, IsPrimed, acc, fixed = TRUE, col = 2:4), e)
evalq(interaction.plot(NewH, IsPrimed, acc), e[e$Induction==0,])
evalq(interaction.plot(IsPrimed,Induction, acc), e)
evalq(interaction.plot(NewH, IsPrimed, acc), e)
evalq(interaction.plot(dur, IsPrimed, acc), e[e$hyp==1,])
evalq(interaction.plot(dur, IsPrimed, acc), e[e$hyp==2,])
evalq(interaction.plot(dur, IsPrimed, acc), e[e$hyp==3,])
evalq(interaction.plot(dur, IsPrimed, acc), e)
evalq(interaction.plot(dur, hyp, acc), e)



## Plot with error bars

e <- evalq(aggregate(list(acc=acc), list(suj=suj, NewH=NewH, IsPrimed=IsPrimed, dur=dur), mean),
           droplevels(data2[data2$block %in% c(3) & data2$NewH %in% c(1,3) 
                            & data2$Induction %in% c(1)
                            ,]))

e <- evalq(aggregate(list(acc=acc), list(suj=suj, NewH=NewH, dur=dur, IsPrimed=IsPrimed, Induction=Induction), mean),
           droplevels(data2[data2$block %in% c(3) & data2$Induction %in% c(0,1) & data2$NewH %in% c(3),]))


grandMean <- mean(e$acc)
#grandMean <- evalq(aggregate(list(Visibility=Visibility), list(dur=dur, suj=suj), mean),e)
partMeans <-evalq(aggregate(list(partMean=acc), list(suj=suj), mean),
                  droplevels(data2[data2$block %in% c(3) & data2$Induction %in% c(0,1)
                                   & data2$NewH %in% c(3)
                                   ,]))

## partMeans are participants means

e <- merge(e, partMeans)
e$nacc <- e$acc-e$partMean+grandMean
## we merge them and remove partMeans from the cells

f <- evalq(aggregate(list(seacc=nacc), list(Induction=Induction, dur=dur), se), e)
f$seacc <- f$seacc*sqrt(8/7) # compute se and apply Morey's 2008 correction (I had a 4X2 design so M=8).
f$macc <- evalq(aggregate(list(macc=acc), list(Induction=Induction, dur=dur), mean), e)$macc


evalq(interaction.plot(dur, Induction, acc, ylim=c(0.4,1) ), e)
j <- 0.0001
xs <- c(17-j,17+j, 33-j,33+j, 67-j,67+j, 84-j,84+j)
y1=f$macc+f$seacc
y0=f$macc-f$seacc

arrows(x0=xs, y0=y0, x1=xs, y1=y1, length=0,
       angle=90,code=3, lty=1, lwd=2)





evalq(interaction.plot(dur,hyp,acc), e[e$IsPrimed==1,])
evalq(lines(dur,hyp),e[e$IsPrimed==0,]) 
#evalq(interaction.plot(dur,hyp,acc), e[e$IsPrimed==0,])

evalq(interaction.plot(hyp, IsPrimed, perf, ylim=c(600,620)), e)
intxplot(perf ~ hyp, groups=IsPrimed, data=e, se=TRUE,
         ylim=c(500,720))
par(new=T)

intxplot(acc ~ dur, groups=hyp, data=e[e$IsPrimed==0,], se=TRUE,
         ylim=c(0.1,1))

intxplot(acc ~ IsPrimed, groups=hyp, data=e, se=TRUE,
         ylim=c(0.4,0.84), xlim=c(0.5,2.5))
lsmeans(ObjVisAcc, pairwise ~ hyp|IsPrimed, adjust="tukey")


#Response Time
ObjVisRT <- lmer(rt~hyp*dur*IsPrimed+(1|suj), droplevels(data2[data2$block %in% c(3) & data2$hyp %in% c(1,2),]))
ObjVisRT2 <- lmer(rt~hyp*IsPrimed+(1|suj), droplevels(data2[data2$block %in% c(3) & data2$hyp %in% c(1,2),]))

AllMOdelsObjVisRT<-anova(ObjVisRT,ObjVisRT2)



ObjVisRT <- lmer(rt~hyp*IsPrimed*dur+(1|suj), droplevels(data2[data2$block %in% c(3),]))
ObjVisRT2 <- lmer(rt~hyp*IsPrimed+(1|suj), droplevels(data2[data2$block %in% c(3),]))


Anova(ObjVisRT2)
lsm.options(pbkrtest.limit = 4901)
lsmeans(ObjVisRT, pairwise ~ dur|hyp, adjust="tukey")

#Response Time Plots and tables
e <- evalq(aggregate(list(perf=perf), list(dur=dur, IsPrimed=IsPrimed, hyp=hyp), mean),
           data2[data2$block %in% c(3),])

evalq(interaction.plot(dur, hyp, perf, ylim=c(530,700), col=c(2,3,4),
                       ylab="RTs", xlab="Peripheral Target Duration"),e)
evalq(interaction.plot(dur, hyp, perf, ylim=c(530,700), col=c(2,3,4),
                       ylab="RTs", xlab="Peripheral Target Duration"), e[e$IsPrimed==0,])


evalq(interaction.plot(dur, hyp, perf, ylim=c(530,700), col=c(2,3,4), legend=FALSE, type = "l", 
                       ylab="RTs", xlab="Peripheral Target Duration"), e[e$IsPrimed==0,])
par(new=T)
evalq(interaction.plot(dur, hyp, perf, ylim=c(530,700), col=c(2,3,4), axes=FALSE, xlab=' ', 
                       ylab=' ', legend=FALSE, type = "l"), e[e$IsPrimed==1,])
legend("topright", "No Hypnosis: Blue, Low: Green, High: Red")
evalq(interaction.plot(dur, hyp, perf), e[e$IsPrimed==1,])
evalq(interaction.plot(IsPrimed, hyp, perf), e)


#for plotting interaction plots with error bars
e <- evalq(aggregate(list(rt=rt), list(suj=suj, swt=swt, ttype=ttype), mean),
           df[evalq(acc==1, df),])
## swt and ttype are my two conditions here

e <- evalq(aggregate(list(acc=acc), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed), mean),
           data2[data2$block %in% c(1,4),])

grandMean <- mean(e$acc)
partMeans <-evalq(aggregate(list(partMean=acc), list(suj=suj), mean),
                  data2[data2$block %in% c(1,4),])

## partMeans are participants means

e <- merge(e, partMeans)
e$nacc <- e$acc-e$partMean+grandMean
## we merge them and remove partMeans from the cells

f <- evalq(aggregate(list(seACC=nacc), list(hyp=hyp, IsPrimed=IsPrimed), se), e)
f$seACC <- f$seACC*sqrt(6/3) # compute se and apply Morey's 2008 correction (I had a 2X2 design so M=4).
f$macc <- evalq(aggregate(list(macc=acc), list(hyp=hyp, IsPrimed=IsPrimed), mean), e)$macc

evalq(interaction.plot(hyp, IsPrimed, acc, ylim=c(0.93,0.98)) ,e)
j <- 0.0001
xs <- c(1-j,2-j,3-j,1+j,2+j,3+j)
y1=f$macc+f$seACC
y0=f$macc-f$seACC

arrows(x0=xs, y0=y0, x1=xs, y1=y1, length=0,
       angle=90,code=3, lty=1, lwd=2)




## swt and ttype are my two conditions here

e <- evalq(aggregate(list(Visibility=Visibility), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed, dur=dur), mean),
           data2[data2$block %in% c(2,4),])

grandMean <- mean(e$Visibility)
partMeans <-evalq(aggregate(list(partMean=Visibility), list(suj=suj), mean),
                  data2[data2$block %in% c(2,4),])

## partMeans are participants means

e <- merge(e, partMeans)
e$nVisibility <- e$Visibility-e$partMean+grandMean
## we merge them and remove partMeans from the cells

f <- evalq(aggregate(list(seVisibility=nVisibility), list(hyp=hyp,IsPrimed=IsPrimed), se), e)
f$seVisibility <- f$seVisibility*sqrt(6/3) # compute se and apply Morey's 2008 correction (I had a 2X2 design so M=4).
f$mVisibility <- evalq(aggregate(list(mVisibility=Visibility), list(hyp=hyp,IsPrimed=IsPrimed), mean), e)$mVisibility

evalq(interaction.plot(hyp, IsPrimed, Visibility, ylim=c(1.6,2.6)) ,e)
j <- 0.01
xs <- c(1-j,2-j,3-j,1+j,2+j,3+j)
y1=f$mVisibility+f$seVisibility
y0=f$mVisibility-f$seVisibility

arrows(x0=xs, y0=y0, x1=xs, y1=y1, length=0,
       angle=90,code=3, lty=1, lwd=2)

## swt and ttype are my two conditions here

e <- evalq(aggregate(list(acc=acc), list(suj=suj, IsPrimed=IsPrimed, hyp=hyp), mean),
           data2[data2$block %in% c(3),])

grandMean <- mean(e$acc)
partMeans <-evalq(aggregate(list(partMean=acc), list(suj=suj), mean),
                  data2[data2$block %in% c(3),])

## partMeans are participants means

e <- merge(e, partMeans)
e$nacc <- e$acc-e$partMean+grandMean
## we merge them and remove partMeans from the cells

f <- evalq(aggregate(list(seacc=nacc), list(IsPrimed=IsPrimed,hyp=hyp), se), e)
f$seacc <- f$seacc*sqrt(6/3) # compute se and apply Morey's 2008 correction (I had a 2X2 design so M=4).
f$macc <- evalq(aggregate(list(macc=acc), list(IsPrimed=IsPrimed,hyp=hyp), mean), e)$macc

evalq(interaction.plot(IsPrimed, hyp, acc, ylim=c(0.45,0.88)) ,e)
j <- 0.001
xs <- c(1-j,2-j,1+j,2+j)
y1=f$macc+f$seacc
y0=f$macc-f$seacc

arrows(x0=xs, y0=y0, x1=xs, y1=y1, length=0,
       angle=90,code=3, lty=1, lwd=2)