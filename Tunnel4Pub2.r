##Analysis Script for the HypnoTunnel Experiment
##By subcomarc, with the resources of Jérôme Sackur
## Needs file "TunnelData.r"

################################################
#Preloading settings, directories and libraries#
################################################

#setwd("C:/Users/lscpuser/Google Drive/R/")
#data <- read.table(file="C:/Users/lscpuser/Google Drive/R/DataTable RAW.txt", header=T)
#source("C:/Users/lscpuser/Google Drive/R/R_Lib_Sackur.r")
setwd("C:/Users/subcomarc/Google Drive/R/")
#data <- read.table(file="C:/Users/subcomarc/Google Drive/R/DataTable RAW.txt", header=T)
source("C:/Users/subcomarc/Google Drive/R/R_Lib_Sackur.r")
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

##################################################################
#Load preprocessed data (for preprocessing see analysesTUNNEL2.r)#
##################################################################

load("TunnelData.r")
#load("TunnelData_2.r")


###################### ###################### ######################
#STATISTICAL ANALYSES#
###################### ###################### ######################

#Prepping the data for stats and adding some convenient columns

data<-data
data$rep <- 1 #for repetition priming
data$rep[data$Target!=data$Prime] <- 0
data$perf <- data$rt #for having rt for all blocks on the same column
data$perf[data$block==2] <- data$rtVis[data$block==2]
data$RTimeVis<-round(data$RTimeVis*1000)
data$RTimeVis[data$block %in% c(1,3)]<-5000
data <- droplevels(data[data$RTimeVis<5001,])
data<- droplevels(data[!(data$block!=2 & data$rt>1200),])
data$IsSeen<-1 #for treating visibility as a binary
data$IsSeen[data$Visibility<2]<-0 #for treating visibility as a binary
#data$EoLate<-0
#data$EoLate[data$dur>33]<-1

data2 <- data[data$dur>0 & data$PrimeNoPrime==1 & data$IsPrimed!=2,] #to have a dataset without the trials without prime

data$block <- factor(data$block)
data$hyp <- factor(data$hyp)
data$dur <- factor(data$dur)
data$Quadrant <- factor(data$Quadrant)
data$IsPrimed <- factor(data$IsPrimed)
data$rep<- factor(data$rep)
#data$EoLate<-factor(data$EoLate)
data2$block <- factor(data2$block)
data2$hyp <- factor(data2$hyp)
data2$dur <- factor(data2$dur)
data2$Quadrant <- factor(data2$Quadrant)
data2$IsPrimed <- factor(data2$IsPrimed)
data2$rep<- factor(data2$rep)
#data2$EoLate<-factor(data2$EoLate)

data3<-data2[data2$dur>33]

contrasts(data$Quadrant)=contr.sum(4)
contrasts(data$block)=contr.sum(4)
contrasts(data$hyp)=contr.sum(3)
contrasts(data$dur)=contr.sum(5)
contrasts(data$rep)=contr.sum(2)
contrasts(data$IsPrimed)=contr.sum(3)
#contrasts(data$EoLate)=contr.sum(2)
contrasts(data2$Quadrant)=contr.sum(4)
contrasts(data2$block)=contr.sum(4)
contrasts(data2$hyp)=contr.sum(3)
contrasts(data2$dur)=contr.sum(4)
contrasts(data2$rep)=contr.sum(2)
contrasts(data2$IsPrimed)=contr.sum(2)
#contrasts(data2$EoLate)=contr.sum(2)

data2<-droplevels(data2[!(data2$hyp==2 & data2$HarSc==5),])
data2<-droplevels(data2[!(data2$hyp==1 & data2$HarSc==7),])

control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))

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

#Block 1

CentrTAcc<-glmer(acc~hyp*dur*IsPrimed*block+(1|suj), family="binomial", data=data2[data2$block %in% c(1,4) 
                                                                             ,])
CentrTAcc2<-glmer(acc~hyp*dur*IsPrimed+(1|suj), family="binomial", data=data2[data2$block %in% c(1,4) 
                                                                             ,])
CentrTAcc3<-glmer(acc~hyp*IsPrimed+(1|suj), family="binomial", data=data2[data2$block %in% c(1,4) 
                                                                               ,])
CentrTAcc<-glmer(acc~hyp*dur*IsPrimed+(1|suj), family="binomial", data=data2[data2$block %in% c(1),])
CentrTAcc2<-glmer(acc~hyp*dur*IsPrimed+(1|suj), family="binomial", data=data2[data2$block %in% c(1) 
                                                                          ,])

CentrTAcc2<-glmer(acc~hyp*dur*IsPrimed+(1|suj), family="binomial", data=data2[data2$block %in% c(1) 
                                                                                     ,])
CentrTAcc3<-glmer(acc~hyp*dur+(1|suj), family="binomial", data=data2[data2$block %in% c(1) 
                                                                          ,])

CentrTAcc<-glmer(acc~HarSc*dur*IsPrimed+(1|suj), family="binomial", data=data2[data2$block %in% c(1,4),])
                                                                             


Anova(CentrTAcc3)
summary(CentrTAcc2)

AllModelsCentrACC<-anova(CentrTAcc,CentrTAcc2,CentrTAcc3)

CentrTAcc<-glmer(acc~hyp*dur*IsPrimed+(1|suj), family="binomial", data=droplevels(data2[data2$block %in% c(4),]))
CentrTAcc2<-glmer(acc~hyp*IsPrimed+(1|suj), family="binomial", data=data2[data2$block %in% c(4) 
                                                                          ,])


test1<-glmer(acc~hyp*dur*block*IsPrimed+(1|suj), family="binomial", data=data2[data2$block %in% c(1,4),])
test4<-glmer(acc~hyp*block*IsPrimed+(1|suj), family="binomial", data=data2[data2$block %in% c(1,4),])
test2<-glmer(acc~hyp*dur*IsPrimed+(1|suj), family="binomial", data=data2[data2$block %in% c(1,4),])
test3<-glmer(acc~hyp*IsPrimed+(1|suj), family="binomial", data=data2[data2$block %in% c(1,4),])
test5<-glmer(acc~hyp*block*dur+(1|suj), family="binomial", data=data2[data2$block %in% c(1,4),])
test6<-glmer(acc~hyp*block+(1|suj), family="binomial", data=data2[data2$block %in% c(1,4),])
test7<-glmer(acc~hyp*dur+(1|suj), family="binomial", data=data2[data2$block %in% c(1,4),])


anova(test1,test2,test3,test4)

anova(test1,test2,test3,test4,test5,test6,test7)


anova(CentrTAcc2,CentrTAcc)

BIC(CentrTAcc); BIC(CentrTAcc2)
BF=exp(-.5*(BIC(CentrTAcc)-BIC(CentrTAcc2)));BF

data2<-data2[data2$block==1 | data2$block==4,]
CentrTAcc<-glmer(acc~IsPrimed*hyp+(1|suj), family="binomial", data=data2)

CentrTAcc<-glmer(acc~hyp*IsPrimed+(1|suj), family="binomial", data=DF)

Anova(CentrTAcc2)
summary(CentrTAcc2)

contrasts(data3$dur)=contr.sum(2)

qqnorm(residuals(CentrTAcc))
plot(fitted(CentrTAcc),residuals(CentrTAcc))

 lsm.options(pbkrtest.limit = 9557)
lsmeans(CentrTAcc2, pairwise ~ IsPrimed|hyp, adjust="tukey")
lsmeans(CentrTAcc, pairwise ~IsPrimed*hyp)


Anova(CentrTAcc)

#Plots

e <- evalq(aggregate(list(acc=acc), list(suj=suj, dur=dur, IsPrimed=IsPrimed, hyp=hyp), mean),
                       data2[data2$block %in% c(1,4),])
e <- evalq(aggregate(list(acc=acc), list(suj=suj, IsPrimed=IsPrimed, hyp=hyp), mean),
           data2[data2$block %in% c(1,4),])

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

CentrTRT<-lmer(rt~hyp*dur*block*IsPrimed+(1|suj), data=data2[data2$block %in% c(1,4),])
CentrTRT2<-lmer(rt~hyp*dur*IsPrimed+(1|suj), data=data2[data2$block %in% c(4),])
CentrTRT3<-lmer(rt~hyp*IsPrimed+(1|suj), data=data2[data2$block %in% c(1,4),])

AllModelsCentrTRT<-anova(CentrTRT,CentrTRT2,CentrTRT3)

CentrTRT<-lmer(rt~hyp*dur*IsPrimed+(1|suj), data=data2[data2$block %in% c(1),])
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

CentrTRT<-lmer(rt~hyp*dur*IsPrimed+(1|suj), data=data2[data2$block %in% c(4),])

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

SubjVisAcc <- glmer(Visibility~hyp*dur*block*IsPrimed+(1|suj), family="poisson", data2[data2$block %in% c(2,4),])
SubjVisAcc2 <- glmer(Visibility~hyp*dur*IsPrimed+(1|suj), family="poisson", data2[data2$block %in% c(2,4),])
SubjVisAcc3 <- glmer(Visibility~hyp*IsPrimed+(1|suj), family="poisson", data2[data2$block %in% c(2,4),])

AllModelsSubjVis<-anova(SubjVisAcc,SubjVisAcc2,SubjVisAcc3)

SubjVisAcc <- lmer(Visibility~hyp*dur*block*IsPrimed+(1|suj), data2[data2$block %in% c(2,4),])

SubjVisAcc2 <- lmer(Visibility~hyp*dur*IsPrimed+(1|suj), data2[data2$block %in% c(2),])


SubjVisAcc <- lmer(Visibility~hyp*dur*IsPrimed+(1|suj), data2[data2$block %in% c(2),])
SubjVisAcc2 <- lmer(Visibility~hyp*IsPrimed+(1|suj), data2[data2$block %in% c(2),])
SubjVisAcc <- lmer(Visibility~hyp*dur*IsPrimed+(1|suj), data2[data2$block %in% c(4),])
SubjVisAcc2 <- lmer(Visibility~hyp*IsPrimed+(1|suj), data2[data2$block %in% c(4),])

test1<-lmer(Visibility~hyp*dur*block*IsPrimed+(1|suj), data2[data2$block %in% c(2,4),])
test2<-lmer(Visibility~hyp*dur*IsPrimed+(1|suj), data2[data2$block %in% c(2,4),])
test3<-lmer(Visibility~hyp*IsPrimed+(1|suj), data2[data2$block %in% c(2,4),])

Anova(test2)

SubjVisAcc2 <- lmer(Visibility~hyp*dur*IsPrimed+(1|suj), data2[data2$block %in% c(4),])

summary(SubjVisAcc2)
Anova(SubjVisAcc2)
anova(SubjVisAcc,SubjVisAcc2)

BIC(SubjVisAcc2); BIC(SubjVisAcc)
BF=exp(-.5*(BIC(SubjVisAcc)-BIC(SubjVisAcc2)));BF
qqnorm(residuals(SubjVisAcc))
plot(fitted(SubjVisAcc),residuals(SubjVisAcc))

Anova(SubjVisAcc2)
lsm.options(pbkrtest.limit = 9557)
lsmeans(SubjVisAcc2, pairwise ~ hyp|dur*IsPrimed, adjust="tukey")

#(Visibility as binary)
SubjVisAcc <- glmer(IsSeen~hyp*dur*IsPrimed+(1|suj), family = "binomial", data2[data2$block %in% c(2,4),])
SubjVisAcc2 <- glmer(IsSeen~hyp*dur*IsPrimed*block+(1|suj), family = "binomial", data2[data2$block %in% c(2,4),])

BIC(SubjVisAcc2); BIC(SubjVisAcc)
BF=exp(-.5*(BIC(SubjVisAcc)-BIC(SubjVisAcc2)));BF
qqnorm(residuals(SubjVisAcc))
plot(fitted(SubjVisAcc),residuals(SubjVisAcc))

Anova(SubjVisAcc)
lsm.options(pbkrtest.limit = 9557)
lsmeans(SubjVisAcc2, pairwise ~ hyp | dur*IsPrimed, adjust="tukey")


#Accuracy plots and tables

e <- evalq(aggregate(list(Visibility=Visibility), list(suj=suj,dur=dur, IsPrimed=IsPrimed, 
                              hyp=hyp, block=block), mean), data[data$block %in% c(2,4),])

e <-droplevels(e[e$block %in% c(2,4),])

e <- evalq(aggregate(list(Visibility=Visibility), 
                     list(suj=suj,dur=dur, IsPrimed=IsPrimed, 
                    hyp=hyp), mean), data2[data2$block %in% c(4),])

evalq(interaction.plot(hyp,IsPrimed,Visibility), e)

intxplot(Visibility ~ dur, groups=IsPrimed, data=e, se=TRUE,
         ylim=c(1,3.5))


evalq(interaction.plot(hyp,block, Visibility), e)

evalq(interaction.plot(dur, hyp, Visibility), e[e$IsPrimed!=2,])

e <- evalq(aggregate(list(Visibility=Visibility), list(suj=suj, dur=dur, IsPrimed=IsPrimed, hyp=hyp), mean),
           data2[data2$block %in% c(2,4),])
evalq(interaction.plot(hyp, IsPrimed, Visibility), e[e$IsPrimed!=2,])

e <- evalq(aggregate(list(IsSeen=IsSeen), list(dur=dur, IsPrimed=IsPrimed, hyp=hyp, block=block), mean),
           data2[data2$block %in% c(2,4),])
evalq(interaction.plot(dur, hyp, IsSeen), e)

#Response Time

SubjVisRT <- lmer(RTimeVis~hyp*dur*block*IsPrimed+(1|suj), data2[data2$block %in% c(2,4),])

Anova(SubjVisRT)

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

BIC(SubjVisRT); BIC(SubjVisRT2)
BF=exp(-.5*(BIC(SubjVisRT)-BIC(SubjVisRT2)));BF
qqnorm(residuals(SubjVisRT2))
plot(fitted(SubjVisRT),residuals(SubjVisRT))

Anova(SubjVisAcc2)
lsm.options(pbkrtest.limit = 9557)
lsmeans(SubjVisAcc2, pairwise ~ hyp| dur*IsPrimed, adjust="tukey")
lsmeans(SubjVisAcc2, pairwise ~ hyp*block | dur*IsPrimed, adjust="tukey")

#(RT Visibility as binary)

#Response Time plots and tables

e <- evalq(aggregate(list(RTimeVis=RTimeVis), list(dur=dur, IsPrimed=IsPrimed, hyp=hyp), mean),
           data2[data2$block %in% c(2,4),])
evalq(interaction.plot(dur, hyp, RTimeVis), e)

e <- evalq(aggregate(list(RTimeVis=RTimeVis), list(dur=dur, IsPrimed=IsPrimed, hyp=hyp, block=block), mean),
           data2[data2$block %in% c(2,4),])
evalq(interaction.plot(hyp, block, RTimeVis), e)

e <- evalq(aggregate(list(RTimeVis=RTimeVis), list(dur=dur, IsSeen=IsSeen, IsPrimed=IsPrimed, hyp=hyp, block=block), mean),
           data2[data2$block %in% c(2,4) & data2$hyp==3,])
evalq(interaction.plot(dur, IsSeen, RTimeVis), e)


#######################
#Objective Visibility
######################


#Anova()
#summary(glht(model20, mcp(Sound="Tukey")))
#lsmeans(mymod, pairwise ~ var1 | var2)
#lsmeans(mymod, pairwise ~ var1 | var2)
#lsmip(lsm, factor1 ~ factor2, type = "response")
#plot(lsm, by = "factor2", intervals = TRUE, type = "response")


#Accuracy

ObjVisAcc <- glmer(acc~hyp*dur*IsPrimed+(1|suj), family='binomial', data2[data2$block %in% c(3),])
ObjVisAcc2 <- glmer(acc~hyp*IsPrimed+(1|suj), family='binomial', data2[data2$block %in% c(3),])

AllModelsObjVisAcc<-anova(ObjVisAcc,ObjVisAcc2)

Anova(ObjVisAcc)
lsmeans(ObjVisAcc, pairwise ~ IsPrimed|dur, adjust="tukey")

ObjVisAcc <- glmer(acc~hyp*dur*IsPrimed+(1|suj), family='binomial', data2[data2$block %in% c(3),])
ObjVisAcc2 <- glmer(acc~hyp*IsPrimed+(1|suj), family='binomial', droplevels(data2[data2$block %in% c(3),]))

anova(ObjVisAcc,ObjVisAcc2)

Anova(ObjVisAcc)

#Accuracy Plots and Tables

e <- evalq(aggregate(list(acc=acc), list(suj=suj, dur=dur, IsPrimed=IsPrimed, hyp=hyp), mean),
           data2[data2$block %in% c(3),])
e <- evalq(aggregate(list(perf=perf), list(suj=suj, IsPrimed=IsPrimed, hyp=hyp), mean),
           data2[data2$block %in% c(3),])
evalq(interaction.plot(IsPrimed, hyp, acc, fixed = TRUE, col = 2:4), e)
evalq(interaction.plot(hyp, IsPrimed, perf), e)
evalq(interaction.plot(dur, IsPrimed, acc), e[e$hyp==1,])
evalq(interaction.plot(dur, IsPrimed, acc), e[e$hyp==2,])
evalq(interaction.plot(dur, IsPrimed, acc), e[e$hyp==3,])
evalq(interaction.plot(dur, IsPrimed, acc), e)
evalq(interaction.plot(dur, hyp, acc), e)

evalq(interaction.plot(dur, hyp, acc, ylim=c(0.2,1), col=c(2,3,4), legend=FALSE, type = "l", 
                       ylab="Accuracy", xlab="Peripheral Target Duration"), e[e$IsPrimed==0,])
par(new=T)
evalq(interaction.plot(dur, hyp, acc, ylim=c(0.2,1), col=c(2,3,4), axes=FALSE, xlab=' ', 
                       ylab=' ', legend=FALSE, type = "l"), e[e$IsPrimed==1,])


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
ObjVisRT <- lmer(rt~hyp*dur*IsPrimed+(1|suj), data2[data2$block %in% c(3),])
ObjVisRT2 <- lmer(rt~hyp*IsPrimed+(1|suj), data2[data2$block %in% c(3),])


ObjVisRT <- lmer(rt~hyp*IsPrimed*dur+(1|suj), droplevels(data2[data2$block %in% c(3),]))
ObjVisRT2 <- lmer(rt~hyp*IsPrimed+(1|suj), droplevels(data2[data2$block %in% c(3),]))

AllMOdelsObjVisRT<-anova(ObjVisRT,ObjVisRT2)

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

########################
#Visibility Correlations
########################
e <- evalq(aggregate(list(acc=acc, Visibility=Visibility), list(suj=suj, dur=dur, IsPrimed=IsPrimed, hyp=hyp, block=block), mean),
           data2[data2$block %in% c(2, 3),])
hyp1<-subset(e[e$block==3,], hyp=="1",select=acc);Vhyp1<-subset(e[e$block==2,], hyp=="1",select=Visibility)
hyp3<-subset(e[e$block==3,], hyp=="3",select=acc);Vhyp3<-subset(e[e$block==2,], hyp=="3",select=Visibility)
hyp2<-subset(e[e$block==3,], hyp=="2",select=acc);Vhyp2<-subset(e[e$block==2,], hyp=="2",select=Visibility)
Dhyp1<-subset(e[e$block==3,], hyp=="1",select=dur);Dhyp2<-subset(e[e$block==3,], hyp=="2",select=dur)
Dhyp3<-subset(e[e$block==3,], hyp=="3",select=dur);
Shyp1<-subset(e[e$block==3,], hyp=="1",select=suj);Shyp2<-subset(e[e$block==3,], hyp=="2",select=suj)
Shyp3<-subset(e[e$block==3,], hyp=="3",select=suj);
Hhyp1<-subset(e[e$block==3,], hyp=="1",select=hyp);Hhyp2<-subset(e[e$block==3,], hyp=="2",select=hyp)
Hhyp3<-subset(e[e$block==3,], hyp=="3",select=hyp);



HiHyp<-as.matrix(cbind(hyp1,Vhyp1,Dhyp1,Shyp1));HiHyp<-as.data.frame(cbind(hyp1,Vhyp1,Dhyp1,Shyp1,Hhyp1))
LoHyp<-as.matrix(cbind(hyp2,Vhyp2,Dhyp2,Shyp2));LoHyp<-as.data.frame(cbind(hyp2,Vhyp2,Dhyp2,Shyp2,Hhyp2))
NoHyp<-as.matrix(cbind(hyp3,Vhyp3,Dhyp3,Shyp3));NoHyp<-as.data.frame(cbind(hyp3,Vhyp3,Dhyp3,Shyp3,Hhyp3))
AllHyp<-rbind(HiHyp,LoHyp,NoHyp)
AllHyp$Visibility<-log(AllHyp$Visibility)
AllHyp$hyp<-as.factor(AllHyp$hyp)


rcorr(HiHyp);rcorr(LoHyp);rcorr(NoHyp);rcorr(AllHyp)

data3<-data2
data3$Visibility<-factor(data3$Visibility)
ObjSubjCorr <- lmer(acc~Visibility*hyp+(1|suj),  AllHyp)
Anova(ObjSubjCorr)
lsmeans(ObjSubjCorr, pairwise ~ hyp, adjust="tukey")

scatterplot(Visibility ~ acc | hyp, data=AllHyp,
            xlab="Objective vis", ylab="Subjective Vis", 
            main="Scatter Plot all hypno", 
            labels=row.names(AllHyp))
plot(AllHyp$acc,AllHyp$Visibility, type="p",)


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