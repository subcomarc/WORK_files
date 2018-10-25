##Analysis Script for the HypnoBlink Experiment
##By subcomarc, with the resources of Jérôme Sackur
## Needs file "blinkData.R"

############### ################################
#Preloading settings, directories and libraries#
############### ################################

setwd("C:/Users/lscpuser/Google Drive/R/")

source("C:/Users/lscpuser/Google Drive/R/R_Lib_Sackur.r")
#setwd("C:/Users/subcomarc/Google Drive/R/")

#source("C:/Users/subcomarc/Google Drive/R/R_Lib_Sackur.r")
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

load("blinkData.r")
data$hypnoC[data$hypnoC==4]<-0 # If you are going to plot, so that no hypno is next to low

####################### ########################################
#Prepping the data for stats and adding some convenient columns#
###################### ######################################### ##################

data <- read.table(file="C:/Users/lscpuser/Google Drive/R/BlinkTable.txt", header=T)
#data <- read.table(file="C:/Users/subcomarc/Google Drive/R/BlinkTable.txt", header=T)

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
data$hypnoC[data$hypnoC==0]<-4
data$lag <- data$AllLags
data$soa <- data$AllSOAs
data$lag <- ordered(data$lag)


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

save(data, file='blinkData.R')



###################### #########################
#All contrasts and dataframes for linear models#
##################### ##########################
data$target1 <- factor(data$target1)
data$target2 <- factor(data$target2)
data$hypnoC <- factor(data$hypnoC)
#data$suj <- factor(data$suj)
data$lag <- factor(data$lag)
data$lagC <- factor(data$lagC)
data$lagComplete <- factor(data$lagComplete)
data$lagCatComplete <- factor(data$lagCatComplete)

contrasts(data$target1)=contr.sum(2)
contrasts(data$target2)=contr.sum(2)
contrasts(data$hypnoC)=contr.sum(4)
contrasts(data$lag)=contr.sum(6)
contrasts(data$lagC)=contr.sum(2)
contrasts(data$lagComplete)=contr.sum(7)
contrasts(data$lagCatComplete)=contr.sum(3)



###################### ###################### ######################
#STATISTICAL ANALYSES#
###################### ###################### ######################


#Remember: accOr= Accuracy Orientation Task. accn= Accuracy Number of Faces Task.
#rtOr= RT orientation task. rtn= RT number of faces task. hypnoC= hypnosis category (no hypno 4, or 0 
#if you plotting)


#Accuracy Inclination task
##########################

AccInc <- glmer(accOr~lagC*target1*target2*hypnoC+(1|suj), family='binomial', data=data)
#qqnorm(residuals(AccInc))
#plot(fitted(AccInc),residuals(AccInc))
Anova(AccInc)
lsm.options(pbkrtest.limit = 9557)
lsmeans(AccInc, pairwise ~ hypnoC*target2, adjust="tukey")

e <- evalq(aggregate(list(accOr=accOr), list(suj=suj, hypnoC=hypnoC, target1=target1, target2=target2, lagC=lagC), 
                     mean),data)
evalq(interaction.plot(hypnoC, lagC, accOr), e)


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
evalq(interaction.plot(target2, hypnoC, accOr), e)



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



#Attentional Blink????????
##########################

e <- evalq(aggregate(list(accOr=accOr), list(suj=suj, hypnoC=hypnoC, target1=target1, target2=target2, lagC=lagC), 
                     mean),data)
evalq(interaction.plot(lagC, target2, accOr), e[e$hypnoC==3,])
evalq(interaction.plot(lagC, hypnoC, accOr), e)



#DIFFUSION MODEL (load the file Fast-dm_diffusion_parameters.R if it exists)
######### #######

## First we plot densities


densLCorr <- density(data$rtOr[data$accOr==1 & data$hypnoC==1])
densMCorr <- density(data$rtOr[data$accOr==1 & data$hypnoC==2])
densHCorr <- density(data$rtOr[data$accOr==1 & data$hypnoC==3])
densNHCorr <- density(data$rtOr[data$accOr==1 & data$hypnoC==4])
densLInc <- density(data$rtOr[data$accOr==0 & data$hypnoC==1])
densMInc <- density(data$rtOr[data$accOr==0 & data$hypnoC==2])
densHInc <- density(data$rtOr[data$accOr==0 & data$hypnoC==3])
densNHInc <- density(data$rtOr[data$accOr==0 & data$hypnoC==4])

M <- max(c(densLCorr$y,densMCorr$y,densHCorr$y,densNHCorr$y, densLInc$y,densMInc$y,densHInc$y,
           densNHInc$y))
#pdf("Densities.pdf", width = 12, height = 6)
plot(densLCorr, col=2, lwd=3, ylim=c(0, M), xlim=c(0, 3500))
lines(densMCorr, col=1, lwd=3); lines(densHCorr, col=3, lwd=3); lines(densNHCorr, col=4, lwd=3) 
lines(densLInc, col=2, lwd=1, lty=2); lines(densMInc, col=1, lwd=1, lty=2); lines(densHInc, col=3, lwd=1, lty=2)
lines(densNHInc, col=4, lwd=1, lty=2)
#dev.off()

###################### ######################
sujs <- unique(data$suj);
i <- 1
df <- data
## ## We recode t1 et t2 levels, because it makes things easier
## ## aftewards.
df$T1 <- paste(df$target1, "t1", sep='')
df[df$nface==1,]$T1 <- 'absent'
df$T2 <- paste(df$target2, "t2", sep='')

## ## We write the data for each subject in different files, in different
## ## folder so as to parallelize by hand.
 for (s in sujs){
#     folder <- paste("C:/Users/lscpuser/Google Drive/R/DM", ceiling(i/12), "/", sep='')
     folder <- "C:/Users/lscpuser/Google Drive/R/DM/"
     temp <- evalq(data.frame(acc=accOr, rt=rtOr/1000, target1=T1,
                              target2=T2, lag=lagC),
                   df[df$suj==s,])
     write.table(temp, file=paste(folder, s, ".txt", sep=''),
                 row.names=FALSE, col.names=FALSE, quote=FALSE)
     i <- i+1
 }
 #fdm <- do.call(rbind, lapply(1:5, function(n){
  #            fdmImp <- readFastDMComplete(file=paste("diffusion/d", n,
              fdmImp <- readFastDMComplete(file="C:/Users/lscpuser/Google Drive/R/DM/blinkFinal.txt",
                                          dvs=c('v','a','t0'),
                                          factors=c('target1', 'target2', 'lag'),
                                          flevels=list(c('absent', 'angryt1', 'neutralt1'),
                                          c('angryt2', 'neutralt2'),
                                          c('early', 'late')))
#         }))


 f <- evalq(aggregate(list(hypnoC=hypnoC), list(suj=suj), unique), data[data$suj>300,])
 #fdm <- merge(fdm, f) #loading it from before because I already have it and it's a waste of time otherwise
 fdmImp <-merge(fdmImp,f)
 fdm <-rbind(fdm,fdmImp) #and now let's bring them all together like a happy family
 save(fdm, file="Fast-dm_diffusion_parameters.R")
###################### ###################### 

load("Fast-dm_diffusion_parameters.R")

 fdm$target1 <- factor(fdm$target1)
 fdm$target2 <- factor(fdm$target2)
 fdm$hypnoC <- factor(fdm$hypnoC)
 fdm$lag <- factor(fdm$lag)

 contrasts(fdm$target1)=contr.sum(3)
 contrasts(fdm$target2)=contr.sum(2)
 contrasts(fdm$hypnoC)=contr.sum(4)
 contrasts(fdm$lag)=contr.sum(2)
 
#Check if you can use lm with the DM parameters (Spoiler: YES YOU CAN)
### ##### 
Test<-lmer(v~hypnoC*target1*target2*lag+(1|suj), data=fdm)
plot(fitted(Test),residuals(Test)) # Residual plot doesn't show patterns, linearity is ok 
qqnorm(residuals(Test)) # Residuals are normal/linear
fdm$v<-as.numeric(as.character(fdm$v)) #Transforming to check for homoskedacity
fdm$logv<-log(fdm$v); fdm$logv[fdm$logv==NaN]<-0.001
Test<-lmer(logv~hypnoC*target1*target2*lag+(1|suj), data=fdm)
plot(fitted(Test),residuals(Test))#Homoskedacity seems ok enough after log transform

Test<-lmer(a~hypnoC*target1*target2*lag+(1|suj), data=fdm)
plot(fitted(Test),residuals(Test)) # Residual plot doesn't show patterns, linearity is ok 
qqnorm(residuals(Test)) # Residuals are normal/linear
fdm$a<-as.numeric(as.character(fdm$a)) #Transforming to check for homoskedacity
fdm$loga<-log(fdm$a); fdm$logv[fdm$loga==NaN]<-0.001
Test<-lmer(loga~hypnoC*target1*target2*lag+(1|suj), data=fdm)
plot(fitted(Test),residuals(Test))#Homoskedacity seems ok after log transform

Test<-lmer(t0~hypnoC*target1*target2*lag+(1|suj), data=fdm)
plot(fitted(Test),residuals(Test)) # Residual plot doesn't show patterns, linearity is ok 
qqnorm(residuals(Test)) # Residuals are normal/linear
fdm$v<-as.numeric(as.character(fdm$t0)) #Transforming to check for homoskedacity
fdm$logv<-log(fdm$t0); fdm$logv[fdm$logt0==NaN]<-0.001
Test<-lmer(logt0~hypnoC*target1*target2*lag+(1|suj), data=fdm)
plot(fitted(Test),residuals(Test))#Homoskedacity seems ok enough after log transform
