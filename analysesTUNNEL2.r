setwd("C:/Users/subcomarc/Google Drive/R/")
setwd("C:/Users/lscpuser/Google Drive/R/")
data <- read.table(file="C:/Users/lscpuser/Google Drive/R/DataTable RAW.txt", header=T)
source("C:/Users/lscpuser/Google Drive/R/R_Lib_Sackur.r")
library('xtable')
library('knitr')
library('lme4')
library('car')
library('multcomp')
options(xtable.floating = FALSE)
options(xtable.timestamp = "")


data <- read.table(file="DataTable RAW.txt", header=T)
data$suj <- factor(data$SubjName)
data$hyp <- floor(data$SubjName/1000)
data$rt <- round(numerize(data$RTime)*1000)
data$acc <- data$Correctness
data$vis <- data$Visibility
data$block <- data$BlockNum
data$dur <- factor(round(data$PrimeDur * 1000), ordered=TRUE)
data$rtVis <- round(numerize(data$RTimeVis)*1000)
data$centralTask <- 1; data$centralTask[data$block %in% c(2, 3)] <- 0
data$periphTask <- 0; data$periphTask[data$block==3] <- 1
data$periphVis <- 0; data$periphVis[data$block %in% c(2,4)] <- 1


# Erreur de codage
data$Hypnotiz[data$suj==3051] <- 3

sujs <- unique(data$suj)

## Be careful:

## hyp:
## 1: high hypnotizability
## 2: low  hypnotizability
## 3: medium (and here, without induction)

## Hypnotiz:
## 1: low
## 2: medium low
## 3: medium high
## 4: high



# Erreur de réponse: les sujets confondent les touches de visibilité
# et de catégorisation
data <- droplevels(data[data$vis!=8,])
data <- droplevels(data[data$acc!=8,])
data <- droplevels(data[data$rt<6999 | data$block==2,])

sum(data$acc[data$centralTask==1]>1)/nrow(data[data$centralTask==1,]) 

data <- droplevels(data[data$centralTask==0 | data$acc<2,])

e <- evalq(aggregate(list(medianRt=rt), list(suj=suj), median), data[data$rt < 10000 & data$centralTask==1,])
e$sd <- evalq(aggregate(rt, list(suj=suj), sd), data [data$rt < 10000 & data$centralTask==1,])$x

data <- merge(data, e)

data$outlier <- data$rt > data$medianRt + 3 * data$sd



e <- evalq(aggregate(list(outlier=outlier), list(suj=suj), mean), data[data$centralTask==1,])
any(e$outlier>.1); max(e$outlier)
data <- droplevels(data[data$centralTask==0 | (data$outlier==0 & data$rt >150) ,])
hist(data$rt[data$centralTask==1])

save(data, file="TunnelData.r")

####
load("TunnelData.r")

##########################################

# Blocs 1 et 4: tâche sur la cible centrale

e <- evalq(aggregate(list(acc=acc), list(hyp=hyp, bloc=BlockNum==1), mean),
           data[data$BlockNum %in% c(1,4),])

e <- evalq(aggregate(list(rt=rt), list(hyp=hyp, bloc=BlockNum==1), median),
           data[data$BlockNum %in% c(1,4),])

par(mfrow=c(2,2))

evalq(interaction.plot(bloc, hyp, acc, lwd=2, lty=1:3, pch=2:4, type='b'), e)
evalq(interaction.plot(bloc, hyp, rt, lwd=2, lty=1:3, pch=2:4, type='b'), e)
e$hyp<-factor(e$hyp)
e$bloc<-factor(e$bloc)
e$sd<-sd(e$acc)
e$sd<-sd(e$rt)
e$sd<-by(e$acc,e$hyp,sd)


###########
opar <- theme_update(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(colour = "black"))
gp <- ggplot(e, aes(bloc, acc, colour=hyp, group=hyp))
gp + geom_line(aes(linetype=hyp), size=.6) + 
  geom_point(aes(shape=hyp), size=3) + 
  geom_errorbar(aes(ymax=acc+sd, ymin=acc-sd), width=.1)
theme_set(opar)

opar <- theme_update(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(colour = "black"))
gp <- ggplot(e, aes(bloc, rt, colour=hyp, group=hyp))
gp + geom_line(aes(linetype=hyp), size=.6) + 
  geom_point(aes(shape=hyp), size=3) + 
  geom_errorbar(aes(ymax=rt+sd, ymin=rt-sd), width=.1)
theme_set(opar)
###########

evalq(plot(density(rt), col=3), data[evalq(hyp==3 & BlockNum==4, data),])
evalq(lines(density(rt)), data[evalq(hyp==1 & BlockNum==4, data),])
evalq(lines(density(rt), col=2), data[evalq(hyp==2 & BlockNum==4, data),])


e <- evalq(aggregate(list(rt=rt), list(suj=suj, hyp=hyp, bloc=BlockNum==1),
                     median),
           data[data$BlockNum %in% c(1,4),])

l <- aov(rt~hyp*bloc+Error(suj/(bloc)), data=e)
summary(l); read.table(l)

e <- evalq(aggregate(list(rt=rt), list(suj=suj, hyp=hyp==1),
                     sd),
           data[data$BlockNum ==4,])

l <- aov(rt~hyp, data=e)
evalq(t.test(rt[hyp], rt[!hyp]), e)

sujs <- levels(suj)
par(mfrow=c(4,6))
lapply(sujs, function(s){evalq(plot(density(rt), main=s),
                               data[evalq(suj==s & BlockNum!=2, data),])})


par(mfrow=c(2,2))

# Visibilité? Blocs 2 et 4

lapply(1:3, function(s){evalq(hist(Visibility, main=s, breaks=0:4, ylim=c(0, 1500)),
                               data[evalq(hyp==s & BlockNum %in% c(2,4) & Visibility!=8, data),])})

sum(data$Visibility[data$BlockNum==4]==8)

e <- evalq(aggregate(list(vis=vis), list(hyp=hyp, dur=dur), mean),
           data[data$block==2,])
evalq(interaction.plot(dur, hyp, vis), e)

e <- evalq(aggregate(list(vis=vis), list(hyp=Hypnotiz, dur=dur),
                     function(x){sum(x>1)}),
           data[data$block %in% c(2,4),])
evalq(interaction.plot(dur, hyp, vis, col=1:4, lty=1:4, lwd=2, pch=2:5, type='b'), e)

data$hiVis <- as.numeric(data$vis>1 )
data$loVis <- as.numeric(data$vis==1 )

library("lme4")

l <- glmer(loVis~hyp*dur+(1|suj), family='binomial', data[data$block %in% c(2,4),])
summary(l)




# Categorization on the prime

e <- evalq(aggregate(list(acc=acc), list(hyp=hyp, dur=dur), mean),
           data[data$block %in% 3,])
evalq(interaction.plot(dur, hyp, acc), e)

e <- evalq(aggregate(list(rt=rt), list(hyp=hyp, dur=dur), median),
           data[data$block %in% 3,])
evalq(interaction.plot(dur, hyp, rt, lwd=2, pch=2:4, lty=2:4, type='b'), e)

e <- evalq(aggregate(list(rt=rt), list(hyp=hyp, prime=Prime), median),
           data[data$block %in% 3,])
evalq(interaction.plot(prime, hyp, rt, lwd=2, pch=2:4, lty=2:4, type='b'), e)



library(lme4)

l <- glmer(acc~hyp*dur+(1|suj), family='binomial', data[data$block==3,])
summary(l)


# Priming? blocks: 1 et 4

# First do we have the distance effect?

e <- evalq(aggregate(list(rt=rt), list(hyp=hyp, target=Target), median),
           data[data$block %in% c(1,4),])
evalq(interaction.plot(target, hyp, rt, lwd=2, pch=2:4, lty=2:4, type='b'), e)

e <- evalq(aggregate(list(rt=rt), list(hyp=hyp, primed=IsPrimed), median),
           data[data$block %in% c(1,4) & data$IsPrimed!=2,])
evalq(interaction.plot(primed, hyp, rt, lwd=2, pch=2:4, lty=2:4, type='b'), e)


data$short <- data$dur %in% c(17,33)
e <- evalq(aggregate(list(rt=rt), list(dur=dur, primed=IsPrimed), median),
           data[data$block %in% c(1,4) & data$IsPrimed!=2,])
evalq(interaction.plot(primed, dur, rt, lwd=2, pch=2:4, lty=2:4, type='b'), e)



e <- evalq(aggregate(list(rt=rt), list(hyp=hyp, primed=IsPrimed), median),
           data[data$block %in% c(1,4) & data$IsPrimed!=2 & data$short,])
evalq(interaction.plot(primed, hyp, rt, lwd=2, pch=2:4, lty=2:4, type='b'), e)

e <- evalq(aggregate(list(rt=rt), list(vis=vis>1, primed=IsPrimed, hyp=hyp),
                     median),
           data[data$block %in% c(1,4) & data$IsPrimed!=2,])

evalq(interaction.plot(hyp, primed, rt, lwd=2, pch=2:4, lty=2:4, type='b',
                       main='Seen'),
      e[e$vis,])
evalq(interaction.plot(hyp, primed, rt, lwd=2, pch=2:4, lty=2:4, type='b',
                       main='Not seen'),
            e[!e$vis,])


# Effets sur les temps réponse en fonction de la visibilité subjective

e <- evalq(aggregate(list(rt=rt), list(hyp=hyp, hiVis=hiVis), median),
           data[data$BlockNum %in% c(1,4),])
evalq(interaction.plot(hiVis, hyp, rt, lwd=2, lty=1:3, pch=2:4, type='b'), e)

# Effets sur le temps de réponse sur la visibilité...

e <- evalq(aggregate(list(rtVis=rtVis), list(suj=suj, hyp=hyp, vis=vis), median),
           data[data$BlockNum %in% c(2,4),])
evalq(interaction.plot(vis, hyp, rtVis, lwd=2, lty=1:3, pch=2:4, type='b'), e)

l <- aov(rtVis~hyp*vis+Error(suj/(vis)), data=e); summary(l)


##########################################
##########################################

##Does Hypnosis affect Visibility ratings? YES, AND IN THE SENSE OF THE SUGGESTION

e <- evalq(aggregate(list(rt=rt, acc=acc, vis=vis), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed), mean),
           data[data$block %in% c(2, 4) & data$IsPrimed<2,]) # Only Block 2 (it's significant when I use mean instead)
e$hyp <- factor(e$hyp); e$IsPrimed <- factor(e$IsPrimed)
l <- aov(vis~hyp*IsPrimed+Error(suj/(IsPrimed)), data=e) #Because I don't test trials were there isn't a prime
summary(l)
par(mfrow=c(2,2))
evalq(hist(vis), data[data$block %in% c(2,4),])
evalq(hist(vis, main='high', freq=FALSE, ylim=c(0,3)),
      data[data$block %in% c(2,4) & data$hyp==1 & data$PrimeNoPrime==0,])
evalq(hist(vis, main='low', freq=FALSE, ylim=c(0,3)),
      data[data$block %in% c(2,4) & data$hyp==2 & data$PrimeNoPrime==0,])
evalq(hist(vis, main='noH', freq=FALSE, ylim=c(0,3)), data[data$block %in% c(2,4) & data$hyp==3 & data$PrimeNoPrime==0,])


library(lme4)

e <- evalq(aggregate(list(vis=vis), list(suj=suj, hyp=hyp, PrimeDur=PrimeDur), mean),
           data[data$block %in% c(2, 4) & data$IsPrimed<2,]) # Only Block 2 (it's significant when I use mean instead)

e$PrimeDur <- factor(round(e$PrimeDur*1000), ordered=TRUE)
e$hyp <- factor(e$hyp)

e$ratio <- e$notSeen/f$notSeen
l <- aov(vis~hyp*PrimeDur+Error(suj/(PrimeDur)), data=e)
l <- aov(ratio~hyp*PrimeDur+Error(suj/(PrimeDur)), data=e)
summary(l); model.tables(l, type='means')
evalq(interaction.plot(PrimeDur, hyp, vis), e)
l0 <- lm(vis~hyp*PrimeDur, data=e[e$hyp!=3,])
l1 <- lm(vis~PrimeDur, data=e[e$hyp!=3,])
BIC(l0); BIC(l1)
BF=exp(-.5*(BIC(l1)-BIC(l0)));BF


e$hyp<-factor(e$hyp)
e$IsPrimed<-factor(e$IsPrimed)
e <- evalq(aggregate(list(notSeen=vis), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed), function(x){sum(x==1)}),
           data[data$block %in% c(2, 4) & data$IsPrimed<2,]) # Only Block 2 (it's significant when I use mean instead)
f <- evalq(aggregate(list(notSeen=vis), list(suj=suj, hyp=hyp), function(x){sum(x<5)}),
           data[data$block %in% c(2, 4) & data$IsPrimed<2,]) # Only Block 2 (it's significant when I use mean instead)
e$ratio <- e$notSeen/f$notSeen
l <- aov(ratio~hyp*IsPrimed, data=e)
xtable(l)
summary(l); model.tables(l, type='means')
pdf(file="HypIsPrimedRatioZoomin.pdf",width=12,height=6)
evalq(interaction.plot(hyp, IsPrimed, ratio), e)
plot(e$ratio,e$hyp)
pdf(file="HypIsPrimedRatio.pdf",width=12,height=6)
intxplot(ratio~hyp, groups=IsPrimed, data=e, se=TRUE, ylim=c(0.16,0.33), xlim=c(0.75,3.5) )
dev.off()
evalq(plot(hyp,ratio), e)
dev.off()
l0 <- lm(ratio~hyp, data=e[e$hyp!=3,])
l1 <- lm(ratio~1, data=e[e$hyp!=3,])
BIC(l0); BIC(l1)
BF=exp(-.5*(BIC(l1)-BIC(l0)));BF
## Again, there is no difference between highs and lows (Why? because
## the Bayes Factor is favoring l1, quite strongly. And l1 is the
## model that assumes no difference between groups.)

## OK, but people could ask us what about the effect of blocks?


e <- evalq(aggregate(list(notSeen=vis), list(suj=suj, hyp=hyp, block=block, IsPrimed=IsPrimed), function(x){sum(x==1)}),
           data[data$block %in% c(2, 4) & data$IsPrimed<2,]) # Only Block 2 (it's significant when I use mean instead)
f <- evalq(aggregate(list(notSeen=vis), list(suj=suj, hyp=hyp, block=block, IsPrimed=IsPrimed), function(x){sum(x<5)}),
           data[data$block %in% c(2, 4) & data$IsPrimed<2,]) # Only Block 2 (it's significant when I use mean instead)
e$ratio <- e$notSeen/f$notSeen
e$hyp<-factor(e$hyp)
e$block<-factor(e$block)
e$IsPrimed<-factor(e$IsPrimed)
l <- aov(ratio~hyp*block*IsPrimed+Error(suj/(block*IsPrimed)), data=e)
summary(l); model.tables(l, type='means')
xtable(l)
evalq(interaction.plot(hyp, block, ratio), e)
l0 <- lm(ratio~hyp*block, data=e)
l1 <- lm(ratio~hyp, data=e)
BIC(l0); BIC(l1)
BF=exp(-.5*(BIC(l1)-BIC(l0)));BF

## Donc il ne faut pas mettre pas l'effet de block.


l <- glmer(notSeen~hyp*IsPrimed+(1|suj), family="poisson", data=e)
summary(l)
ee <- aggregate(e$notSeen, list(hyp=e$hyp), mean)

l0 <- glmer(notSeen~hyp+(1|suj), family="poisson", data=e[e$hyp!=3,])
l1 <- glmer(notSeen~(1|suj), family="poisson", data=e[e$hyp!=3,])
summary(l1)
BIC(l0); BIC(l1)
BF=exp(.5*(BIC(l1)-BIC(l0)));BF

## It seems that we can conclude that there is no effect of group...


e <- evalq(aggregate(list(rt=rt, acc=acc, vis=vis), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed), median),
           data[data$BlockNum=="2",]) # Only Block 2 (it's significant when I use mean instead)

e <- evalq(aggregate(list(rt=rt, acc=acc, vis=vis), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed), median),
           data[data$BlockNum=="4",]) # Only Block 4 (significant both mean and median)

e <- evalq(aggregate(list(rt=rt, acc=acc, vis=vis), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed), mean),
           data[( (data$BlockNum=="2") | (data$BlockNum=="4")) & data$IsPrimed<2,]) # Both visibility blocks (significant both)

e$hyp<-factor(e$hyp); e$IsPrimed <- factor(e$IsPrimed)

l <- aov(vis~hyp, data=e[e$IsPrimed<2,]) #Because I don't test trials were there isn't a prime
l <- aov(vis~hyp*IsPrimed+Error(suj/(IsPrimed)), data=e) #Because I don't test trials were there isn't a prime

summary(l)
model.tables(l, type='means')
evalq(interaction.plot(hyp, IsPrimed, vis), e)







TukeyHSD(l) # Only the diffs between High and No hypno are significant after tukey test
evalq(plot(hyp,vis),e[e$IsPrimed<2,]) #turn to means before plot for nicer results

##Do Hypnosis and Prime Duration interact for Visibility ratings? YES

e <- evalq(aggregate(list(rt=rt, acc=acc, vis=vis), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed, PrimeDur=PrimeDur), median),
           data[data$BlockNum=="2",]) # Only Block 2 (PrimeDur Main effect and that's it)

e <- evalq(aggregate(list(rt=rt, acc=acc, vis=vis), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed, PrimeDur=PrimeDur), median),
           data[data$BlockNum=="4",]) # Only Block 4 (significant both mean and median)

e <- evalq(aggregate(list(rt=rt, acc=acc, vis=vis), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed, PrimeDur=PrimeDur), median),
           data[(data$BlockNum=="2") | (data$BlockNum=="4"),]) # Both visibility blocks (NOT SIGNIFICANT)

e <- evalq(aggregate(list(vis=vis), list(suj=suj, hyp=hyp, PrimeDur=PrimeDur, block=block, IsPrimed=IsPrimed), mean),
           data[(data$BlockNum=="2") | (data$BlockNum=="4"),]) # Both visibility blocks (NOT SIGNIFICANT)

e$hyp<-factor(e$hyp)
e$PrimeDur<-factor(e$PrimeDur)
e$block<-factor(e$block)

l <- aov(vis~hyp*block*PrimeDur+Error(suj/(PrimeDur*block)), data=e[e$IsPrimed<3,])
summary(l)
xtable(l)
l0 <- lm(vis~hyp*PrimeDur*block, data=e[e$IsPrimed<3,])
l1 <- lm(vis~PrimeDur, data=e[e$IsPrimed<3,])
BIC(l0); BIC(l1)
BF=exp(-.5*(BIC(l1)-BIC(l0)));BF

xtable(l)
model.tables(l, type='means')

TukeyHSD(l, "hyp") # Why wont tukey work??
evalq(interaction.plot(PrimeDur, hyp, vis),e[e$IsPrimed<3,]) #for just block 4
pdf(file="HypXPrimeDur.pdf",width=12,height=6)
intxplot(vis~PrimeDur, groups=hyp, data=e, se=TRUE, ylim=range(e$vis))
dev.off()

###########
e$bloc<-factor(e$bloc)
e$hyp<-factor(e$hyp)
e$PrimeDur<-factor(e$PrimeDur)

e$sd<-1
Allsd<-by(e$vis,e$hyp,sd)
e[e$hyp=="1",]$sd<-Allsd[1]
e[e$hyp=="2",]$sd<-Allsd[2]
e[e$hyp=="3",]$sd<-Allsd[3]
Allsd2<-by(e[e$hyp=="1",]$vis,e[e$hyp=="1",]$PrimeDur,sd)
Allsd2<-by(e[e$hyp=="2",]$vis,e[e$hyp=="2",]$PrimeDur,sd)
Allsd2<-by(e[e$hyp=="3",]$vis,e[e$hyp=="3",]$PrimeDur,sd)


opar <- theme_update(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(colour = "black"))
gp <- ggplot(e[e$hyp=="3",], aes(PrimeDur, vis, colour=hyp, group=hyp))
gp + geom_line(aes(linetype=hyp), size=2) + 
  geom_point(aes(shape=hyp), size=4) + 
  geom_errorbar(aes(ymax=vis+sd, ymin=vis-(sd), width=.1)
theme_set(opar)

##########

##Do Hypnosis and Prime Duration interact for RT for the visibility task?
##Only a Hypnosis main effect, no prime~hyp interaction

e <- evalq(aggregate(list(rtVis=rtVis, acc=acc, vis=vis), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed, PrimeDur=PrimeDur), median),
           data[data$BlockNum=="2",]) # Only Block 2 (hyp main effect, but something's not right maybe)

e <- evalq(aggregate(list(rtVis=rtVis, acc=acc, vis=vis), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed, PrimeDur=PrimeDur), median),
           data[data$BlockNum=="4",]) # Only Block 4 (NOT SIGNIFICANT)

e <- evalq(aggregate(list(rtVis=rtVis, acc=acc, vis=vis), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed, PrimeDur=PrimeDur), median),
           data[(data$BlockNum=="2") | (data$BlockNum=="4"),]) # Both visibility blocks (hyp main effect, but something's not right maybe)

e$hyp<-factor(e$hyp)
e$PrimeDur<-factor(e$PrimeDur)

l <- aov(rtVis~hyp*PrimeDur+Error(suj/(PrimeDur)), data=e[e$IsPrimed<3,]) 
l <- aov((rtVis~hyp), data=e[e$IsPrimed<3,]) 
summary(l)
model.tables(l, type='means')

TukeyHSD(l, "hyp") # 1-3 2-3 VERY SIGNIFICANT. DIFF BETWEEN 1 and 2 NOT AT ALL
evalq(interaction.plot(PrimeDur, hyp, rtVis),e[e$IsPrimed<3,]) #
evalq(plot(hyp,rtVis, ylim=c(700,1000)),e[e$IsPrimed<3,]) #


##Does Hypnosis affect Accuracy and RTs on the main task? NOT AT ALL

e <- evalq(aggregate(list(rt=rt, acc=acc), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed), mean),
           data[data$BlockNum=="1",]) # Only Block 1
e <- evalq(aggregate(list(rt=rt, acc=acc), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed), mean),
           data[data$BlockNum=="4",]) # Only Block 4 
e <- evalq(aggregate(list(rt=rt, acc=acc, vis=vis), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed), mean),
           data[(data$BlockNum=="1") | (data$BlockNum=="4"),]) # Both

e$hyp<-factor(e$hyp)

l <- aov(acc~hyp, data=e[e$IsPrimed<3,]) # NO HYPNOSIS~ACC INTERACTION FOR THE MAIN TASK
l <- aov(rt~hyp, data=e[e$IsPrimed<3,]) # NO HYPNOSIS~RT INTERACTION FOR THE MAIN TASK

summary(l)

data$cong <- data$IsPrimed
e <- evalq(aggregate(list(rt=rt, acc=acc), list(suj=suj, hyp=hyp, cong=cong, dur=dur, block=block), median),
           droplevels(data[data$block %in%c(1,4) & data$cong!=2,]))
e$hyp<-factor(e$hyp); e$bloc <- factor(e$block); e$cong <- factor(e$cong); e$dur <- factor(e$dur)
l <- aov(acc~hyp*cong*dur*block+Error(suj/(cong*dur*block)), data=e)
summary(l)
evalq(interaction.plot(hyp, cong, acc), e)

l <- aov(rt~hyp*cong*dur*block+Error(suj/(cong*dur*block)), data=e)
summary(l)
evalq(interaction.plot(dur, block, rt, main='incongruent'), e[e$cong==0,])
evalq(interaction.plot(dur, block, rt, main='congruent'), e[e$cong==1,])

data2<-data[data$dur>0,]

data$log<-log(data$acc)
data$block <- factor(data$block)
data$hyp <- factor(data$hyp)
data$dur <- factor(data$dur)
data$cong <- factor(data$cong)
data$Quadrant <- factor(data$Quadrant)

contrasts(data$Quadrant)=contr.sum(4)
contrasts(data$block)=contr.sum(4)
contrasts(data$hyp)=contr.sum(3)
contrasts(data$dur)=contr.sum(5)
contrasts(data$cong)=contr.sum(2)

data2$log<-log(data2$acc)
data2$block <- factor(data2$block)
data2$hyp <- factor(data2$hyp)
data2$dur <- factor(data2$dur)
data2$cong <- factor(data2$cong)
data2$Quadrant <- factor(data2$Quadrant)

contrasts(data2$Quadrant)=contr.sum(4)
contrasts(data2$block)=contr.sum(4)
contrasts(data2$hyp)=contr.sum(3)
contrasts(data2$dur)=contr.treatment(4)
contrasts(data2$cong)=contr.sum(2)

l <- glmer(acc~hyp*cong+(1|suj), family='binomial',
           data=droplevels(data[data$block %in%c(1,4) & data$cong!=2,]))

l <- glmer(acc~hyp*dur+Quadrant+(1|suj), family='binomial',
           data=droplevels(data[data$block %in%c(3) & data$hyp!=2,]))

l <- glmer(acc~hyp*dur+(1|suj), family='binomial',
           data=droplevels(data2[data2$block %in%c(3),]))
summary(glht(l, mcp(hyp ="Tukey")))

hyps <- c(1,2,3)
par(mfrow=c(2,2))
lapply(hyps, function(h){
           plot(density(data$rt[data$hyp==h & data$cong==0 & data$block%in%c(1,4)]), col=2)
           lines(density(data$rt[data$hyp==h & data$cong==1 & data$block%in%c(1,4)]), col=1)
})

evalq(plot(density(rt), col=3), data[data$hyp==3 & data$block%in%c(1,4),])
evalq(lines(density(rt), col=1), data[data$hyp==1 & data$block%in%c(1,4),])
evalq(lines(density(rt), col=2), data[data$hyp==2 & data$block%in%c(1,4),])


##Does Hypnosis affect Accuracy and RTs for the task on the prime? NOT AT ALL

e <- evalq(aggregate(list(rt=rt, acc=acc), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed), mean),
           data[data$BlockNum=="3",]) # Only Block 3 has the task on the prime

e$hyp<-factor(e$hyp)

l <- aov(acc~hyp, data=e[e$IsPrimed<2,]) # NO HYPNOSIS~ACC INTERACTION FOR THE MAIN TASK
l <- aov(rt~hyp, data=e[e$IsPrimed<2,]) # NO HYPNOSIS~RT INTERACTION FOR THE MAIN TASK

summary(l) 


##Do we have Priming? NOT AT ALL IT WOULD SEEM (IN PROGRESS)

e <- evalq(aggregate(list(rt=rt, acc=acc, vis=vis), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed, PrimeDur=PrimeDur>.04), mean),
           data[(data$BlockNum==4) & (data$IsPrimed<2),]) # Only Block 1 

e <- evalq(aggregate(list(rt=rt, acc=acc, vis=vis), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed, PrimeDur=PrimeDur), mean),
           data[(data$BlockNum=="4") & (data$IsPrimed<2),]) # Only Block 4 

e <- evalq(aggregate(list(rt=rt, acc=acc, vis=vis), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed, PrimeDur=PrimeDur, BlockNum=BlockNum), mean),
           data[((data$BlockNum=="1") | (data$BlockNum=="4")) & (data$IsPrimed<2),]) # Both 

e$hyp<-factor(e$hyp)
e$PrimeDur<-factor(e$PrimeDur)
e$IsPrimed<-factor(e$IsPrimed)

l <- aov((rt~IsPrimed), data=e) 
l <- aov(rt~hyp*IsPrimed*PrimeDur*BlockNum+Error(suj/(PrimeDur*IsPrimed*BlockNum)), data=e) 

summary(l)
model.tables(l, type='means')
evalq(interaction.plot(IsPrimed, PrimeDur, rt), e)



TukeyHSD(l, "hyp") # 1-3 2-3 VERY SIGNIFICANT. DIFF BETWEEN 1 and 2 NOT AT ALL
evalq(interaction.plot(PrimeDur, hyp, rtVis),e[e$IsPrimed<3,]) #
evalq(plot(hyp,rtVis, ylim=c(700,1000)),e[e$IsPrimed<3,]) #

##d-prime and criterion (blocks 2 and 4)
data$dur <- factor(round(1000*data$PrimeDur), ordered=TRUE)
data$seen <- as.numeric(data$vis>1)


df <- evalq(aggregate(list(hit=seen), list(suj=suj, hyp=hyp), mean),
            data[data$block %in% c(2, 4) & data$IsPrimed<2,])

dffa <- evalq(aggregate(list(fa=seen), list(suj=suj, hyp=hyp), mean),
            data[data$block %in% c(2, 4) & data$IsPrimed==2,])
df <- merge(df, dffa)


sum(df$fa==0); sum(df$fa==1); 
any(df$fa==0); any(df$fa==1); 
df$fa[df$fa==0] <- .01
sum(df$hit==0); sum(df$hit==1); 
any(df$hit==0); any(df$hit==1); 
df$hit[df$hit==0] <- .01
df$hit[df$hit==1] <- .99
df$dprime <- dprime(df$hit, df$fa)
df$beta <- beta(df$hit, df$fa)
evalq(interaction.plot(dur, hyp, dprime), df)
evalq(interaction.plot(dur, hyp, beta), df)

aggregate(df$dprime, list(hyp=df$hyp), mean)
            
f <- evalq(aggregate(list(notSeen=vis), list(suj=suj, hyp=hyp, PrimeDur=PrimeDur), function(x){sum(x<5)}),
           data[data$block %in% c(2, 4) & data$IsPrimed<2,]) # Only Block 2 (it's significant when I use mean instead)



library("grt")
library("psyphy")
dataDprime <-data[(data$BlockNum=="2") | (data$BlockNum=="4"),] # Both

#Create Seen (1) and Unseen (0)
dataDprime$Seen<-0;
for (i in (1:nrow(dataDprime))){
  if ((dataDprime$vis[i]>1) & (dataDprime$IsPrimed[i]<2)){dataDprime$Seen[i]<-1}
}

# Create False Alarm rate (1)

dataDprime$FalseAlarms<-0
for (i in (1:nrow(dataDprime))){
  if ((dataDprime$vis[i]>1) & (dataDprime$PrimeNoPrime[i]==0)){dataDprime$FalseAlarms[i]<-1}
}

#calculate dprime
dprimeH<-dataDprime[dataDprime$hyp==1,]
dprimeL<-dataDprime[dataDprime$hyp==2,]
dprimeN<-dataDprime[dataDprime$hyp==3,]  

dprimeH16<-dprimeH[dprimeH$PrimeDur==0.0167,]
dprimeH33<-dprimeH[dprimeH$PrimeDur==0.0334,]
dprimeH67<-dprimeH[dprimeH$PrimeDur==0.0668,]
dprimeH84<-dprimeH[dprimeH$PrimeDur==0.0835,]

dprimeL16<-dprimeL[dprimeL$PrimeDur==0.0167,]
dprimeL33<-dprimeL[dprimeL$PrimeDur==0.0334,]
dprimeL67<-dprimeL[dprimeL$PrimeDur==0.0668,]
dprimeL84<-dprimeL[dprimeL$PrimeDur==0.0835,]

dprimeN16<-dprimeN[dprimeN$PrimeDur==0.0167,]
dprimeN33<-dprimeN[dprimeN$PrimeDur==0.0334,]
dprimeN67<-dprimeN[dprimeN$PrimeDur==0.0668,]
dprimeN84<-dprimeN[dprimeN$PrimeDur==0.0835,]
    
  
dpH16<-(sum(dprimeH16$Seen))/nrow(dprimeH16)
dpH33<-(sum(dprimeH33$Seen))/nrow(dprimeH33)
dpH67<-(sum(dprimeH67$Seen))/nrow(dprimeH67)
dpH84<-(sum(dprimeH84$Seen))/nrow(dprimeH84)

dpL16<-(sum(dprimeL16$Seen))/nrow(dprimeL16)
dpL33<-(sum(dprimeL33$Seen))/nrow(dprimeL33)
dpL67<-(sum(dprimeL67$Seen))/nrow(dprimeL67)
dpL84<-(sum(dprimeL84$Seen))/nrow(dprimeL84)

dpN16<-(sum(dprimeN16$Seen))/nrow(dprimeN16)
dpN33<-(sum(dprimeN33$Seen))/nrow(dprimeN33)
dpN67<-(sum(dprimeN67$Seen))/nrow(dprimeN67)
dpN84<-(sum(dprimeN84$Seen))/nrow(dprimeN84)


                              

        
##Does Hypnosis affect/improves/eliminates Priming?

##Does Hypnosis interact with the single/double task binome?




e <- evalq(aggregate(list(rt=rt), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed), median),
                                                                      data[data$BlockNum=="1",])
e <- evalq(aggregate(list(rt=rt), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed), median),
           data[data$BlockNum=="2",])
e <- evalq(aggregate(list(rt=rt, acc=acc, vis=vis), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed), mean),
           data[data$BlockNum=="2",])
e <- evalq(aggregate(list(rt=rt), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed, vis=vis), median),
           data[data$BlockNum=="4",])
e$vis<-factor(e$vis)
e$hyp<-factor(e$hyp)
e$vis<-as.numeric(as.character(e$vis))

e <- evalq(aggregate(list(rt=rt), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed), median),
           data[data$BlockNum=="4",])

e <- evalq(aggregate(list(rt=rt), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed), median),
           data[(data$BlockNum=="1") | (data$BlockNum=="4"),])

l <- aov(vis~hyp, data=e[e$IsPrimed<2,])
l <- aov(rt~hyp*IsPrimed+Error(suj/(IsPrimed)), data=e)
l <- aov(rt~hyp*vis+Error(suj/(vis)), data=e)
summary(l)
model.tables(l, type='means')
evalq(interaction.plot(hyp,vis), e)
evalq(plot(hyp,vis),e[e$IsPrimed<2,])


##Testing accuracy and response time on block 3 (main task on prime)

df <- droplevels(data[data$block==3 & data$IsPrimed<3 & data$hyp<3 ,])
df <- droplevels(data[data$block==3 & data$IsPrimed<3,])

e <- evalq(aggregate(list(rt=rt, acc=acc), list(suj=suj, hyp=hyp, dur=dur), mean),
           df[df$IsPrimed<2 & df$dur!=17,])
e <- evalq(aggregate(list(rt=rt, acc=acc), list(suj=suj, hyp=hyp, PrimeDur=PrimeDur), mean),
           df[df$IsPrimed<2,])
e$hyp <- factor(e$hyp)
e$dur<-factor(e$dur)

l <- aov(acc~hyp*dur+Error(suj/(dur)), data=e)
l <- aov(rt~hyp*dur+Error(suj/(dur)), data=e)
summary(l); model.tables(l, type='means')
xtable(l)

df$hyp <- factor(df$hyp)
df$dur<-as.numeric(as.character(factor(df$dur)))

contrasts(df$hyp) <- contr.sum(3)
contrasts(df$dur) <- contr.sum(5)

lm<-glmer(acc~hyp*PrimeDur+(1|suj), data=df, family='binomial')

summary(lm)

pdf(file="HypObjectiveVis.pdf",width=12,height=6)
intxplot(acc~dur, groups=hyp, data=df, se=TRUE, ylim=c(0.4,1), xlim=c(0.75,5.5)  )
dev.off()

require(nlme)
# Below call does not work for glht, thus we created the interaction term in the data frame
#model <- lme(response ~ group*time, random = ~ 1 | subject / time, dat)
model <- lme(acc ~ hyp*PrimeDur, random = ~ 1 | suj / PrimeDur, e)
require(multcomp)
summary(glht(model, linfct=mcp(grtim="Tukey")), test = adjusted(type = "bonferroni"))

evalq(interaction.plot(dur, hyp, acc), e)
evalq(interaction.plot(PrimeDur, hyp, acc), e)


library(lme4)

df$hyp <- factor(df$hyp)
contrasts(df$hyp) <- contr.sum(3)
contrasts(df$dur) <- contr.sum(4)

l <- glmer(acc~hyp*dur+(1|suj), data=df, family='binomial')
summary(l)

evalq(interaction.plot(dur, hyp, acc), e)

##Test for priming
## Definition of repetitions

data$rep <- 1
data$rep[data$Target!=data$Prime] <- 0

data$perf <- data$rt
data$perf[data$block==2] <- data$rtVis[data$block==2]

e <- evalq(aggregate(list(perf=perf, acc=acc), list(suj=suj, block=block, hyp=hyp, IsPrimed=IsPrimed,dur=dur), mean), data[data$IsPrimed<2 & data$block!=2,])

e$asinacc <- asin(e$acc)
e$IsPrimed<-factor(e$IsPrimed)
e$block<-factor(e$block)
e$dur<-factor(e$dur)
e$hyp<-factor(e$hyp)

l <- aov(perf~hyp*IsPrimed*block*dur+Error(suj/(IsPrimed*block*dur)), data=e[(e$block==1 | e$block==3 | e$block==4) & e$hyp!=3 ,])
l <- aov(acc~hyp*IsPrimed*dur+Error(suj/(IsPrimed*dur)), data=e[(e$block==1 | e$block==3 | e$block==4) & e$hyp!=3 ,])
l <- aov(acc~hyp*IsPrimed*dur+Error(suj/(IsPrimed*dur)), data=e[e$block==3,])
l <- aov(asinacc~hyp*IsPrimed*block+Error(suj/(IsPrimed*block)), data=e)
l <- aov(asinacc~hyp*IsPrimed+Error(suj/(IsPrimed)), data=e)
summary(l)
par(mfrow=c(2,2))
evalq(interaction.plot(dur, IsPrimed, acc), e[e$block==1 | e$block==3 | e$block==4  ,])
evalq(interaction.plot(block, IsPrimed, acc), e)
evalq(interaction.plot(hyp, IsPrimed, acc), e)
evalq(interaction.plot(block, dur, acc), e)
evalq(interaction.plot(dur, IsPrimed, acc), e[e$block==3,])
evalq(interaction.plot(hyp, IsPrimed, perf), e[e$block==4 | e$block==1, ])
l <- aov(perf~hyp*IsPrimed*block+Error(suj/(IsPrimed*block)), data=e[e$dur!=17,])
summary(l)
model.tables(l, type='means')


data$hyp2<-data$hyp
data$IsPrimed2<-data$IsPrimed
data$block2<-data$block


data$hyp2<-factor(data$hyp2)
data$IsPrimed2<-factor(data$IsPrimed2)
data$block2<-factor(data$block2)


library(lme4)
l1 <- glmer(acc~hyp2*IsPrimed2*block2+(1|suj), family='binomial', data=data[data$IsPrimed<2 & data$block!=2,])

l1 <- glmer(acc~hyp2*IsPrimed2*block2+(1|suj), family='binomial', data=data[data$IsPrimed<2 & data$block!=2,])
l0 <- glmer(acc~IsPrimed2*block2+(1|suj), family='binomial', data=data[data$IsPrimed<2 & data$block!=2,])
BIC(l1); BIC(l0)
BF=exp(-.5*(BIC(l1)-BIC(l0)));BF
summary(l0)


data$subli <- 1
data$subli[data$dur>17] <- 0
data$subli <- factor(data$subli)

e <- evalq(aggregate(list(perf=perf, acc=acc), list(suj=suj, block=block, hyp=hyp, IsPrimed=IsPrimed, subli=subli), mean), data[data$IsPrimed<2 & data$block!=2,])

e$asinacc <- asin(e$acc)

l <- aov(acc~hyp*IsPrimed*subli+Error(suj/(IsPrimed*subli)), data=e[e$block=="3" | e$block=="3",])
l <- aov(perf~hyp*IsPrimed*block*subli+Error(suj/(IsPrimed*block*subli)), data=e)

summary(l)

l1 <- glmer(acc~IsPrimed*block+(1|suj), family='binomial', data=data[data$subli==0 & data$IsPrimed<2 & data$block!=2,])
l0 <- glmer(acc~block+(1|suj), family='binomial', data=data[data$subli==0 & data$IsPrimed<2 & data$block!=2,])
BIC(l1); BIC(l0)
BF=exp(.5*(BIC(l1)-BIC(l0)));BF
summary(l0)

evalq(interaction.plot(subli, IsPrimed, acc), e[e$block=="3" | e$block=="3",])
evalq(interaction.plot(block, IsPrimed, perf), e)
evalq(interaction.plot(block, dur, perf), e)
evalq(interaction.plot(dur, IsPrimed, perf), e[e$block==4,])
evalq(interaction.plot(hyp, IsPrimed, perf), e[e$block==4 | e$block==1, ])


#####
##Calculate Correlations

Highs<-evalq(aggregate(list(Correctness=Correctness, Visibility=Visibility), 
                       list(suj=suj, hyp=hyp, bloc=BlockNum, Dur=PrimeDur), mean), data[data$BlockNum %in% c(2,3) & data$hyp=="1",])
Lows<-evalq(aggregate(list(Correctness=Correctness, Visibility=Visibility), 
                      list(suj=suj, hyp=hyp, bloc=BlockNum, Dur=PrimeDur), mean), data[data$BlockNum %in% c(2,3) & data$hyp=="2",])
NoHyp<-evalq(aggregate(list(Correctness=Correctness, Visibility=Visibility), 
                       list(suj=suj, hyp=hyp, bloc=BlockNum, Dur=PrimeDur), mean), data[data$BlockNum %in% c(2,3) & data$hyp=="3",])

HighsCorr<-as.data.frame(cbind(Highs$Visibility[Highs$bloc=="2"], Highs$Correctness[Highs$bloc=="3"],
                               Highs$hyp[Highs$bloc=="2"], Highs$suj[Highs$bloc=="2"], Highs$Dur[Highs$bloc=="2"]))
LowsCorr<-as.data.frame(cbind(Lows$Visibility[Lows$bloc=="2"], Lows$Correctness[Lows$bloc=="3"],
                              Lows$hyp[Lows$bloc=="2"], Lows$suj[Lows$bloc=="2"],Lows$Dur[Lows$bloc=="2"]))
NoHypCorr<-as.data.frame(cbind(NoHyp$Visibility[NoHyp$bloc=="2"], NoHyp$Correctness[NoHyp$bloc=="3"],
                               NoHyp$hyp[NoHyp$bloc=="2"], NoHyp$suj[NoHyp$bloc=="2"], NoHyp$Dur[NoHyp$bloc=="2"]))

AllVis<-as.data.frame(rbind(HighsCorr,LowsCorr,NoHypCorr))
library(plyr)
library(Hmisc)
rename(AllVis, c("V1"="subjective", "V2"="objective", "V3"="hyp", "V4"="suj"))

RLOW<-rcorr(as.matrix(LowsCorr), type="pearson")
RHIGH<-rcorr(as.matrix(HighsCorr), type="pearson")
RNOH<-rcorr(as.matrix(NoHypCorr), type="pearson")
#http://www.danielsoper.com/statcalc/calculator.aspx?id=104
RALL<-rcorr(as.matrix(AllVis), type="pearson")

library(car) 
scatterplot(V1 ~ V2 | V3, data=AllVis,
            xlab="Objective vis", ylab="Subjective Vis", 
            main="Scatter Plot all hypno", 
            labels=row.names(AllVis))
scatterplot(V1 ~ V2 | V3, data=rbind(HighsCorr,LowsCorr),
                        xlab="Objective vis", ylab="Subjective Vis", 
                        main="Scatter Plot High + Low", 
                        ) 
scatterplot(V1 ~ V2 | V3, data=LowsCorr,
            xlab="Objective vis", ylab="Subjective Vis", 
            main="Scatter Plot Low", 
            ) 
scatterplot(V1 ~ V2 | V3, data=HighsCorr,
            xlab="Objective vis", ylab="Subjective Vis", 
            main="Scatter Plot High", 
            ) 
scatterplot(V1 ~ V2 | V3, data=NoHypCorr,
            xlab="Objective vis", ylab="Subjective Vis", 
            main="Scatter Plot No H", 
            ) 

plot(NoHypCorr$V1,NoHypCorr$V2)
plot(NoHypCorr$V1,NoHypCorr$V2)
plot(NoHypCorr$V1,NoHypCorr$V2)

splom(AllVis$V1,AllVis$V2)
pairs(AllVis$V1,AllVis$V2)

library("lme4")

AllVis$V2<-factor(AllVis$V2)
AllVis$V1<-factor(AllVis$V1)
AllVis$V3<-factor(AllVis$V3)
AllVis$V4<-factor(AllVis$V4)
AllVis$V5<-factor(AllVis$V5)


AllVis$V1<-as.numeric(as.character((AllVis$V1)))
AllVis$V2<-as.numeric(as.character((AllVis$V2)))
AllVis$V3<-as.numeric(as.character((AllVis$V3)))
AllVis$V4<-as.numeric(as.character((AllVis$V4)))
AllVis$V5<-as.numeric(as.character((AllVis$V5)))

l <- lmer(V1~V2+(1|V4), AllVis)
summary(l)
print(l, correlation=TRUE)