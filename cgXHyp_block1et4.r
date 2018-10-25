load("TunnelData.r")

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
data$IsSeen<-1 #for treating visibility as a binary
data$IsSeen[data$Visibility<2]<-0 #for treating visibility as a binary

data2 <- data[data$dur>0,] #to have a dataset without the trials without prime

data$block <- factor(data$block)
data$hyp <- factor(data$hyp)
data$dur <- factor(data$dur)
data$Quadrant <- factor(data$Quadrant)
data$IsPrimed <- factor(data$IsPrimed)
data$rep<- factor(data$rep)
data2$block <- factor(data2$block)
data2$hyp <- factor(data2$hyp)
data2$dur <- factor(data2$dur)
data2$Quadrant <- factor(data2$Quadrant)
data2$IsPrimed <- factor(data2$IsPrimed)
data2$rep<- factor(data2$rep)

contrasts(data$Quadrant)=contr.sum(4)
contrasts(data$block)=contr.sum(4)
contrasts(data$hyp)=contr.sum(3)
contrasts(data$dur)=contr.sum(5)
contrasts(data$rep)=contr.sum(2)
contrasts(data$IsPrimed)=contr.sum(3)
contrasts(data2$Quadrant)=contr.sum(4)
contrasts(data2$block)=contr.sum(4)
contrasts(data2$hyp)=contr.sum(3)
contrasts(data2$dur)=contr.sum(4)
contrasts(data2$rep)=contr.sum(2)
contrasts(data2$IsPrimed)=contr.sum(2)

data$target <- data$Target; data$prime <- data$Prime

## Priming En fait, on a bien un effet de priming de congruence, sur
## l'accuracy (rien sur les RTs), mais qui intéragit avec l'hypnose
## (modèle l3 plus bas): les résistants sont d'une part moins bons, et
## particulièrement pour les essais incongruents... Donc c'est bon:
## d'un sens l'hypnose fait disparaître l'effet de priming.

## Ça marche en gros pour l'ies aussi!! donc c'est cool.


df <- data[(data$block==1 | data$block==4) & data$PrimeNoPrime==1,]
df <- data[(data$block==4) & data$PrimeNoPrime==1,]
df <- data[(data$block==1) & data$PrimeNoPrime==1,]
df$logRt <- log(df$rt)
df$invRt <- 1/(df$rt)
df$cg <- factor(df$IsPrimed)
df$block <- factor(df$block)
df$dur <- as.numeric(as.character(df$dur))

DF <- droplevels(df[df$cg!=2 & df$rt<1200,])

DF$dur <- factor(DF$dur)

e <- evalq(aggregate(list(err=1-acc), list(suj=suj, cg=IsPrimed, dur=dur, hyp=hyp,
                                       block=block),
                     mean), droplevels(df[df$rt<1500 & df$IsPrimed!=2,]))
g <- evalq(aggregate(list(acc=acc), list(suj=suj, cg=IsPrimed, dur=dur, hyp=hyp,
                                       block=block),
                     mean), droplevels(df[df$rt<1500 & df$IsPrimed!=2,]))


f <- aggregate(e$err, list(cg=e$cg, dur=e$dur, hyp=e$hyp), mean)


par(mfrow=c(2,2))
My <- max(f$x); my <- min(f$x)
evalq(interaction.plot(dur, cg, x, col=1:2, ylim=c(my, My), main='Highs'),
      f[f$hyp==1,])
evalq(interaction.plot(dur, cg, x, col=1:2, ylim=c(my, My), main='Lows'),
      f[f$hyp==2,])
evalq(interaction.plot(dur, cg, x, col=1:2, ylim=c(my, My), main='Med without'),
      f[f$hyp==3,])

evalq(interaction.plot(hyp, cg, x, col=1:2, ylim=c(my, My), main='Hyp x Cg (no dur)'), f)



e <- evalq(aggregate(list(err=1-acc), list(suj=suj, cg=IsPrimed, dur=dur, hyp=hyp,
                                       block=block),
                     mean), droplevels(df[df$rt<1500 & df$IsPrimed!=2,]))

l <- aov(err~cg*hyp+block+dur+Error(suj/(cg+dur+block)), data=e); summary(l)

e <- evalq(aggregate(list(err=1-acc), list(suj=suj, cg=IsPrimed, hyp=hyp),mean),
           droplevels(df[df$rt<1500 & df$IsPrimed!=2,]))

l <- aov(err~cg*hyp+Error(suj/(cg)), data=e); summary(l)



library(lme4)
## l0 <- glmer(acc~cg*dur+(1|suj), family='binomial', optimizer= 'bobyqa',
##             data=DF)
## l1 <- glmer(acc~cg+dur+(1|suj),family='binomial', optimizer= 'bobyqa', data=DF)
## l2 <- glmer(acc~dur+(1|suj),family='binomial', optimizer= 'bobyqa', data=DF)

## l3 <- glmer(acc~cg*hyp*block+(1|suj), family='binomial', optimizer= 'bobyqa',
##             data=DF)
## summary(l3)

l4 <- glmer(acc~cg*hyp+(1|suj), family='binomial', optimizer= 'bobyqa',
            data=DF)

Anova(l4)

summary(l4)

DF <- droplevels(df[df$hyp!=3 & df$cg!=2 & df$rt<1200,])
contrasts(DF$hyp) <- contr.sum(2)

l5 <- glmer(acc~cg*hyp+(1|suj), family='binomial', optimizer= 'bobyqa',
            data=DF)
summary(l5)

