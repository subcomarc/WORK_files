data <- read.table(file="C:/Users/lscpuser/Google Drive/R/BlinkTable.txt", header=T)
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
hist(data$rtOr[data$rtOr < 2000])

data$rtn <- data$RTimeHM*1000
hist(data$rtn[data$rtn<20000])

## 3: high: 2: medium; 1: low
data$hypno <- data$HLRAWcol
data$hypnoC <- data$HLcol
data$lag <- data$AllLags
data$soa <- data$AllSOAs

data$lag <- ordered(data$lag)
data$target1 <- factor(data$target1)
data$target2 <- factor(data$target2)
data$suj <- factor(data$suj)

## outliers

e <- evalq(aggregate(list(accn=accn, accOr=accOr, hypnoC=hypnoC), list(suj=suj), mean), data)
evalq(plot(accn, accOr, col=hypnoC,type='p', xlim=c(0, 1), ylim=c(0,1)), e)
l <- lm(accOr~accn, e)
abline(l)

e$acc <- (e$accOr+e$accn)/2
hist(e$acc); mean(e$acc); sd(e$acc)

sum(e$acc < .6)
bad <- as.vector(e$suj[e$acc < .6])

e <- e[!e$suj %in% bad,]
min(e$acc); hist(e$acc)

data <- droplevels(data[!data$suj %in% bad,])

sum(data$rtOr > 4000)

hist(data$rtOr[data$rtOr<4000])

e <- evalq(aggregate(list(medianRtOr=rtOr, hypno=hypno, hypnoC=hypnoC), list(suj=suj), median),
           data[data$rtOr < 10000,])
evalq(plot(hypno, medianRtOr, col=hypnoC,type='p'), e)
boxplot(e$medianRtOr)
sd(e$medianRtOr)
mean(e$medianRtOr)+3*sd(e$medianRtOr)

bad <- as.vector(e$suj[e$medianRtOr > 1600])
data <- droplevels(data[!data$suj %in% bad,])

e <- evalq(aggregate(list(medianRtOr=rtOr, hypno=hypno, hypnoC=hypnoC), list(suj=suj), median),
           data[data$rtOr < 6000,])
e$sd <- evalq(aggregate(rtOr, list(suj=suj), sd), data)$x

data <- merge(data, e)

data$outlier <- data$rtOr > data$medianRtOr + 3 * data$sd

e <- evalq(aggregate(list(outlier=outlier), list(suj=suj), mean), data)
any(e$outlier>.1); max(e$outlier)
data <- droplevels(data[!data$outlier & data$rtOr>150,])
hist(data$rtOr[data$rtOr < 300])
data$hypnoC <- ordered(data$hypnoC)
hist(data$rtOr)

e <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj, lag=lag), mean),
           data[data$nface==2,])
f <- evalq(aggregate(list(accOr=accOr, accn=accn), list(lag=lag), mean), e)
eprime <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj), mean), data[data$nface==1,])
fprime <- mean(eprime$accOr)

e <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj, lag=lag, target1=target1, target2=target2), mean),
           data[data$nface==2,])
f <- evalq(aggregate(list(accOr=accOr, accn=accn), list(lag=lag, target1=target1, target2=target2), mean), e)

f$acc <- (f$accOr+f$accn)/2

par(mfrow=c(2,2))
evalq(interaction.plot(lag, target1, accOr, main="ACC OR T2 angry", ylim=c(.6, .90)), f[f$target2=="angry",])
evalq(interaction.plot(lag, target1, accOr, main="ACC OR T2 neutral", ylim=c(.6, .90)), f[f$target2=="neutral",])

evalq(interaction.plot(lag, target1, accn, main="ACC N T2 angry", ylim=c(.4, 1)), f[f$target2=="angry",])
evalq(interaction.plot(lag, target1, accn, main="ACC NT2 neutral", ylim=c(.4, 1)), f[f$target2=="neutral",])

e$lag <- ordered(e$lag)
e$target1 <- factor(e$target1)
e$target2 <- factor(e$target2)
e$suj <- factor(e$suj)

data$lagC <- "early"
data$lagC[data$lag>5] <- "late"
data$lagC <- factor(data$lagC)
e <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj, lagC=lagC, target1=target1, target2=target2), mean),
           data[data$nface==2,])
l <- aov(accn~lagC*target1*target2+Error(suj/(lagC*target1*target2)), data=e)
summary(l)

e <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj, lagC=lagC, target1=target1, target2=target2, hypnoC=hypnoC), mean),
           data[data$nface==2,])
l <- aov(accn~lagC*target1*target2*hypnoC+Error(suj/(lagC*target1*target2)), data=e)
summary(l)

e <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj, target1=target1, target2=target2, hypnoC=hypnoC), mean),
           data[data$nface==2,])
l <- aov(accn~target1*target2*hypnoC+Error(suj/(target1*target2)), data=e)
summary(l)




l <- aov(accOr~lag*target1*target2+Error(suj/(lag*target1*target2)), data=e)
summary(l)
model.tables(l, type='means')
f <- evalq(aggregate(list(accOr=accOr), list(target1=target1, target2=target2), mean), e)
evalq(interaction.plot(target1, target2, accOr), e)

e <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj, lag=lag, target1=target1, target2=target2, hypnoC=hypnoC), mean),
           data[data$nface==2,])

l <- aov(accOr~lag*target1*target2*hypnoC+Error(suj/(lag*target1*target2)), data=e)
summary(l)
model.tables(l, type='means')

e <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj, target1=target1, target2=target2, hypnoC=hypnoC), mean),
           data[data$nface==2,])
l <- aov(accn~target1*target2*hypnoC+Error(suj/(target1*target2)), data=e)
summary(l)
model.tables(l, type='means')
par(mfrow=c(2,2))
evalq(interaction.plot(target1, target2, acc0r, main="low", ylim=c(.67, .77)), e[e$hypnoC==1,])
evalq(interaction.plot(target1, target2, accOr, main="medium", ylim=c(.67, .77)), e[e$hypnoC==2,])
evalq(interaction.plot(target1, target2, accOr, main="high", ylim=c(.67, .77)), e[e$hypnoC==3,])
ef <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj, target1=target1, target2=target2), mean),
           data[data$nface==2,])
evalq(interaction.plot(target1, target2, accOr, main="global", ylim=c(.67, .77)), ef)

par(mfrow=c(2,2))
evalq(interaction.plot(target1, target2, accn, main="low", ylim=c(.7, .82)), e[e$hypnoC==1,])
evalq(interaction.plot(target1, target2, accn, main="medium", ylim=c(.7, .82)), e[e$hypnoC==2,])
evalq(interaction.plot(target1, target2, accn, main="high", ylim=c(.7, .82)), e[e$hypnoC==3,])
ef <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj, target1=target1, target2=target2), mean),
           data[data$nface==2,])
evalq(interaction.plot(target1, target2, accOr, main="global", ylim=c(.67, .77)), ef)


library(lme4)
contrasts(data$target1)<-contr.sum(2)
contrasts(data$target2)<-contr.sum(2)
contrasts(data$hypnoC)<-contr.sum(3)
contrasts(data$lagC)<-contr.sum(2)
l <- glmer(accn~lagC*target1*target2*hypnoC+(1|suj), family='binomial', data=data)
summary(l)


l <- glmer(accOr~target2*HLRAWcol+(1|suj), family='binomial', data=data)
summary(l)

l <- glmer(accOr~target2*HLRAWcol*lagC+(1|suj), family='binomial', data=data)
summary(l)

e <- evalq(aggregate((blinkIndex=accOr), list(suj=suj, harvard=HLRAWcol), mean), data[data$target2=='angry' & data$target1=="angry" & data$hypnoC==2,])
e$neutral <- evalq(aggregate((blinkIndex=accOr), list(suj=suj, harvard=HLRAWcol), mean), data[data$target2=='neutral' & data$target1=="angry" & data$hypnoC==2,])$x
e$blinkI <- e$x-e$neutral
plot(e$harvard, e$blinkI, type="p")
l <- lm(e$blinkI~e$harvard); abline(l)


e <- evalq(aggregate(list(rtOr=rtOr), list(suj=suj, lag=lag, target1=target1, target2=target2), mean),
           data[data$nface==2,])
f <- evalq(aggregate(list(rtOr=rtOr), list(lag=lag, target1=target1, target2=target2), mean), e)

par(mfrow=c(2,2))
evalq(interaction.plot(lag, target1, rtOr, main="T2 angry", ylim=c(750, 1100)), f[f$target2=="angry",])
evalq(interaction.plot(lag, target1, rtOr, main="T2 neutral", ylim=c(750, 1100)), f[f$target2=="neutral",])



e <- evalq(aggregate(list(rtOr=rtOr), list(suj=suj, lag=lag, target1=target1, target2=target2, hypnoC=hypnoC), median),
           data[data$nface==2,])

l <- aov(rtOr~lag*target1*target2*hypnoC+Error(suj/(lag*target1*target2)), data=e)
summary(l)
model.tables(l, type='means')



