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
data$hypnoC <- ordered(data$hypnoC)
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
           data[data$rtOr < 6000,])
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

hist(data$rtOr)

e <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj, lag=lag), mean),
           data[data$nface==2,])
f <- evalq(aggregate(list(accOr=accOr, accn=accn), list(lag=lag), mean), e)
eprime <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj), mean), data[data$nface==1,])
fprime <- mean(eprime$accOr)

e <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj, lag=lag, target1=target1, target2=target2), mean),
           data[data$nface==2,])
f <- evalq(aggregate(list(accOr=accOr, accn=accn), list(lag=lag, target1=target1, target2=target2), mean), e)

par(mfrow=c(2,2))
evalq(interaction.plot(lag, target1, accOr, main="T2 angry", ylim=c(.6, .90)), f[f$target2=="angry",])
evalq(interaction.plot(lag, target1, accOr, main="T2 neutral", ylim=c(.6, .90)), f[f$target2=="neutral",])

e$lag <- ordered(e$lag)
e$target1 <- factor(e$target1)
e$target2 <- factor(e$target2)
e$suj <- factor(e$suj)


l <- aov(accOr~lag*target1*target2+Error(suj/(lag*target1*target2)), data=e)
summary(l)
model.tables(l, type='means')

e <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj, lag=lag, target1=target1, target2=target2, hypnoC=hypnoC), mean),
           data[data$nface==2,])

l <- aov(accOr~lag*target1*target2*hypnoC+Error(suj/(lag*target1*target2)), data=e)
summary(l)
model.tables(l, type='means')



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



