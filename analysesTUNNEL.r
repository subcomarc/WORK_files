setwd("C:/Users/lscpuser/Google Drive/R/")
data <- read.table(file="C:/Users/lscpuser/Google Drive/R/DataTable RAW.txt", header=T)
source("C:/Users/lscpuser/Google Drive/R/R_Lib_Sackur.r")
data$suj <- factor(data$SubjName)
data$hyp <- floor(data$SubjName/1000)
data$rt <- round(numerize(data$RTime)*1000)
data$acc <- data$Correctness
data$vis <- data$Visibility
data$block <- data$BlockNum
data$dur <- round(data$PrimeDur * 1000)
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
data <- droplevels(data[data$rt<3000 | data$block==2,])

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

e <- evalq(aggregate(list(rt=rt, acc=acc, vis=vis), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed), median),
           data[data$BlockNum=="2",]) # Only Block 2 (it's significant when I use mean instead)

e <- evalq(aggregate(list(rt=rt, acc=acc, vis=vis), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed), median),
           data[data$BlockNum=="4",]) # Only Block 4 (significant both mean and median)

e <- evalq(aggregate(list(rt=rt, acc=acc, vis=vis), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed), median),
           data[(data$BlockNum=="2") | (data$BlockNum=="4"),]) # Both visibility blocks (significant both)

e$hyp<-factor(e$hyp)

l <- aov(vis~hyp, data=e[e$IsPrimed<2,]) #Because I don't test trials were there isn't a prime

summary(l)
model.tables(l, type='means')

TukeyHSD(l) # Only the diffs between High and No hypno are significant after tukey test
evalq(plot(hyp,vis),e[e$IsPrimed<2,]) #turn to means before plot for nicer results

##Do Hypnosis and Prime Duration interact for Visibility ratings? YES

e <- evalq(aggregate(list(rt=rt, acc=acc, vis=vis), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed, PrimeDur=PrimeDur), median),
           data[data$BlockNum=="2",]) # Only Block 2 (PrimeDur Main effect and that's it)

e <- evalq(aggregate(list(rt=rt, acc=acc, vis=vis), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed, PrimeDur=PrimeDur), median),
           data[data$BlockNum=="4",]) # Only Block 4 (significant both mean and median)

e <- evalq(aggregate(list(rt=rt, acc=acc, vis=vis), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed, PrimeDur=PrimeDur), median),
           data[(data$BlockNum=="2") | (data$BlockNum=="4"),]) # Both visibility blocks (NOT SIGNIFICANT)

e$hyp<-factor(e$hyp)
e$PrimeDur<-factor(e$PrimeDur)

l <- aov(vis~hyp*PrimeDur+Error(suj/(PrimeDur)), data=e[e$IsPrimed<3,]) 
summary(l)
model.tables(l, type='means')

TukeyHSD(l, "hyp") # Why wont tukey work??
evalq(interaction.plot(PrimeDur, hyp, vis),e[e$IsPrimed<3,]) #for just block 4

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

##Does Hypnosis affect Accuracy and RTs for the task on the prime? NOT AT ALL

e <- evalq(aggregate(list(rt=rt, acc=acc), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed), mean),
           data[data$BlockNum=="3",]) # Only Block 3 has the task on the prime

e$hyp<-factor(e$hyp)

l <- aov(acc~hyp, data=e[e$IsPrimed<2,]) # NO HYPNOSIS~ACC INTERACTION FOR THE MAIN TASK
l <- aov(rt~hyp, data=e[e$IsPrimed<2,]) # NO HYPNOSIS~RT INTERACTION FOR THE MAIN TASK

summary(l) 


##Do we have Priming? NOT AT ALL IT WOULD SEEM (IN PROGRESS)

e <- evalq(aggregate(list(rt=rt, acc=acc, vis=vis), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed, PrimeDur=PrimeDur), mean),
           data[(data$BlockNum=="1") & (data$IsPrimed<2),]) # Only Block 1 

e <- evalq(aggregate(list(rt=rt, acc=acc, vis=vis), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed, PrimeDur=PrimeDur), mean),
           data[(data$BlockNum=="4") & (data$IsPrimed<2),]) # Only Block 4 

e <- evalq(aggregate(list(rt=rt, acc=acc, vis=vis), list(suj=suj, hyp=hyp, IsPrimed=IsPrimed, PrimeDur=PrimeDur), mean),
           data[((data$BlockNum=="1") | (data$BlockNum=="4")) & (data$IsPrimed<2),]) # Both 

e$hyp<-factor(e$hyp)
e$PrimeDur<-factor(e$PrimeDur)
e$IsPrimed<-factor(e$IsPrimed)

l <- aov((rt~IsPrimed), data=e) 
l <- aov(rt~hyp*IsPrimed*PrimeDur+Error(suj/(PrimeDur*IsPrimed)), data=e) 

summary(l)
model.tables(l, type='means')

TukeyHSD(l, "hyp") # 1-3 2-3 VERY SIGNIFICANT. DIFF BETWEEN 1 and 2 NOT AT ALL
evalq(interaction.plot(PrimeDur, hyp, rtVis),e[e$IsPrimed<3,]) #
evalq(plot(hyp,rtVis, ylim=c(700,1000)),e[e$IsPrimed<3,]) #

##d-prime and criterion

library("grt")



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

