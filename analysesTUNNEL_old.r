data <- read.table(file="DataTable RAW.txt", header=T)


data$suj <- factor(data$SubjName)
data$hyp <- factor(round(data$SubjName/1000))
data$rt <- round(numerize(data$RTime)*1000)
data$acc <- data$Correctness
data$vis <- data$Visibility
data$block <- data$BlockNum
data$dur <- round(data$PrimeDur * 1000)
data$rtVis <- round(numerize(data$RTimeVis)*1000)

# Erreur de codage
# data$Hypnotiz[data$suj==3051] <- 3
data <- droplevels(data[data$suj!=0,])

# Erreur de réponse: les sujets confondent les touches de visibilité
                                        # et de catégorisation*

nrow(data[data$vis==8 | data$acc==8,])/nrow(data[data$vis!=8,])
data <- droplevels(data[data$vis!=8 & data$acc!=8,])
sujs <- unique(data$suj)

hist(data$rt[data$block!=2 & data$rt<6000])

# sujets outliers

e <- evalq(aggregate(list(acc=acc), list(suj=suj), mean), data[data$block==1,])
evalq(boxplot(rt~suj, ylim=c(0, 10000)), data[data$block==3,])
lines(sujs, e$acc*2000, col=2, lwd=2)



# Quel doit être le cut-off?
sum(data$rt[data$block!=2]>4000)/length(data$rt[data$block==2])



data <- droplevels(data[data$rt<4000 | data$block==2,])



# Blocs 1 et 4: tâche sur la cible centrale

e <- evalq(aggregate(list(acc=acc), list(hyp=hyp, bloc=BlockNum==1), mean),
           data[data$BlockNum %in% c(1,4),])
evalq(interaction.plot(bloc, hyp, acc, lwd=2, lty=1:3, pch=2:4, type='b'), e)
e <- evalq(aggregate(list(acc=acc), list(suj=suj, hyp=hyp, bloc=BlockNum==1), mean),
           data[data$BlockNum %in% c(1,4),])

l <- aov(acc~hyp*bloc+Error(suj/(bloc)), data=e); summary(l)

e <- evalq(aggregate(list(rt=rt), list(suj=suj, hyp=hyp, bloc=BlockNum==1), median),
           data[data$BlockNum %in% c(1,4),])
l <- aov(rt~hyp*bloc+Error(suj/(bloc)), data=e); summary(l)
evalq(interaction.plot(bloc, hyp, rt, lwd=2, lty=1:3, pch=2:4, type='b'), e)


par(mfrow=c(2,2))

evalq(plot(density(rt), col=2), data[evalq(hyp==2 & BlockNum==4, data),])
evalq(lines(density(rt)), data[evalq(hyp==1 & BlockNum==4, data),])
evalq(lines(density(rt), col=3), data[evalq(hyp==3 & BlockNum==4, data),])

e <- evalq(aggregate(list(rt=rt), list(suj=suj, hyp=hyp, bloc=BlockNum==1),
                     median),
           data[data$BlockNum %in% c(1,4),])

l <- aov(rt~hyp*bloc+Error(suj/(bloc)), data=e)
summary(l); model.tables(l, type='mean')

e <- evalq(aggregate(list(rt=rt), list(suj=suj, hyp=hyp==1),
                     sd),
           data[data$BlockNum ==4,])

l <- aov(rt~hyp, data=e)
evalq(t.test(rt[hyp], rt[!hyp]), e)

sujs <- levels(data$suj)
par(mfrow=c(4,6))
dump <- lapply(sujs, function(s){evalq(plot(density(rt), main=s),
                               data[evalq(suj==s & BlockNum!=2, data),])})



# Visibilité? Blocs 2 et 4
par(mfrow=c(2,2))
lapply(1:3, function(s){evalq(hist(Visibility, main=paste("hyp", s), breaks=0:4, ylim=c(0, 1),
                                   freq=FALSE),
                               data[evalq(hyp==s & BlockNum %in% c(2,4) & Visibility!=8, data),])})

sum(data$Visibility[data$BlockNum==4]==8)

e <- evalq(aggregate(list(vis=vis), list(hyp=hyp, dur=dur), mean),
           data[data$block==2,])
evalq(interaction.plot(dur, hyp, vis), e)

e <- evalq(aggregate(list(vis=vis), list(hyp=Hypnotiz, dur=dur),
                     function(x){sum(x>1)}),
           data[data$block %in% c(2,4),])
f <- evalq(aggregate(list(vis=vis), list(hyp=Hypnotiz, dur=dur),
                     function(x){sum(x==x)}),
           data[data$block %in% c(2,4),])
e$vis <- e$vis/f$vis

e <- evalq(aggregate(list(vis=vis), list(hyp=Hypnotiz, dur=dur),
                     mean),
           data[data$block %in% c(2,4),])
e <- evalq(aggregate(list(rtVis=rtVis), list(hyp=hyp, dur=dur),
                     mean),
           data[data$block %in% c(2,4),])

evalq(interaction.plot(dur, hyp, vis, col=1:4, lty=1:4, lwd=2, pch=2:5, type='b'), e)
evalq(interaction.plot(dur, hyp, rtVis, col=1:4, lty=1:4, lwd=2, pch=2:5, type='b'), e)



data$hiVis <- as.numeric(data$vis>1 )
data$loVis <- as.numeric(data$vis==1 )

l <- glmer(loVis~hyp*dur+(1|suj), family='binomial', data[data$block %in% c(2,4),])
summary(l)




# Categorization on the prime

e <- evalq(aggregate(list(acc=acc), list(hyp=hyp, dur=dur), mean),
           data[data$block %in% 3,])
evalq(interaction.plot(dur, hyp, acc, col=1:4, lty=1:4, lwd=2, pch=2:5, type='b'), e)

e <- evalq(aggregate(list(acc=acc), list(suj=suj, hyp=hyp, dur=dur), mean),
           data[data$block %in% 3,])
e$acc <- asin(e$acc)
l <- aov(acc~hyp*dur+Error(suj/(dur)), data=e); summary(l); model.tables(l, type='mean')

evalq(t.test(acc[dur==17], mu=.5), e)




e <- evalq(aggregate(list(rt=rt), list(suj=suj, hyp=hyp, dur=dur), median),
           data[data$block %in% 3,])
evalq(interaction.plot(dur, hyp, rt, lwd=2, pch=2:4, lty=2:4, type='b', col=1:3), e)
l <- aov(rt~dur*hyp+Error(suj/(dur)), data=e); summary(l)



e <- evalq(aggregate(list(rt=rt), list(hyp=hyp, prime=Prime), median),
           data[data$block %in% 3,])
evalq(interaction.plot(prime, hyp, rt, lwd=2, pch=2:4, lty=2:4, type='b'), e)



library(lme4)
data$scaledDur <- scale(data$dur)

l <- glmer(acc~hyp*scaledDur+(1|suj), family='binomial', data[data$block==3,])
summary(l)


# Priming? blocks: 1 et 4

# First do we have the distance effect?

e <- evalq(aggregate(list(rt=rt), list(hyp=hyp, target=Target), median),
           data[data$block %in% c(1,4),])
evalq(interaction.plot(target, hyp, rt, lwd=2, pch=2:4, lty=2:4, type='b'), e)

e <- evalq(aggregate(list(rt=rt), list(hyp=hyp, primed=IsPrimed), median),
           data[data$block %in% c(1,4) & data$IsPrimed!=2,])
evalq(interaction.plot(primed, hyp, rt, lwd=2, pch=2:4, lty=2:4, type='b'), e)
e <- evalq(aggregate(list(rt=rt), list(suj=suj, hyp=hyp, primed=IsPrimed), median),
           data[data$block %in% c(1,4) & data$IsPrimed!=2,])

l <- aov(rt~primed*hyp+Error(suj/(primed)), data=e); summary(l)
model.tables(l, type='means')


data$short <- data$dur %in% c(17,33)
e <- evalq(aggregate(list(rt=rt), list(dur=dur, primed=IsPrimed), median),
           data[data$block %in% c(1) & data$IsPrimed!=2,])
evalq(interaction.plot(primed, dur, rt, lwd=2, pch=2:7, lty=2:7, type='b'), e)

e <- evalq(aggregate(list(rt=rt), list(dur=dur), median),
           data[data$block %in% c(1,4) & data$IsPrimed!=2,])




e <- evalq(aggregate(list(rt=rt), list(dur=dur, hyp=hyp, primed=IsPrimed), median),
           data[data$block %in% c(1,4) & data$IsPrimed!=2,])

par(mfrow=c(2,2))
evalq(interaction.plot(primed, dur, rt, lwd=2, pch=2:6, lty=2:6, type='b', main="hyp: 1", ylim=c(500, 590)),
      e[e$hyp==1,])
evalq(interaction.plot(primed, dur, rt, lwd=2, pch=2:6, lty=2:6, type='b', main="hyp: 2", ylim=c(500, 590)), e[e$hyp==2,])
evalq(interaction.plot(primed, dur, rt, lwd=2, pch=2:6, lty=2:6, type='b', main="hyp: 3", ylim=c(500, 590)), e[e$hyp==3,])

e <- evalq(aggregate(list(rt=rt, acc=acc), list(vis=vis>1, primed=IsPrimed, hyp=hyp),
                     mean),
           data[data$block %in% c(4) & data$IsPrimed!=2,])

evalq(interaction.plot(hyp, primed, acc, lwd=2, pch=2:4, lty=2:4, type='b',
                       main='Seen'),
      e[e$vis,])
evalq(interaction.plot(hyp, primed, acc, lwd=2, pch=2:4, lty=2:4, type='b',
                       main='Not seen'),
            e[!e$vis,])


evalq(plot(rt, rtVis), data[data$block==4,])
l <- lm(rt~rtVis, data[data$block==4,])
summary(l)

e <- evalq(aggregate(rt, list(suj=suj, hyp=hyp, vis=vis>1, IsPrimed=IsPrimed), mean), data[data$block==4 & data$IsPrimed!=2,])
l <- aov(x~vis*hyp*IsPrimed+Error(suj/(vis*IsPrimed)), data=e)
summary(l)

l <- aov(x~hyp*IsPrimed+Error(suj/(IsPrimed)), data=e[e$vis==TRUE,])
summary(l)

e <- evalq(aggregate(acc, list(suj=suj, hyp=hyp, vis=vis>1, IsPrimed=IsPrimed), mean), data[data$block==4 & data$IsPrimed!=2,])
l <- aov(x~vis*hyp*IsPrimed+Error(suj/(vis*IsPrimed)), data=e)
summary(l)

l <- aov(x~hyp*IsPrimed+Error(suj/(IsPrimed)), data=e[e$vis==TRUE,])
summary(l)



# Effets sur les temps réponse en fonction de la visibilité subjective

e <- evalq(aggregate(list(rt=rt), list(hyp=hyp, hiVis=hiVis), median),
           data[data$BlockNum %in% c(1,4),])
evalq(interaction.plot(hiVis, hyp, rt, lwd=2, lty=1:3, pch=2:4, type='b'), e)

# Effets sur le temps de réponse sur la visibilité...

e <- evalq(aggregate(list(rtVis=rtVis), list(suj=suj, hyp=hyp, vis=vis), median),
           data[data$BlockNum %in% c(2,4),])
evalq(interaction.plot(vis, hyp, rtVis, lwd=2, lty=1:3, pch=2:4, type='b'), e)

l <- aov(rtVis~hyp*vis+Error(suj/(vis)), data=e); summary(l)
