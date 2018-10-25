data <- read.table(file="DataTable.csv", header=T)
data$suj <- factor(data$SubjName)
data$hyp <- round(data$SubjName/1000)
data$rt <- round(numerize(data$RTime)*1000)
data$acc <- data$Correctness
data$vis <- data$Visibility
data$block <- data$BlockNum
data$dur <- round(data$PrimeDur * 1000)
data$rtVis <- round(numerize(data$RTimeVis)*1000)

# Erreur de codage
data$Hypnotiz[data$suj==3051] <- 3


# Erreur de réponse: les sujets confondent les touches de visibilité
# et de catégorisation
data <- droplevels(data[data$vis!=8,])
data <- droplevels(data[data$rt<3000 | data$block==2,])



# Blocs 1 et 4: tâche sur la cible centrale

e <- evalq(aggregate(list(acc=acc), list(hyp=hyp, bloc=BlockNum==1), mean),
           data[data$BlockNum %in% c(1,4),])

e <- evalq(aggregate(list(rt=rt), list(hyp=hyp, bloc=BlockNum==1), median),
           data[data$BlockNum %in% c(1,4),])

par(mfrow=c(2,2))

evalq(interaction.plot(bloc, hyp, rt, lwd=2, lty=1:3, pch=2:4, type='b'), e)

evalq(plot(density(rt)), data[evalq(hyp==1 & BlockNum==4, data),])
evalq(lines(density(rt), col=2), data[evalq(hyp==2 & BlockNum==4, data),])
evalq(lines(density(rt), col=3), data[evalq(hyp==3 & BlockNum==4, data),])

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
