#+TITLE:Analyse blink hypnose 
#+AUTHOR: Jérôme Sackur
#+email: jerome.sackur@gmail.com
#+INFOJS_OPT: 
#+BABEL: :session *R* :cache yes :results output graphics :exports both :tangle yes 
-----

  # Import data and outlier cleanup

  #Don't run that again and again: load the data file (blinkData.R) instead.

#2016 - 03 - 18 

#One change compared to what we did on the 2016 - 03 - 17: Exclusion of
#participants based on slow global median RT is done if participant's
#median RT is 3 sd above the median of the medians of the group. This
#excludes only one participant. Basic pattern of results does not
#change.


#+BEGIN_SRC R :tangle yes
  data <- read.table(file="C:/Users/subcomarc/Google Drive/R/BlinkTable.txt", header=T)
  pilot <- read.table(file="C:/Users/subcomarc/Google Drive/R/TablePilotAB.txt", header=T)
  data <- droplevels(data[data$SubName!=101,])
  data$suj <- data$SubName
  sujs <- unique(data$suj); length(sujs)
  pilot$suj <- pilot$SubjName


  ntrial <- data$WhichTrial
  nblock <- data$WhichBlock

  data$target1 <- 'neutral'
  data$target1[data$FaceCombination== 3 | data$FaceCombination == 4] <- "angry"
  data$target2 <- 'neutral'
  data$target2[data$FaceCombination== 1 | data$FaceCombination == 4 | data$FaceCombination == 5] <- "angry"
  pilot$target1 <- 'neutral'
  pilot$target1[pilot$FaceCombination== 3 | pilot$FaceCombination == 4] <- "angry"
  pilot$target2 <- 'neutral'
  pilot$target2[pilot$FaceCombination== 1 | pilot$FaceCombination == 4 | pilot$FaceCombination == 5] <- "angry"


  data$nface <- 2
  data$nface[data$FaceCombination>4] <- 1

  data$accn <- data$FacesCorrect

  data$accOr <- data$AllAnswers


  pilot$nface <- 2
  pilot$nface[pilot$FaceCombination>4] <- 1

  pilot$accn <- pilot$FacesAcc

  pilot$accOr <- pilot$Accuracy

  data$rtOr <- data$RTimeWS*1000
  hist(data$rtOr[data$rtOr < 10000])

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
  limRt <- mean(e$medianRtOr)+3*sd(e$medianRtOr)

  bad <- as.vector(e$suj[e$medianRtOr > limRt])
  data <- droplevels(data[!data$suj %in% bad,])

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
  data$lagC <- factor(data$lagC)

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
  data$lagCatComplete <- factor(data$lagCatComplete)


  ### For the pilot


  e <- evalq(aggregate(list(accn=accn, accOr=accOr), list(suj=suj), mean), pilot)
  evalq(plot(accn, accOr, type='p', xlim=c(0, 1), ylim=c(0,1)), e)
  l <- lm(accOr~accn, e)
  abline(l)

  e$acc <- (e$accOr+e$accn)/2
  hist(e$acc); mean(e$acc); sd(e$acc)

  sum(e$acc < .6)
  bad <- as.vector(e$suj[e$acc < .6])

  e <- e[!e$suj %in% bad,]
  min(e$acc); hist(e$acc)

  pilot <- droplevels(pilot[!pilot$suj %in% bad,])

  pilot$lag <- pilot$Lags
  pilot$lagC <- "early"
  pilot$lagC[pilot$lag>=5] <- "late"
  pilot$lagC <- factor(pilot$lagC)

  ## Fake lag 0 for no T1

  pilot$lagComplete <- as.numeric(pilot$lag)
  pilot$lagComplete[pilot$nface==1] <- 0
  pilot$lagComplete <- ordered(pilot$lagComplete)

  pilot$lagCatComplete <- "early"
  pilot$lagCatComplete[pilot$lag>=5] <- "late"
  pilot$lagCatComplete[pilot$nface==1] <- "a_noT1" # a for alphebetical order for plotting purposes
  pilot$lagCatComplete <- factor(pilot$lagCatComplete)



  save(pilot, file='blinkPilot.R')
  save(data, file='blinkData.R')
#+END_SRC



#+BEGIN_SRC R
  load('blinkData.R')
  load('blinkPilot.R')
#+END_SRC


  # Basic plots
#+BEGIN_SRC R
  max(evalq(table(target1, target2, lagC, lag, suj), data))
  e <- evalq(mget.vaTer(rtOr, accOr, list(suj=suj, target1=target1, target2=target2, lag=lagC, allags=lag)), data)
  f <- evalq(aggregate(list(hypnoC=hypnoC), list(suj=suj), unique), data)
  e <- merge(e, f)

  l <- aov(v~lag*target1*target2*hypnoC+Error(suj/(lag*target1*target2)), data=e)
  summary(l)
  l <- aov(a~lag*target1*target2*hypnoC+Error(suj/(lag*target1*target2)), data=e)
  summary(l)


  e <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj, lag=lag), mean),
             data[data$nface==2,])
  f <- evalq(aggregate(list(accOr=accOr, accn=accn), list(lag=lag), mean), e)
  eprime <- evalq(aggregate(list(accOr=accOr, accn=accn, rtOr=rtOr), list(suj=suj), mean), data[data$nface==1,])
  eprime2 <- evalq(aggregate(list(accOr=accOr, accn=accn, rtOr=rtOr), list(suj=suj), mean), data[data$nface==2,])
  fprimeACCOR <- mean(eprime$accOr)
  fprimeRTOR <- mean(eprime$rtOr)
  fprimeSDACCOR<-sd(eprime$accOr)
  fprimeSDRTOR<-sd(eprime$rtOr)
  fACCOR <- mean(eprime2$accOr)
  fRTOR <- mean(eprime2$rtOr)
  fSDACCOR<- sd(eprime2$accOr)
  fSDRTOR<- sd(eprime2$rtOr)
  
  bartable<- rbind(fACCOR, fprimeACCOR)  ## get the cross tab
  pdf(file="MeanAccacrossPROPER.pdf",width=12,height=6)
  barplot(bartable, main="Accuracy across facenum", xlim=c(0,20), 
          ylab="Accuracy", ylim=c(0,1), beside = TRUE, legend = FALSE)  ## plot 
  dev.off()
  
  bartable<- rbind(fRTOR, fprimeRTOR)  ## get the cross tab
  pdf(file="MeanRTacrossPROPER.pdf",width=12,height=6)
  barplot(bartable, main="Mean RT across facenum", names.arg = c("RT","RT1Face"), xlim=c(0,20), 
          ylab="accuracy", ylim=c(600,900), beside = TRUE, legend = FALSE)  ## plot 
  dev.off()
  
  h<-evalq(aggregate(list(accOr=accOr, accn=accn, rtOr=rtOr), list(suj=suj, lag=lag, nface=nface, target1=target1, target2=target2), mean), 
           data)
  l<- aov(accOr~lag*nface*target1*target2+Error(suj/lag*target1*target2), data=h)
  l<- aov(rtOr~lag*nface*target1*target2+Error(suj/lag*target1*target2), data=h)
  summary(l)
  e <- evalq(aggregate(list(accOr=accOr, accn=accn, rtOr=rtOr), list(suj=suj, lag=lagC, allags=lag,
                                target1=target1, target2=target2,
                                hypnoC=hypnoC), mean),
             data[data$nface==2,])

  f <- evalq(aggregate(list(accOr=accOr, accn=accn, rtOr=rtOr), list(lag=lag, allags=allags, target1=target1,
                                             target2=target2, hypnoC=hypnoC), mean), e)

  
  l <- aov(accOr~lag*target1*target2*hypnoC+Error(suj/(lag*target1*target2)), data=e)
  summary(l)
  model.tables(l, type='means')

  ##mine
  e <- evalq(aggregate(list(accOr=accOr), list(suj=suj, nface=nface, target1=target1, target2=target2,
                                                                     hypnoC=hypnoC), mean),
            #data[data$nface==2,])
              data)
e <- evalq(aggregate(list(accOr=accOr, accn=accn, rtOr=rtOr), list(suj=suj, lag=lagC, nface=nface, allags=lag,
                                                                   target1=target1, target2=target2,
                                                                   hypnoC=hypnoC), mean),
           data[data$nface==1,])
f <- evalq(aggregate(list(accOr=accOr, accn=accn, rtOr=rtOr), list(lag=lag, allags=allags,
                                                                   target2=target2, hypnoC=hypnoC), mean), e)
  f <- evalq(aggregate(list(accOr=accOr, accn=accn, rtOr=rtOr), list(lag=lag, allags=allags, target1=target1,
                                                                     target2=target2, hypnoC=hypnoC), mean), e)
  l <- aov(accOr~lag*target2*hypnoC+Error(suj/(lag*target2)), data=e)
  summary(l)
l <- aov(rtOr~lag*target2*hypnoC+Error(suj/(lag*target2)), data=e)
summary(l)
  
hh <- aov(accOr~nface*hypnoC*target2+Error(suj/(target2*nface)), data=e)
summary(hh)

  tt <- aov(rtOr~allags*target1*target2*hypnoC+Error(suj/(lag*target1*target2)), data=e)
  summary(tt)

  evalq(interaction.plot(hypnoC,nface,accOr), e)
  evalq(interaction.plot(hypnoC,target2,accOr), e)
  evalq(interaction.plot(target2, lag, rtOr), e)
  
  library('xtable')
  library('knitr')
  options(xtable.floating = FALSE)
  options(xtable.timestamp = "")
  
  xtable(l)
  xtable(tt)
  model.tables(tt, type='means')
  h<-model.tables(l, type='means')
  
  ##mine
  
  
  par(mfrow=c(2,2))
  evalq(interaction.plot(hypnoC, target2, accOr), f)
  evalq(interaction.plot(target1, target2, accOr), f)
  evalq(interaction.plot(lag, target2, rtOr), f)


  par(mfrow=c(3,2))

 ##mine
  ##AccOr
  pdf("AccOrT1_T2ALL.pdf", width=12, height=6)
  evalq(interaction.plot(allags, target1, accOr, main="Acc OR T1 (T2 all) All H", ylim=c(.6, .90)),
        f)
  dev.off()
  pdf("AccOrT1_T2angry.pdf", width=12, height=6)
  evalq(interaction.plot(allags, target1, accOr, main="Acc OR T1 (T2 angry) All H", ylim=c(.6, .90)),
        f[f$target2=="angry",])
  dev.off()
  pdf("AccOrT1_T2neutral.pdf", width=12, height=6)
  evalq(interaction.plot(allags, target1, accOr, main="Acc OR T1 (T2 neutral) All H", ylim=c(.6, .90)),
        f[f$target2=="neutral",])
  dev.off()
  pdf("AccOrT2_T1ALL.pdf", width=12, height=6)
  evalq(interaction.plot(allags, target2, accOr, main="Acc OR T2 (T1 all) All H", ylim=c(.6, .90)),
        f)
  dev.off()
  pdf("AccOrT2_T1angry.pdf", width=12, height=6)
  evalq(interaction.plot(allags, target2, accOr, main="Acc OR T2 (T1 angry) All H", ylim=c(.6, .90)),
        f[f$target1=="angry",])
  dev.off()
  pdf("AccOrT2_T1neutral.pdf", width=12, height=6)
  evalq(interaction.plot(allags, target2, accOr, main="Acc OR T2 (T1 neutral) All H", ylim=c(.6, .90)),
        f[f$target1=="neutral",])
  dev.off()
  
  ##RTOr
  pdf("RTOrT1_T2ALL.pdf", width=12, height=6)
  evalq(interaction.plot(allags, target1, rtOr, main="RT OR T1 (T2 all) All H", ylim=c(600,1000)),
        f)
  dev.off()
  pdf("RTOrT1_T2angry.pdf", width=12, height=6)
  evalq(interaction.plot(allags, target1, rtOr, main="RT OR T1 (T2 angry) All H", ylim=c(600,1000)),
        f[f$target2=="angry",])
  dev.off()
  pdf("RTOrT1_T2neutral.pdf", width=12, height=6)
  evalq(interaction.plot(allags, target1, rtOr, main="RT OR T1 (T2 neutral) All H", ylim=c(600,1000)),
        f[f$target2=="neutral",])
  dev.off()
  pdf("RTOrT2_T1ALL.pdf", width=12, height=6)
  evalq(interaction.plot(allags, target2, rtOr, main="RT OR T2 (T1 all) All H", ylim=c(600,1000)),
        f)
  dev.off()
  pdf("RTOrT2_T1angry.pdf", width=12, height=6)
  evalq(interaction.plot(allags, target2, rtOr, main="RT OR T2 (T1 angry) All H", ylim=c(600,1000)),
        f[f$target1=="angry",])
  dev.off()
  pdf("RTOrT2_T1neutral.pdf", width=12, height=6)
  evalq(interaction.plot(allags, target2, rtOr, main="RT OR T2 (T1 neutral) All H", ylim=c(600,1000)),
        f[f$target1=="neutral",])
  dev.off()
  
  
  evalq(interaction.plot(allags, target1, accOr, main="Acc OR T1 (T2 angry) All H", ylim=c(.6, .90)),
        f[f$target2=="angry",])
  evalq(interaction.plot(lag, target2, accOr, main="Acc OR T1 (T2 angry) All H"),
        f[data$nface==1,])
  evalq(interaction.plot(lag, target2, rtOr, main="Acc OR T1 (T2 angry) All H"),
        f[data$nface==1,])
  
  
  evalq(interaction.plot(allags, target2, rtOr, main="ACC OR T2 angry - LOW"),
        data[data$nface==1,])
  evalq(plot(allags,rtOr),e)
  evalq(plot(allags,rtOr),e)
  
  
##mine 
  
  evalq(interaction.plot(allags, target1, accOr, main="ACC OR T2 angry - LOW", ylim=c(.6, .90)),
        f[f$target2=="angry" & f$hypnoC==1,])
  evalq(interaction.plot(lag, target1, accOr, main="ACC OR T2 neutral - LOW", ylim=c(.6, .90)),
        f[f$target2=="neutral" & f$hypnoC==1,])

  evalq(interaction.plot(lag, target1, accOr, main="ACC OR T2 angry - MEDIUM", ylim=c(.6, .90)),
        f[f$target2=="angry" & f$hypnoC==2,])
  evalq(interaction.plot(lag, target1, accOr, main="ACC OR T2 neutral - MEDIUM", ylim=c(.6, .90)),
        f[f$target2=="neutral" & f$hypnoC==2,])

  evalq(interaction.plot(lag, target1, accOr, main="ACC OR T2 angry - HIGH", ylim=c(.6, .90)),
        f[f$target2=="angry" & f$hypnoC==3,])
  evalq(interaction.plot(lag, target1, accOr, main="ACC OR T2 neutral - HIGH", ylim=c(.6, .90)),
        f[f$target2=="neutral" & f$hypnoC==3,])


  e <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj, lag=lagC, target1=target1,
                                             target2=target2, hypnoC=hypnoC), mean),
             data[data$nface==2 & data$target2=='neutral',])
  f <- evalq(aggregate(list(accOr=accOr, accn=accn), list(lag=lag, target1=target1,
                                             target2=target2, hypnoC=hypnoC), mean), e)
  par(mfrow=c(2,3))
  evalq(interaction.plot(lag, target1, accOr, main="ACC OR T2 Neutral - LOW", ylim=c(.6, .90)),
        f[f$hypnoC==1,])
  evalq(interaction.plot(lag, target1, accOr, main="ACC OR T2 Neutral - Medium", ylim=c(.6, .90)),
        f[f$hypnoC==2,])
  evalq(interaction.plot(lag, target1, accOr, main="ACC OR T2 Neutral - high", ylim=c(.6, .90)),
        f[f$hypnoC==3,])
  l <- aov(accOr~lag*target1*hypnoC+Error(suj/(lag*target1)), data=e)
  summary(l)

  e <- evalq(aggregate(list(rt=rtOr), list(suj=suj, lag=lag, target1=target1,
                                             target2=target2, hypnoC=hypnoC), median),
             data[data$nface==2,])
  f <- evalq(aggregate(list(rt=rt), list(lag=lag, target1=target1,
                                             target2=target2, hypnoC=hypnoC), mean), e)

  l <- aov(rt~lag*target1*target2*hypnoC+Error(suj/(lag*target1*target2)), data=e)
  summary(l)
  model.tables(l, type='means')
  pdf("T2HypRT.pdf", width = 12, height = 6)
  evalq(interaction.plot(hypnoC, target2, rt), f)
  dev.off()

  par(mfrow=c(3,2))
  evalq(interaction.plot(lag, target1, rt, main="ACC OR T2 angry - LOW", ylim=c(550, 1050)),
        f[f$target2=="angry" & f$hypnoC==1,])
  evalq(interaction.plot(lag, target1, rt, main="ACC OR T2 neutral - LOW", ylim=c(550, 1050)),
        f[f$target2=="neutral" & f$hypnoC==1,])

  evalq(interaction.plot(lag, target1, rt, main="ACC OR T2 angry - MEDIUM", ylim=c(550, 1050)),
        f[f$target2=="angry" & f$hypnoC==2,])
  evalq(interaction.plot(lag, target1, rt, main="ACC OR T2 neutral - MEDIUM", ylim=c(550, 1050)),
        f[f$target2=="neutral" & f$hypnoC==2,])

  evalq(interaction.plot(lag, target1, rt, main="ACC OR T2 angry - HIGH", ylim=c(550, 1050)),
        f[f$target2=="angry" & f$hypnoC==3,])
  evalq(interaction.plot(lag, target1, rt, main="ACC OR T2 neutral - HIGH", ylim=c(550, 1050)),
        f[f$target2=="neutral" & f$hypnoC==3,])

  data$logRt <- log(data$rtOr)
  e <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj, lag=lag, target1=target1,
                                             target2=target2, hypnoC=hypnoC), mean),
             data[data$nface==2,])
  f <- evalq(aggregate(list(accOr=accOr, accn=accn), list(lag=lag, target1=target1,
                                             target2=target2, hypnoC=hypnoC), mean), e)

  ert <- evalq(aggregate(list(rt=logRt), list(suj=suj, lag=lag, target1=target1,
                                             target2=target2, hypnoC=hypnoC), median),
             data[data$nface==2,])
  frt <- evalq(aggregate(list(rt=rt), list(lag=lag, target1=target1,
                                             target2=target2, hypnoC=hypnoC), mean), ert)

  e <- merge(e,ert)
  e$ies <- e$rt/e$accOr

  par(mfrow=c(3,2))
  evalq(interaction.plot(lag, target1, ies, main="ACC OR T2 angry - LOW", ylim=c(7, 11.5)),
        e[e$target2=="angry" & e$hypnoC==1,])
  evalq(interaction.plot(lag, target1, ies, main="ACC OR T2 neutral - LOW", ylim=c(7, 11.5)),
        e[e$target2=="neutral" & e$hypnoC==1,])

  evalq(interaction.plot(lag, target1, ies, main="ACC OR T2 angry - MEDIUM", ylim=c(7, 11.5)),
        e[e$target2=="angry" & e$hypnoC==2,])
  evalq(interaction.plot(lag, target1, ies, main="ACC OR T2 neutral - MEDIUM", ylim=c(7, 11.5)),
        e[e$target2=="neutral" & e$hypnoC==2,])

  evalq(interaction.plot(lag, target1, ies, main="ACC OR T2 angry - HIGH", ylim=c(7, 11.5)),
        e[e$target2=="angry" & e$hypnoC==3,])
  evalq(interaction.plot(lag, target1, ies, main="ACC OR T2 neutral - HIGH", ylim=c(7, 11.5)),
        e[e$target2=="neutral" & e$hypnoC==3,])


  l <- aov(ies~lag*target1*target2*hypnoC+Error(suj/(lag*target1*target2)), data=e[e$hypnoC!=2,])
  summary(l)
  model.tables(l, type='means')

  evalq(interaction.plot(hypnoC, target2, ies), e)



  l <- aov(accn~lag*target1*target2*hypnoC+Error(suj/(lag*target1*target2)), data=e)
  summary(l)
  model.tables(l, type='means')

  evalq(interaction.plot(hypnoC, target2, accOr), e)

  par(mfrow=c(3,2))
  evalq(interaction.plot(lag, target1, accn, main="ACC N T2 angry - LOW", ylim=c(.4, 1)),
        f[f$target2=="angry" & f$hypnoC==1,])
  evalq(interaction.plot(lag, target1, accn, main="ACC N T2 neutral - LOW", ylim=c(.4, 1)),
        f[f$target2=="neutral" & f$hypnoC==1,])
  evalq(interaction.plot(lag, target1, accn, main="ACC N T2 angry - MEDIUM", ylim=c(.4, 1)),
        f[f$target2=="angry" & f$hypnoC==2,])
  evalq(interaction.plot(lag, target1, accn, main="ACC N T2 neutral - MEDIUM", ylim=c(.4, 1)),
        f[f$target2=="neutral" & f$hypnoC==2,])
  evalq(interaction.plot(lag, target1, accn, main="ACC N T2 angry - HIGH", ylim=c(.4, 1)),
        f[f$target2=="angry" & f$hypnoC==3,])
  evalq(interaction.plot(lag, target1, accn, main="ACC N T2 neutral - HIGH", ylim=c(.4, 1)),
        f[f$target2=="neutral" & f$hypnoC==3,])




  evalq(interaction.plot(lag, target1, accn, main="ACC N T2 angry", ylim=c(.4, 1)), f[f$target2=="angry",])
  evalq(interaction.plot(lag, target1, accn, main="ACC NT2 neutral", ylim=c(.4, 1)), f[f$target2=="neutral",])

  e <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj, lag=lag), mean),
             pilot[pilot$nface==2,])
  f <- evalq(aggregate(list(accOr=accOr, accn=accn), list(lag=lag), mean), e)
  eprime <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj), mean), pilot[pilot$nface==1,])
  fprime <- mean(eprime$accOr)

  e <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj, lag=lag, target1=target1, target2=target2), mean),
             pilot[pilot$nface==2,])
  f <- evalq(aggregate(list(accOr=accOr, accn=accn), list(lag=lag, target1=target1, target2=target2), mean), e)


  par(mfrow=c(2,2))
  evalq(interaction.plot(lag, target1, accOr, main="ACC OR T2 angry", ylim=c(.6, .90)), f[f$target2=="angry",])
  evalq(interaction.plot(lag, target1, accOr, main="ACC OR T2 neutral", ylim=c(.6, .90)), f[f$target2=="neutral",])

  evalq(interaction.plot(lag, target1, accn, main="ACC N T2 angry", ylim=c(.4, 1)), f[f$target2=="angry",])
  evalq(interaction.plot(lag, target1, accn, main="ACC NT2 neutral", ylim=c(.4, 1)), f[f$target2=="neutral",])

#+END_SRC
  

  # Basic anovas

#+BEGIN_SRC R

  e <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj, lag=lag, target1=target1, target2=target2), mean),
             data[data$nface==2,])
  l <- aov(accOr~lag*target1*target2+Error(suj/(lag*target1*target2)), data=e)
  summary(l)
  evalq(interaction.plot(target1, target2, accOr), e)
  # So there is a full recovery of the blink with an angry T2. We can
  # thus look for blink effects in the subset of the data with T2
  # neutral

  e <- evalq(aggregate(list(rt=rtOr, accOr=accOr), list(suj=suj, lag=lag, target1=target1, target2=target2), mean),
             data[data$nface==2,])
  l <- aov(rt~lag*target1*target2+Error(suj/(lag*target1*target2)), data=e)
  summary(l)
  evalq(interaction.plot(target1, target2, accOr), e)
  # So there is a full recovery of the blink with an angry T2. We can

  data$hypnoC <- ordered(data$hypnoC)

  e <- evalq(aggregate(list(rt=rtOr), list(suj=suj, lag=lag, target1=target1,
                                             hypnoC=hypnoC), median),
             data[data$nface==2 & data$target2=="neutral",])
  
  e <- evalq(aggregate(list(rt=rtOr, accOr=accOr), list(suj=suj, lag=lag, target1=target1, target2=target2,
                                           hypnoC=hypnoC), median),
             data[data$nface==2,])


  l <- aov(rt~lag*target1*hypnoC+Error(suj/(lag*target1)), data=e)
  l <- aov(accOr~lag*target2*hypnoC+Error(suj/(lag*target2)), data=e)
  summary(l)
  evalq(interaction.plot(hypnoC, target2, rt), e)

  e <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj, target1=target1,
                                             hypnoC=hypnoC), mean),
             data[data$nface==2 & data$target2=="neutral",])


  l <- aov(accn~target1*hypnoC+Error(suj/(target1)), data=e)
  summary(l)
  evalq(interaction.plot(hypnoC, target1, accOr), e)

  e <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj, target1=target1, target2=target2,
                                             hypnoC=hypnoC), mean),
             data)


  l <- aov(accn~target2*target1*hypnoC+Error(suj/(target2*target1)), data=e)
  summary(l)
  evalq(interaction.plot(hypnoC, target1, accn), e)





  e <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj, target2=target2,
                                             hypnoC=hypnoC), mean),
             data)

  l <- aov(accn~target2*hypnoC+Error(suj/(target2)), data=e)
  summary(l)
  evalq(interaction.plot(hypnoC, target2, accOr), e)

  contrasts(data$target1)<-contr.sum(2)
  data$hypnoC <- ordered(data$hypnoC)

  data$hyp <- ordered(data$HLRAWcol)

  l <- glmer(accOr~target1*hyp+(1|suj), data=data, family='binomial')
  summary(l)


  ## For the plots, we add the "lag 0 = no T1" of lagComplete
  e <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj, lag=lagCatComplete,
                                             target1=target1, hypnoC=hypnoC), mean),
             data[data$target2=="neutral",])
  par(mfrow=c(2,2))
  evalq(interaction.plot(lag, target1, accOr, main='low', ylim=c(.4, 1)), e[e$hypnoC==1,])
  evalq(interaction.plot(lag, target1, accOr, main='medium', , ylim=c(.4, 1)), e[e$hypnoC==2,])
  evalq(interaction.plot(lag, target1, accOr, main='high', , ylim=c(.4, 1)), e[e$hypnoC==3,])



  e <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj, target1=target1, target2=target2, hypnoC=hypnoC), mean),
             data[data$nface==2,])
  l <- aov(accOr~target1*target2*hypnoC+Error(suj/(target1*target2)), data=e)
  summary(l)
  model.tables(l, type='means')
  evalq(interaction.plot(hypnoC, target2, accOr), e)
  evalq(interaction.plot(hypnoC, target1, accOr), e)

  #mine
  
  bartable<- evalq(rbind(target2, hypnoC, accOr), e)  ## get the cross tab
  pdf(file="MeanAccacrossPROPER.pdf",width=12,height=6)
  barplot(bartable, main="Accuracy across facenum", xlim=c(0,20), 
          ylab="Accuracy",  beside = TRUE, legend = TRUE) ## plot 
  dev.off()
  
  e <- evalq(aggregate(list(accOr=accOr, accn=accn, rtOr=rtOr), list(suj=suj, lag=lag, target1=target1, target2=target2, hypnoC=hypnoC), mean),
             data[data$nface==2,])
  l <- aov(rtOr~target1*target2*lag*hypnoC+Error(suj/(target1*target2*lag)), data=e)
  

  
  summary(l)
  
  
  #mine
  

#+END_SRC

  # We do the same on the pilot and compare pilot and data

#+BEGIN_SRC R
  pilot$target1 <- factor(pilot$target1)
  pilot$target2 <- factor(pilot$target2)
  pilot$lag <- ordered(pilot$lag)
  pilot$suj <- factor(pilot$suj)


  e <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj, lag=lag, target1=target1, target2=target2), mean),
             pilot[pilot$nface==2,])
  l <- aov(accOr~lag*target1*target2+Error(suj/(lag*target1*target2)), data=e)
  summary(l)
  evalq(interaction.plot(target1, target2, accOr), e)
  # So there is a full recovery of the blink with an angry T2. We can
  # thus look for blink effects in the subset of the pilot with T2
  # neutral



  e <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj, lag=lagC, target1=target1,
                                                          target2=target2), mean),
             pilot[pilot$nface==2,])

  l <- aov(accOr~lag*target1*target2+Error(suj/(lag*target1*target2)), data=e)
  summary(l)
  evalq(interaction.plot(lag, target1, accOr), e)


  ## For the plots, we add the "lag 0 = no T1" of lagComplete
  e <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj, lag=lagCatComplete,
                                             target1=target1), mean),
             pilot[pilot$target2=="neutral",])
  par(mfrow=c(2,2))
  evalq(interaction.plot(lag, target1, accOr, main='No hypnosis', ylim=c(.6, .9)), e)

  f <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj, lag=lagCatComplete,
                                             target1=target1), mean),
             data[data$target2=="neutral" & data$hypnoC==3,])


  evalq(interaction.plot(lag, target1, accOr, main='Hypnosis, Highs', ylim=c(.6, .9)), f)


  e <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj, lag=lagC,
                                             target1=target1), mean),
             pilot[pilot$target2=="neutral" & pilot$nface==2,])

  f <- evalq(aggregate(list(accOr=accOr, accn=accn), list(suj=suj, lag=lagC,
                                             target1=target1), mean),
             data[data$target2=="neutral" & data$hypnoC==3 & data$nface==2,])

  e$condition <- "noHyp"
  f$condition <- "Hyp"
  g <- rbind(e,f)

  l <- aov(accOr~lag*target1*condition+Error(suj/(lag*target1)), data=g)
  summary(l)
  l <- aov(accOr~lag*target1+Error(suj/(lag*target1)), data=g[g$condition=="noHyp",])
  summary(l)
  l <- aov(accOr~lag*target1+Error(suj/(lag*target1)), data=g[g$condition=="Hyp",])
  summary(l)

  library(lme4)

  e <- evalq(data.frame(suj=suj, lag=lagC, target1=target1, accOr=accOr, condition='Hyp'),
             data[data$target2=="neutral" & data$hypnoC==3 & data$nface==2,])
  f <- evalq(data.frame(suj=suj, lag=lagC, target1=target1, accOr=accOr, condition='noHyp'),
             pilot[data$target2=="neutral" & data$nface==2,])
  g <- rbind(e,f)
  contrasts(g$lag)<-contr.sum(2)
  contrasts(g$target1)<-contr.sum(2)
  contrasts(g$condition)<-contr.sum(2)

  l <- glmer(accOr~lag*target1*condition+(1|suj), family='binomial', data=g)
  summary(l)

  l <- glmer(accOr~lag*target1+(1|suj), family='binomial', data=e)
  summary(l)

  l <- glmer(accOr~lag*target1+(1|suj), family='binomial', data=f)
  summary(l)

#+END_SRC

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



# Diffusion models

#+BEGIN_SRC R
  load("blinkData.R")

  ## First we plot densities


  densLCorr <- density(data$rtOr[data$accOr==1 & data$hypnoC==1])
  densMCorr <- density(data$rtOr[data$accOr==1 & data$hypnoC==2])
  densHCorr <- density(data$rtOr[data$accOr==1 & data$hypnoC==3])
  densLInc <- density(data$rtOr[data$accOr==0 & data$hypnoC==1])
  densMInc <- density(data$rtOr[data$accOr==0 & data$hypnoC==2])
  densHInc <- density(data$rtOr[data$accOr==0 & data$hypnoC==3])

  M <- max(c(densLCorr$y,densMCorr$y,densHCorr$y, densLInc$y,densMInc$y,densHInc$y))
pdf("Densities.pdf", width = 12, height = 6)
  plot(densLCorr, col=2, lwd=3, ylim=c(0, M), xlim=c(0, 3500))
  lines(densMCorr, col=1, lwd=3); lines(densHCorr, col=3, lwd=3)
  lines(densLInc, col=2, lwd=1, lty=2); lines(densMInc, col=1, lwd=1, lty=2); lines(densHInc, col=3, lwd=1, lty=2)
dev.off()
  ## ## Okay. So it kind of seem to make sense to do diffusion models.


  ## sujs <- unique(data$suj);
  ## i <- 1
  ## df <- data
  ## ## We recode t1 et t2 levels, because it makes things easier
  ## ## aftewards.
  ## df$T1 <- paste(df$target1, "t1", sep='')
  ## df[df$nface==1,]$T1 <- 'absent'
  ## df$T2 <- paste(df$target2, "t2", sep='')

  ## ## We write the date for each subject in different files, in different
  ## ## folder so as to parallelize by hand.
  ## for (s in sujs){
  ##     folder <- paste("diffusion/d", ceiling(i/12), "/", sep='')
  ##     temp <- evalq(data.frame(acc=accOr, rt=rtOr/1000, target1=T1,
  ##                              target2=T2, lag=lagC),
  ##                   df[df$suj==s,])
  ##     write.table(temp, file=paste(folder, s, ".txt", sep=''),
  ##                 row.names=FALSE, col.names=FALSE, quote=FALSE)
  ##     i <- i+1
  ## }



  ## source("/home/jerome/science/0librairies/R_Lib_Sackur.r")

  ## fdm <- do.call(rbind, lapply(1:4, function(n){
  ##              fdmImp <- readFastDMComplete(file=paste("diffusion/d", n,
  ##                                               '/blinkFinal.out', sep=''),
  ##                                           dvs=c('v','a','t0'),
  ##                                           factors=c('target1', 'target2', 'lag'),
  ##                                           flevels=list(c('absent', 'angryt1', 'neutralt1'),
  ##                                               c('angryt2', 'neutralt2'),
  ##                                               c('early', 'late')))
  ##          }))


  ## f <- evalq(aggregate(list(hypnoC=hypnoC), list(suj=suj), unique), data)
  ## fdm <- merge(fdm, f)
  ## save(fdm, file="Fast-dm_diffusion_parameters.R")

  load("Fast-dm_diffusion_parameters.R")

  l1 <- aov(v~target1*target2*hypnoC*lag+Error(suj/(target1*target2*lag)), data=fdm)
  
  
  
  summary(l1)
  pdf("DiffVLagT1.pdf", width = 12, height = 6)
  evalq(interaction.plot(lag, target1, v), fdm)
  dev.off()
  model.tables(l1, type='means')

  l2 <- aov(a~target1*target2*hypnoC*lag+Error(suj/(target1*target2*lag)), data=fdm)
  summary(l2)
  model.tables(l1, type='means')

  afdm <- evalq(aggregate(list(a=a), list(target1=target1, target2=target2, lag=lag), mean), fdm)
  par(mfrow=c(2,2))
  evalq(plot(a~target1), afdm); 
  evalq(plot(a~target2), afdm); 

  pdf("DiffAT2Hypno.pdf", width = 12, height = 6)
  evalq(interaction.plot(hypnoC, target2, a), fdm)
  dev.off()
  evalq(interaction.plot(target1, target2, a), fdm)
  evalq(interaction.plot(lag, target2, a), fdm)

  l2 <- aov(a~target1*target2*lag*hypnoC+Error(suj/(target1*target2*lag)),
           data=fdm[fdm$hypnoC!=3,])
  summary(l2)
  evalq(interaction.plot(hypnoC, target2, a), fdm)


  l3 <- aov(t0~target1*target2*lag*hypnoC+Error(suj/(target1*target2*lag)), data=fdm)
  summary(l3)
  evalq(interaction.plot(lag, target2, t0), fdm)
  evalq(interaction.plot(lag, target1, t0), fdm)
  evalq(interaction.plot(hypnoC, target1, t0), fdm)
  

  model.tables(l, type='means')

  library('xtable')
  library('knitr')
  options(xtable.floating = FALSE)
  options(xtable.timestamp = "")
  
  xtable(l1)
  xtable(l2)
  xtable(l3)
  
  model.tables(tt, type='means')
  h<-model.tables(l, type='means')
  
  
  
  par(mfrow=c(2,2))
  evalq(interaction.plot(hypnoC, target2, t0,
                         main="Early Angry T1", ylim=c(.22, .32)),
        fdm[evalq(lag=="early" & target1=="angryt1", fdm),])
  evalq(interaction.plot(hypnoC, target2, t0,
                         main="Early Neutral T1", ylim=c(.22, .32)),
        fdm[evalq(lag=="early" & target1=="neutralt1", fdm),])
  evalq(interaction.plot(hypnoC, target2, t0,
                         main="Late Angry T1", ylim=c(.22, .32)),
        fdm[evalq(lag=="late" & target1=="angryt1", fdm),])
  evalq(interaction.plot(hypnoC, target2, t0,
                         main="Late Neutral T1", ylim=c(.22, .32)),
        fdm[evalq(lag=="late" & target1=="neutralt1", fdm),])




#+END_SRC
