e <- evalq(aggregate(list(accOr=accOr, accn=accn, rtOr=rtOr), list(suj=suj, target1=target1, target2=target2, hypnoC=hypnoC), mean),
           data[data$nface==2,])
l <- aov(accOr~target1*target2*hypnoC+Error(suj/(target1*target2)), data=e)
summary(l)
model.tables(l, type='means')
evalq(interaction.plot(hypnoC, target1, accOr), e)
pdf("HypnoCTarget2Accor.pdf", height = 6, width = 12)
evalq(interaction.plot(hypnoC, target2, accOr), e)

dev.off()

AngryH1<-evalq(mean(accOr), e[e$target2=="angry" & e$hypnoC=="1",])
AngryH2<-evalq(mean(accOr), e[e$target2=="angry" & e$hypnoC=="2",])
AngryH3<-evalq(mean(accOr), e[e$target2=="angry" & e$hypnoC=="3",])
NeutralH1<-evalq(mean(accOr), e[e$target2=="neutral" & e$hypnoC=="1",])
NeutralH2<-evalq(mean(accOr), e[e$target2=="neutral" & e$hypnoC=="2",])
NeutralH3<-evalq(mean(accOr), e[e$target2=="neutral" & e$hypnoC=="3",])

AccPlot<-structure(list(Low=c(AngryH1, NeutralH1), Medium=c(AngryH2,NeutralH2), High=c(AngryH3,NeutralH3)))
AccPlot<-data.frame(AccPlot)
pdf("HypnoCTarget2AccorBAR.pdf", height = 6, width = 12)
barplot(as.matrix(AccPlot), ylim=c(0,1), beside=TRUE)
dev.off()


AngryH1<-evalq(mean(rtOr), e[e$target2=="angry" & e$hypnoC=="1",])
AngryH2<-evalq(mean(rtOr), e[e$target2=="angry" & e$hypnoC=="2",])
AngryH3<-evalq(mean(rtOr), e[e$target2=="angry" & e$hypnoC=="3",])
NeutralH1<-evalq(mean(rtOr), e[e$target2=="neutral" & e$hypnoC=="1",])
NeutralH2<-evalq(mean(rtOr), e[e$target2=="neutral" & e$hypnoC=="2",])
NeutralH3<-evalq(mean(rtOr), e[e$target2=="neutral" & e$hypnoC=="3",])

RTPlot<-structure(list(Low=c(AngryH1, NeutralH1), Medium=c(AngryH2,NeutralH2), High=c(AngryH3,NeutralH3)))
RTPlot<-data.frame(RTPlot)
pdf("HypnoCTarget2RTBAR.pdf", height = 6, width = 12)
barplot(as.matrix(RTPlot), ylim=c(0,1000), beside=TRUE)
dev.off()

pdf("HypnoCTarget2RTBARPROPER.pdf", height = 6, width = 12)
barplot(as.matrix(RTPlot), ylim=c(600,950), beside=TRUE)
dev.off()
