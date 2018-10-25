e <- evalq(aggregate(list(accOr=accOr, accn=accn, rtOr=rtOr), list(suj=suj, lag=lagC, allags=lag,
                                                                   target1=target1, target2=target2,
                                                                   hypnoC=hypnoC), mean),
           data[data$nface==2,])

f <- evalq(aggregate(list(accOr=accOr, accn=accn, rtOr=rtOr), list(lag=lag, allags=allags, target1=target1,
                                                                   target2=target2, hypnoC=hypnoC), mean), e)


l <- aov(accOr~lag*target1*target2*hypnoC+Error(suj/(lag*target1*target2)), data=e)
l <- aov(accOr~lag*target1+Error(suj/(lag*target1)), data=e)
l <- aov(accOr~allags*target2+Error(suj/(allags*target2)), data=e[e$allags==2 | e$allags==8,])

l <- aov(accOr~allags*target2+Error(suj/(allags*target2)), data=e[(e$allags==2 | e$allags==8) & e$hypnoC==1,])
l <- aov(accOr~allags*target2+Error(suj/(allags*target2)), data=e[(e$allags==2 | e$allags==8) & e$hypnoC==2,])
l <- aov(accOr~allags*target2+Error(suj/(allags*target2)), data=e[(e$allags==2 | e$allags==8) & e$hypnoC==3,])


summary(l)
model.tables(l, type='means')

evalq(interaction.plot(lag, target2, accOr),
      e[e$allags==2 | e$allags==8,])

evalq(interaction.plot(lag, target2, accOr),
      e[(e$allags==2 | e$allags==8) & e$hypnoC==1,])
evalq(interaction.plot(lag, target2, accOr),
      e[(e$allags==2 | e$allags==8) & e$hypnoC==2,])
evalq(interaction.plot(lag, target2, accOr),
      e[(e$allags==2 | e$allags==8) & e$hypnoC==3,])

par(mfrow=c(2,3))
pdf(file="BlinkPerHypno.pdf",width=12,height=6)
evalq(interaction.plot(lag, target2, accOr),
      e[e$allags==2 | e$allags==8,])

evalq(interaction.plot(lag, target2, accOr),
      e[(e$allags==2 | e$allags==8) & e$hypnoC==1,])
evalq(interaction.plot(lag, target2, accOr),
      e[(e$allags==2 | e$allags==8) & e$hypnoC==2,])
evalq(interaction.plot(lag, target2, accOr),
      e[(e$allags==2 | e$allags==8) & e$hypnoC==3,])
dev.off()