plotMorey <- function(MainFrame, DF, x, part, xFactor1, xFactor2, Morey, errpos,...) 
  {
  MainFrame=data
  DF=e
  x=e$accOr
  datapart=data$accOr
  part=data$suj
  xFactor1=e$lagC
  xFactor2=e$target2
  Morey=3
  errpos=0.001
  
  grandMean <- mean(x)
  partMeans <-evalq(aggregate(list(partMean=datapart), list(part=part), mean),
                    MainFrame)
  DF <- merge(DF, partMeans)
  DF$nax <- DF$x-DF$partMean+grandMean
  f <- evalq(aggregate(list(sex=nax), list(xFactor1=xFactor1, xFactor2=xFactor2), se), DF)
  f$sex <- f$sex*sqrt/(Morey) # (6/3) compute se and apply Morey's 2008 correction (For a 2X2 design M=4).
  f$mx <- evalq(aggregate(list(mx=x), list(xFactor1=xFactor1, xFactor2=xFactor2), mean), DF)$mx
  coordmax<-length(levels(DF$xFactor1))
  coordmax<-c(1:coordmax)
  xs <- coordmax-errpos
  y1=f$mx+f$sex
  y0=f$mx-f$sex
  pdf("myplot.pdf")
  evalq(interaction.plot(xFactor1, xFactor2, x) ,DF)
  arrows(x0=xs, y0=y0, x1=xs, y1=y1, length=0,
         angle=90,code=3, lty=1, lwd=2)
  dev.off()  
}

