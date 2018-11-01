load('data/HarvardData.R')

# We first exclude participants that we don't trust
dataRaw <- data
data <- droplevels(data[data$trust==1,])
N <- length(data$hf)

# Score Distribution un filtered
data$score <- apply(data[,1:12], 1, sum)
h <- hist(data$score, breaks=0:12)

h$counts/N
cumsum(h$counts/N)

h <- hist(data$score, breaks=c(0, 3, 7, 10, 12))
round(h$counts/N*100)
cumsum(h$counts/N)


# Score Distribution filtered

criterion <- 3
data$scoreFiltered <- apply(data[,1:12] & data[,15:26] > criterion, 1, sum)
h <- hist(data$scoreFiltered, breaks=0:12)

h$counts/N
cumsum(h$counts/N)

h <- hist(data$scoreFiltered, breaks=c(0, 2, 6, 9, 12))
h$counts/N*100
cumsum(h$counts/N)

filtered <- read.table(file="data/HypnoFilteredResponses.txt",
                   col.names=c('hf', 'ec', 'lhd', 'raf', 'fi', 'lar', 'mh', 'dhs', 'f', 'eo', 'la',
                               'amn', 'hyp', 'trust'))
filtered <- droplevels(filtered[filtered$trust==1,])
filtered$score <- apply(filtered[,1:12], 1, sum)

h <- hist(filtered$score, breaks=0:12)

h$counts/N
cumsum(h$counts/N)

h <- hist(data$scoreFiltered, breaks=c(0, 3, 7, 10, 12))
h$counts/N*100
cumsum(h$counts/N)

# Item pass rate

# Item scale correlation

library(lme4)

l <- glm(hf ~ I(ec+lhd+raf+ fi+ lar+ mh+ dhs+ f+ eo +la+amn), family='binomial', data=data)
summary(l)
l <- glm(amn ~ I(hf+ec+lhd+raf+ fi+ lar+ mh+ dhs+ f+ eo +la), family='binomial', data=data)
summary(l)

library(psych)

alpha(data[,1:12] & data[,15:26] > criterion)
