load('C:/Users/lscpuser/Downloads/HarvardData.R')

# We first exclude participants that we don't trust
dataRaw <- data
data <- droplevels(data[data$trust==1,])
N <- length(data$hf)

###JEROME'S CODE##
###JEROME'S CODE##
###JEROME'S CODE##
# Score Distribution un filtered
data$score <- apply(data[,1:12], 1, sum)
h <- hist(data$score, breaks=0:12)

h$counts/N
cumsum(h$counts/N)

h <- hist(data$score, breaks=c(0, 3, 7, 10, 12))
round(h$counts/N*100)
cumsum(h$counts/N)


# Score Distribution filtered

criterion <- 2
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
###JEROME'S CODE##
###JEROME'S CODE##
###JEROME'S CODE##

# Check if the mean sample difference in hypnotizability between filtered and unfiltered scores is significant

PureScores<-data$score
FilteredScores<-data$scoreFiltered
Scores<-c(FilteredScores, PureScores)
MeansFiltered <-mean(FilteredScores)                     # Mean sample hyp for filtered and unfiltered
SDFiltered<-sd(FilteredScores)
MeanPure<-mean(PureScores)
SDPure<-sd(PureScores)
IsFiltered<-c(rep("Filtered",length(data$scoreFiltered)), rep("NotFiltered",length(data$score))) #Is the difference between them significant?
Data4aov<-data.frame(Scores,IsFiltered)
aovtest=aov(Scores~IsFiltered, data=Data4aov)       
summary(aovtest)
boxplot(Scores~IsFiltered, data=Data4aov)

# Now check this same filtering effect on suggestion success

FilteredData<-data[,1:12] & data[,15:26] > criterion #Prepare data
PureData<-data[,1:12]
FilteredData<-FilteredData+0
FilteredData<-data.frame(FilteredData)
AllData<-rbind(FilteredData,PureData)
Data4logit<-data.frame(cbind(AllData,IsFiltered))
Data4logit$IsFiltered<-factor(Data4logit$IsFiltered)

hflogit<-glm(hf ~ IsFiltered, data=Data4logit, family = "binomial") 
summary(hflogit)
corrhflogit<-0.0729*12 #And do it: Head falling not s
eclogit<-glm(ec ~ IsFiltered, data=Data4logit, family = "binomial") 
summary(eclogit)
correclogit<-0.078828*12 # Eye closure not s
lhdlogit<-glm(lhd ~ IsFiltered, data=Data4logit, family = "binomial") 
summary(lhdlogit)
corrlhdlogit<-0.111*12 # Hand lowering not s
raflogit<-glm(raf ~ IsFiltered, data=Data4logit, family = "binomial")
summary(raflogit)
corrraflogit<-0.0351*12 # Arm immobilization not s
filogit<-glm(fi ~ IsFiltered, data=Data4logit, family = "binomial") 
summary(filogit)
corrfilogit<-0.00132*12 # Finger locking IS S
larlogit<-glm(lar ~ IsFiltered, data=Data4logit, family = "binomial")
summary(larlogit)
corrlarlogit<-0.0245*12 # # Arm rigidity NOT S
mhlogit<-glm(mh ~ IsFiltered, data=Data4logit, family = "binomial")
summary(mhlogit)
corrmhlogit<-0.0783*12 # Magnetic hands is not S
dhslogit<-glm(dhs ~ IsFiltered, data=Data4logit, family = "binomial")
summary(dhslogit)
corrdhslogit<-0.00575*12 # Head shake IS S
flogit<-glm(f ~ IsFiltered, data=Data4logit, family = "binomial") 
summary(flogit)
corrflogit<- 0.478*12 # Fly is NOT AT ALL S (which is interesting)
eologit<-glm(eo ~ IsFiltered, data=Data4logit, family = "binomial")# Eye opening IS S
summary(eologit)
correologit<- 0.0106*12 # Eye opening not S
lalogit<-glm(la ~ IsFiltered, data=Data4logit, family = "binomial")
summary(lalogit)
corrlalogit<- 4.67e-07*12 # Left Ankle is "WAY TOO S"
amnlogit<-glm(amn ~ IsFiltered, data=Data4logit, family = "binomial")
summary(amnlogit)
corramnlogit<- 0.0346*12 # Amnesia not S

# Item pass rate for filtered and unfiltered

SumPure<-sapply(PureData,sum)
SumFiltered<-sapply(FilteredData,sum)
PercentPure<-sapply(SumPure, function(x) (x/length(PureData$hf))*100)
PercentFiltered<-sapply(SumFiltered, function(x) (x/length(FilteredData$hf))*100)
# Grouped Bar Plot (based on the bar() wrapper I got online)
#bar(dv = Perc, 
 #factors = c(FilterYN,FilterYN), 
 #dataframe = Table4BarPlot, 
 #errbar = FALSE, 
 #ylim=c(0, 140))  #I increased the upper y-limit to accommodate the legend. 
#better to use ggplot2



# Item scale correlation for unfiltered

library(lme4)

hfcorr <- glm(hf ~ I(ec+lhd+raf+ fi+ lar+ mh+ dhs+ f+ eo +la+amn), family='binomial', data=data)
summary(hfcorr)
amncorr <- glm(amn ~ I(hf+ec+lhd+raf+ fi+ lar+ mh+ dhs+ f+ eo +la), family='binomial', data=data)
summary(amncorr)
eccorr <- glm(ec ~ I(hf+lhd+raf+ fi+ lar+ mh+ dhs+ f+ eo +la+amn), family='binomial', data=data)
summary(eccorr)
lhdcorr <- glm(lhd ~ I(hf+ec+raf+ fi+ lar+ mh+ dhs+ f+ eo +la+amn), family='binomial', data=data)
summary(lhdcorr)
rafcorr <- glm(raf ~ I(hf+ec+lhd+ fi+ lar+ mh+ dhs+ f+ eo +la+amn), family='binomial', data=data)
summary(rafcorr)
ficorr <- glm(fi ~ I(hf+ec+lhd+ raf+ lar+ mh+ dhs+ f+ eo +la+amn), family='binomial', data=data)
summary(ficorr)
larcorr <- glm(lar ~ I(hf+ec+lhd+ raf+ fi+ mh+ dhs+ f+ eo +la+amn), family='binomial', data=data)
summary(larcorr)
mhcorr <- glm(mh ~ I(hf+ec+lhd+ raf+ fi+ lar+ dhs+ f+ eo +la+amn), family='binomial', data=data)
summary(mhcorr)
dhscorr <- glm(dhs ~ I(hf+ec+lhd+ raf+ fi+ lar+ mh+ f+ eo +la+amn), family='binomial', data=data)
summary(dhscorr)
fcorr <- glm(f ~ I(hf+ec+lhd+ raf+ fi+ lar+ mh+ dhs+ eo +la+amn), family='binomial', data=data)
summary(fcorr)
eocorr <- glm(eo ~ I(hf+ec+lhd+ raf+ fi+ lar+ mh+ dhs+ f +la+amn), family='binomial', data=data)
summary(eocorr)
lacorr <- glm(la ~ I(hf+ec+lhd+ raf+ fi+ lar+ mh+ dhs+ f +eo+amn), family='binomial', data=data)
summary(lacorr)

# Item scale correlation for filtered


hfcorr <- glm(hf ~ I(ec+lhd+raf+ fi+ lar+ mh+ dhs+ f+ eo +la+amn), family='binomial', data=FilteredData)
summary(hfcorr)
amncorr <- glm(amn ~ I(hf+ec+lhd+raf+ fi+ lar+ mh+ dhs+ f+ eo +la), family='binomial', data=FilteredData)
summary(amncorr)
eccorr <- glm(ec ~ I(hf+lhd+raf+ fi+ lar+ mh+ dhs+ f+ eo +la+amn), family='binomial', data=FilteredData)
summary(eccorr)
lhdcorr <- glm(lhd ~ I(hf+ec+raf+ fi+ lar+ mh+ dhs+ f+ eo +la+amn), family='binomial', data=FilteredData)
summary(lhdcorr)
rafcorr <- glm(raf ~ I(hf+ec+lhd+ fi+ lar+ mh+ dhs+ f+ eo +la+amn), family='binomial', data=FilteredData)
summary(rafcorr)
ficorr <- glm(fi ~ I(hf+ec+lhd+ raf+ lar+ mh+ dhs+ f+ eo +la+amn), family='binomial', data=FilteredData)
summary(ficorr)
larcorr <- glm(lar ~ I(hf+ec+lhd+ raf+ fi+ mh+ dhs+ f+ eo +la+amn), family='binomial', data=FilteredData)
summary(larcorr)
mhcorr <- glm(mh ~ I(hf+ec+lhd+ raf+ fi+ lar+ dhs+ f+ eo +la+amn), family='binomial', data=FilteredData)
summary(mhcorr)
dhscorr <- glm(dhs ~ I(hf+ec+lhd+ raf+ fi+ lar+ mh+ f+ eo +la+amn), family='binomial', data=FilteredData)
summary(dhscorr)
fcorr <- glm(f ~ I(hf+ec+lhd+ raf+ fi+ lar+ mh+ dhs+ eo +la+amn), family='binomial', data=FilteredData)
summary(fcorr)
eocorr <- glm(eo ~ I(hf+ec+lhd+ raf+ fi+ lar+ mh+ dhs+ f +la+amn), family='binomial', data=FilteredData)
summary(eocorr)
lacorr <- glm(la ~ I(hf+ec+lhd+ raf+ fi+ lar+ mh+ dhs+ f +eo+amn), family='binomial', data=FilteredData)
summary(lacorr)

#The alpha (Kuder-Richardson equivallent since it's binary data) for total scale reliability

library(psych)

# Unfiltered

alpha(PureData)

pqvalues<-sapply(PercentPure, function(x) (x/100)*(1-(x/100)))              #by hand
sumpq<-sum(pqvalues)
KR20<-12/11*(1-((sumpq)/(var(PureScores))))

# Filtered

alpha(FilteredData)

pqvalues<-sapply(PercentFiltered, function(x) (x/100)*(1-(x/100)))              #by hand
sumpq<-sum(pqvalues)
KR20<-12/11*(1-((sumpq)/(var(FilteredScores))))

#IT CHECKS OUT BOTH WAYS, QUITE NICELY.

# Rank order correlations for our scale and the reference samples (optional but nice)

