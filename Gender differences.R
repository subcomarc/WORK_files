load('C:/Users/lscpuser/Downloads/HarvardData.R')
Gender <-read.table(file='C:/Users/lscpuser/Google Drive/R/HarvardGender/Gender.txt', header = T) # 1 female, 0 male
data<-cbind(data,Gender)

# Gender differences (bootstrapping)
criterion <- 2
gender<-data.frame(data[,47]) # 1 female, 0 male (this is the gender info that I attached)
filtscores<-data.frame(apply(data.frame((data[,1:12] & data[,15:26] > criterion)+0), 1, sum))
FilteredData<-cbind(filtscores, gender) #Prepare data
colnames(FilteredData)<-c("Hypnotizability","Gender")
FilteredData$Gender<-factor(FilteredData$Gender)
scores<-data.frame(apply(data[,1:12], 1, sum))
gender<-data.frame(data[,47])
PureData<-cbind(scores, gender)
colnames(PureData)<-c("Hypnotizability","Gender")
PureData$Gender<-factor(PureData$Gender)

Malegrouplength<-length(PureData[PureData$Gender=="0",]$Hypnotizability)
Femalegrouplength<-length(PureData[PureData$Gender=="1",]$Hypnotizability)
MeanMale<-mean(PureData[PureData$Gender=="0",]$Hypnotizability)
MeanFemale<-mean(PureData[PureData$Gender=="1",]$Hypnotizability)
Difference<-MeanFemale-MeanMale
B = 1000
#Uncorrected
MaleSamples<-with(PureData, matrix(sample(Hypnotizability[Gender=="0"], size = Malegrouplength * B, 
                                          replace = TRUE), B, Malegrouplength))

FemaleSamples<-with(PureData, matrix(sample(Hypnotizability[Gender=="1"], size = Femalegrouplength * B,
                                            replace = TRUE), B, Femalegrouplength))

AllFMeans = apply(FemaleSamples, 1, mean)
AllMMeans = apply(MaleSamples, 1, mean)
boot.stat =  AllFMeans - AllMMeans

checkwOMean <- lapply(boot.stat, function(x) x>Difference)
checkwOMean<-data.frame(checkwOMean)
checkwOMean<-checkwOMean+0
sum(checkwOMean)

ggplot(data.frame(x = boot.stat), aes(x = x)) + geom_density()

xbars = with(PureData, by(Hypnotizability, Gender, mean))
me = 2 * sd(boot.stat)
(xbars[2] - xbars[1]) + c(-1, 1) * me

round((xbars[2] - xbars[1]) + c(-1, 1) * me, 1)

  #Significance (Permutation test)

pool<-PureData
pool<-FilteredData

AllMeansDiff<-NULL

for (i in 1 : B) {
  
  resample <- sample (c(1:nrow (pool)), nrow(pool))
  
  FemalePerm = pool[resample,1][1 : (Femalegrouplength+1)]
  MalePerm = pool[resample,1][(Femalegrouplength+2) : nrow(pool)]
  AllMeansDiff[i] = mean (FemalePerm) - mean (MalePerm) 
  
}

pvalue <- (sum (abs (AllMeansDiff) >= abs(Difference)) + 1)/ (B+1)

#Corrected

MaleSamples<-with(FilteredData, matrix(sample(Hypnotizability[Gender=="0"], size = Malegrouplength * B, 
                                              replace = TRUE), B, Malegrouplength))

FemaleSamples<-with(FilteredData, matrix(sample(Hypnotizability[Gender=="1"], size = Femalegrouplength * B,
                                                replace = TRUE), B, Femalegrouplength))

AllFMeans = apply(FemaleSamples, 1, mean)
AllMMeans = apply(MaleSamples, 1, mean)
boot.stat = AllMMeans - AllFMeans

checkwOMean <- lapply(boot.stat, function(x) x>Difference)
checkwOMean<-data.frame(checkwOMean)
checkwOMean<-checkwOMean+0
sum(checkwOMean)

ggplot(data.frame(x = boot.stat), aes(x = x)) + geom_density()

xbars = with(PureData, by(Hypnotizability, Gender, mean))
me = 2 * sd(boot.stat)
(xbars[2] - xbars[1]) + c(-1, 1) * me

round((xbars[2] - xbars[1]) + c(-1, 1) * me, 1)

  #Significance
