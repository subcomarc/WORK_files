##Analysis Script for the BPCO Pilot Experiment
##By subcomarc, with the resources of Jérôme Sackur
##Data provided by François LARUE and the CHB team
## Needs file "BPCOData.R"

############### ################################
#Preloading settings, directories and libraries#
############### ################################

setwd("/home/subcomarc/Google_Drive/R/")
source("/home/subcomarc/Google_Drive/R/R_Lib_Sackur.r")

library("nlme")
library('xtable')
library('knitr')
library('lme4')
library('car')
library('multcomp')
library('lsmeans')
library("Hmisc")
#library("HH")
library("openxlsx")
library("xlsx")
library("ggplot2")


options(xtable.floating = FALSE)
options(xtable.timestamp = "")


##################### ###################################################
#Load preprocessed data (if this is the first time, DO PREPPING INSTEAD)#
#################### ####################################################

load("/home/subcomarc/Google_Drive/R/BPCOData.R")
load("/home/subcomarc/Google_Drive/R/BPCODeltaData.R")

####################### ########################################
#Prepping the data for stats and adding some convenient columns#
###################### ######################################### ##################

# Load Raw data from Excel file and prepare a "data" dataframe to work with without jeopardizing the original #
# dataRAW <- read.xlsx(file="/home/subcomarc/Google_Drive/R/BPCO.xlsx", 1) # First wave of participants
# dataRAW <- read.xlsx(file="/home/subcomarc/Google_Drive/R/BPCO2.xlsx", 1) # Second wave (6 more)
dataRAW <- read.xlsx(file="/home/subcomarc/Google_Drive/R/BPCO3.xlsx", 1) # Second wave (20 total)
data <- dataRAW

# Rename into informative colnames #

data$patient <- data$Identifiant.numérique.du.patient
data$patient <- as.numeric(data$patient)
data$age <- data$Age
data$date <- data$Date
data$gender <- data$Sexe
data$hypnosis <- data$SÉANCE.D.HYPNOSE..
data$preborg <- data$Échelle.de.Borg.au.début
data$postborg <- data$Échelle.de.Borg.postérieure
data$preFR <- data$Prise.FR.de.repos..après.2.minutes.de.la.demande.de.rester.au.calme.
data$postFR <- data$Prise.FR.postérieure..juste.après.la.procedure.
data$preSO <- data$Prise.SO.de.repos..après.2.minutes.de.la.demande.de.rester.au.calme.
data$postSO <- data$Prise.SO.postérieure..juste.après.la.procedure.
data$preIASTAcalme <- data$Échelle.IASTA.réduite.au.début...Je.me.sens.calme..
data$preIASTAtendu <- data$Échelle.IASTA.réduite.au.début...Je.suis.tendu.e...
data$preIASTAderange <- data$Échelle.IASTA.réduite.au.début...Je.suis.dérangé.e...
data$preIASTAdetendu <- data$Échelle.IASTA.réduite.au.début...Je.suis.détendu.e...
data$preIASTAsatisf <- data$Échelle.IASTA.réduite.au.début...Je.me.sens.satisfait.e...
data$preIASTApreocc <- data$Échelle.IASTA.réduite.au.début...Je.suis.préoccupé.e...
data$postIASTAcalme <- data$Échelle.IASTA.réduite.postérieure...Je.me.sens.calme..
data$postIASTAtendu <- data$Échelle.IASTA.réduite.postérieure...Je.suis.tendu.e...
data$postIASTAderange <- data$Échelle.IASTA.réduite.postérieure...Je.suis.dérangé.e...
data$postIASTAdetendu <- data$Échelle.IASTA.réduite.postérieure...Je.suis.détendu.e...
data$postIASTAsatisf <- data$Échelle.IASTA.réduite.postérieure...Je.me.sens.satisfait.e...
data$postIASTApreocc <- data$Échelle.IASTA.réduite.postérieure...Je.suis.préoccupé.e...


# Reshape dataframe, change inputs and do some more naming #

data$hypnosis <- gsub ("Oui", 1, data$hypnosis)
data$hypnosis <- gsub ("Non \\(c'est l'étape contrôle\\)", 0, data$hypnosis)

for (i in 1:ncol(data)) {
  data[[i]] <- gsub ("Pas du tout", 1, data[[i]])
  data[[i]] <- gsub ("Un peu", 2, data[[i]])
  data[[i]] <- gsub ("Modérément", 3, data[[i]])
  data[[i]] <- gsub ("Beaucoup", 4, data[[i]])
}

for (i in 29:ncol(data)) { #gotta start by the index of the first col to have char vals we need to turn to nums 
  data[[i]] <- as.numeric(data[[i]])
  data[[i]] <- as.numeric(data[[i]])
  data[[i]] <- as.numeric(data[[i]])
  data[[i]] <- as.numeric(data[[i]])
}
#initialize final df
e <- data.frame(gender = character(),
                hypnosis = numeric(),
                borg = numeric(),
                FR = numeric(),
                SO = numeric(),
                IASTAcalme = numeric(),
                IASTAtendu = numeric(),
                IASTAderange = numeric(),
                IASTAdetendu = numeric(),
                IASTAsatisf = numeric(),
                IASTApreocc = numeric(),
                stringsAsFactors=FALSE)

#start by adding pre-values for the hypnosis condition (step 1 of 4)

temp <- evalq(aggregate(list(borg = preborg, FR = preFR, SO = preSO,
                             IASTAcalme = preIASTAcalme, IASTAtendu = preIASTAtendu, IASTAderange = preIASTAderange,
                             IASTAdetendu = preIASTAdetendu, IASTAsatisf = preIASTAsatisf, IASTApreocc = preIASTApreocc), 
                        list(patient = patient), mean), data[data$hypnosis == 1,]) # I can use mean cause there's only one value
temp$hypnosis <- 1 
temp$whenexp <- 1

e <- rbind(e, temp) #put it in the initialized df

#now add post-values for the hypnosis condition (step 2 of 4)

temp <- evalq(aggregate(list(borg = postborg, FR = postFR, SO = postSO,
                             IASTAcalme = postIASTAcalme, IASTAtendu = postIASTAtendu, IASTAderange = postIASTAderange,
                             IASTAdetendu = postIASTAdetendu, IASTAsatisf = postIASTAsatisf, IASTApreocc = postIASTApreocc), 
                        list(patient = patient), mean), data[data$hypnosis == 1,]) # I can use mean cause there's only one value
temp$hypnosis <- 1 
temp$whenexp <- 2

e <- rbind(e, temp) #put it in the initialized df

#now add pre-values for the no hypnosis condition (step 3 of 4)
temp <- evalq(aggregate(list(borg = preborg, FR = preFR, SO = preSO,
                             IASTAcalme = preIASTAcalme, IASTAtendu = preIASTAtendu, IASTAderange = preIASTAderange,
                             IASTAdetendu = preIASTAdetendu, IASTAsatisf = preIASTAsatisf, IASTApreocc = preIASTApreocc), 
                        list(patient = patient), mean), data[data$hypnosis == 0,]) # I can use mean cause there's only one value
temp$hypnosis <- 0 
temp$whenexp <- 1

e <- rbind(e, temp) #put it in the initialized df

#now add post-values for the no hypnosis condition (step 4 of 4)

temp <- evalq(aggregate(list(borg = postborg, FR = postFR, SO = postSO,
                             IASTAcalme = postIASTAcalme, IASTAtendu = postIASTAtendu, IASTAderange = postIASTAderange,
                             IASTAdetendu = postIASTAdetendu, IASTAsatisf = postIASTAsatisf, IASTApreocc = postIASTApreocc), 
                        list(patient = patient), mean), data[data$hypnosis == 0,]) # I can use mean cause there's only one value
temp$hypnosis <- 0 
temp$whenexp <- 2

e <- rbind(e, temp) #put it in the initialized df

# Calculate IASTA single value per participant per condition #
#INSTRUCTIONS: To calculate the total STAI score (range 20 - 80):
# reverse scoring of the positive items (calm, relaxed, content) so 1=4, 2=3, 3=2 and 4=1;
# sum all six scores;
# multiply total score by 20/6;
# refer to Spielberger’s manuals to interpret scores (a ‘normal’ score is approx. 34 - 36)

# reverse scoring of the positive items (calm, relaxed, content) so 1=4, 2=3, 3=2 and 4=1;
for (i in 1:nrow(e)) {
  if (e$IASTAcalme[[i]] == 1){
    e$IASTAcalme[[i]] <- 4} else if (e$IASTAcalme[[i]] == 2){
      e$IASTAcalme[[i]] <- 3} else if (e$IASTAcalme[[i]] == 3){
        e$IASTAcalme[[i]] <- 2} else {
          e$IASTAcalme[[i]] <- 1
          }
}

for (i in 1:nrow(e)) {
  if (e$IASTAdetendu[[i]] == 1){
    e$IASTAdetendu[[i]] <- 4} else if (e$IASTAdetendu[[i]] == 2){
      e$IASTAdetendu[[i]] <- 3} else if (e$IASTAdetendu[[i]] == 3){
        e$IASTAdetendu[[i]] <- 2} else {
          e$IASTAdetendu[[i]] <- 1
        }
}

for (i in 1:nrow(e)) {
  if (e$IASTAsatisf[[i]] == 1){
    e$IASTAsatisf[[i]] <- 4} else if (e$IASTAsatisf[[i]] == 2){
      e$IASTAsatisf[[i]] <- 3} else if (e$IASTAsatisf[[i]] == 3){
        e$IASTAsatisf[[i]] <- 2} else {
          e$IASTAsatisf[[i]] <- 1
        }
}

# sum all six scores and multiply total score by 20/6;
e$IASTAscore <- 1

for (i in 1:nrow(e)){
  e$IASTAscore[[i]] <- (e$IASTAcalme[[i]] + e$IASTAtendu[[i]] + e$IASTAderange[[i]] + e$IASTAdetendu[[i]] +
                          e$IASTAsatisf[[i]] + e$IASTApreocc[[i]]) * (20 / 6)
}
##Quick sneak-peak
hist(e[e$hypnosis == 1 & e$whenexp ==1,]$IASTAscore)
hist(e[e$hypnosis == 1 & e$whenexp ==2,]$IASTAscore)
hist(e[e$hypnosis == 0 & e$whenexp ==1,]$IASTAscore)
hist(e[e$hypnosis == 0 & e$whenexp ==2,]$IASTAscore)
##

# Add one more column: SO/FR. We'll call it respqual

e$respqual <- 1

for (i in 1:nrow(e)){
  e$respqual[[i]] <- (e$SO[[i]] / e$FR[[i]])
}

# Recreate the dataframes only with clean columns

data <- e

# Check for outliers # (1 and 3 for f are the "before" for both conditions)

# f <- evalq(aggregate(list(FR=FR, SO=SO, borg=borg, IASTAscore=IASTAscore), list(whenexp=whenexp, hypnosis=hypnosis), mean), e)
# temp <- evalq(aggregate(list(FR=FR, SO=SO, borg=borg, IASTAscore=IASTAscore), list(whenexp=whenexp, hypnosis=hypnosis), sd), e)
# f$FRsd <- temp$FR
# f$SOsd <- temp$SO
# f$borgsd <- temp$borg
# f$IASTAsd <- temp$IASTAscore

# limFR <- f$FR+3*f$FRsd #check mean RT
# bad <- as.vector(e$patient[e$FR > limFR[1]]) # 102 at the base for hypnosis?
# limSO <- f$SO+3*f$SOsd #check mean RT
# bad <- as.vector(e$patient[e$SO > limSO[1]])
# limborg <- f$borg+3*f$borgsd #check mean RT
# bad <- as.vector(e$patient[e$borg > limborg[1]])
# limIASTA <- f$IASTAscore+3*f$IASTAsd #check mean RT
# bad <- as.vector(e$patient[e$IASTAscore > limIASTA[1]])

# e <- droplevels(e[!e$patient %in% bad,]) # To eliminate subjects based on outliers
#data <- droplevels(data[!data$patient %in% c(102),]) # To eliminate subjects based on outliers
#data <- droplevels(data[!data$patient %in% c(108),]) # To eliminate subjects based on outliers


# Add a dataframe (BPCODeltaData) with 4 new columns: SODelta, FRDelta, borgDelta & IASTADelta

temp1 <- e[e$whenexp ==1,]
temp2 <- e[e$whenexp ==2,]
#temp1 <- data[data$whenexp ==1,] #to change it without having to recreate the whole thing, but don't forget to save it afterwards!
#temp2 <- data[data$whenexp ==2,]

temp1$SODelta <- temp2$SO - temp1$SO
temp1$FRDelta <- temp2$FR - temp1$FR
temp1$borgDelta <- temp2$borg - temp1$borg
temp1$IASTADelta <- temp2$IASTAscore - temp1$IASTAscore

temp1$SODeltaPercent <- ceiling(((temp2$SO - temp1$SO)/temp1$SO)*100)
temp1$FRDeltaPercent <- ceiling(((temp2$FR - temp1$FR)/temp1$FR)*100)
#temp1$borgDeltaPercent <- ceiling(((temp2$borg - temp1$borg)/temp1$borg)*100)
temp1$IASTADeltaPercent <- ceiling(((temp2$IASTAscore - temp1$IASTAscore)/temp1$IASTAscore)*100)

  
deltadata <- evalq(data.frame(patient, SODelta, FRDelta, borgDelta, IASTADelta, 
                              SODeltaPercent, FRDeltaPercent, IASTADeltaPercent, hypnosis), temp1)

# Eliminate abreactions (FR raises dramatically after H) # OPTIONAL, CURRENTLY NOT DOING IT! JUST MARKING THEM

deltadata$FReffective <- 1
deltadata$SOeffective <- 1
deltadata$Borgeffective <- 1
deltadata$IASTAeffective <- 1
#deltadata$effective[deltadata$FRDelta < 0 & deltadata$hypnosis == 1] <- 0
deltadata$FReffective[deltadata$FRDelta > 0] <- 0
deltadata$SOeffective[deltadata$SODelta > 0] <- 0
deltadata$Borgeffective[deltadata$borgDelta > 0] <- 0
deltadata$IASTAeffective[deltadata$IASTADelta > 0] <- 0

deltadata$FReffective[deltadata$FRDelta == 0] <- 2
deltadata$SOeffective[deltadata$SODelta == 0] <- 2
deltadata$Borgeffective[deltadata$borgDelta == 0] <- 2
deltadata$IASTAeffective[deltadata$IASTADelta == 0] <- 2

deltadata$FReffective <- factor(deltadata$FReffective, labels = c("Up", "Down", "Same"))
deltadata$SOeffective <- factor(deltadata$SOeffective, labels = c("Up", "Down", "Same"))
deltadata$Borgeffective <- factor(deltadata$Borgeffective, labels = c("Up", "Down", "Same"))
deltadata$IASTAeffective <- factor(deltadata$IASTAeffective, labels = c("Up", "Down", "Same"))

data$FReffective <- 1
data$SOeffective <- 1
data$Borgeffective <- 1
data$IASTAeffective <- 1

data[data$hypnosis==1 & data$whenexp==1,]$FReffective <- deltadata[deltadata$hypnosis==1,]$FReffective
data[data$hypnosis==1 & data$whenexp==2,]$FReffective <- deltadata[deltadata$hypnosis==1,]$FReffective
data[data$hypnosis==0 & data$whenexp==1,]$FReffective <- deltadata[deltadata$hypnosis==0,]$FReffective
data[data$hypnosis==0 & data$whenexp==2,]$FReffective <- deltadata[deltadata$hypnosis==0,]$FReffective

data[data$hypnosis==1 & data$whenexp==1,]$SOeffective <- deltadata[deltadata$hypnosis==1,]$SOeffective
data[data$hypnosis==1 & data$whenexp==2,]$SOeffective <- deltadata[deltadata$hypnosis==1,]$SOeffective
data[data$hypnosis==0 & data$whenexp==1,]$SOeffective <- deltadata[deltadata$hypnosis==0,]$SOeffective
data[data$hypnosis==0 & data$whenexp==2,]$SOeffective <- deltadata[deltadata$hypnosis==0,]$SOeffective

data[data$hypnosis==1 & data$whenexp==1,]$Borgeffective <- deltadata[deltadata$hypnosis==1,]$Borgeffective
data[data$hypnosis==1 & data$whenexp==2,]$Borgeffective <- deltadata[deltadata$hypnosis==1,]$Borgeffective
data[data$hypnosis==0 & data$whenexp==1,]$Borgeffective <- deltadata[deltadata$hypnosis==0,]$Borgeffective
data[data$hypnosis==0 & data$whenexp==2,]$Borgeffective <- deltadata[deltadata$hypnosis==0,]$Borgeffective

data[data$hypnosis==1 & data$whenexp==1,]$IASTAeffective <- deltadata[deltadata$hypnosis==1,]$IASTAeffective
data[data$hypnosis==1 & data$whenexp==2,]$IASTAeffective <- deltadata[deltadata$hypnosis==1,]$IASTAeffective
data[data$hypnosis==0 & data$whenexp==1,]$IASTAeffective <- deltadata[deltadata$hypnosis==0,]$IASTAeffective
data[data$hypnosis==0 & data$whenexp==2,]$IASTAeffective <- deltadata[deltadata$hypnosis==0,]$IASTAeffective

data$FReffective <- factor(data$FReffective, labels = c("Up", "Down", "Same"))
data$SOeffective <- factor(data$SOeffective, labels = c("Up", "Down", "Same"))
data$Borgeffective <- factor(data$Borgeffective, labels = c("Up", "Down", "Same"))
data$IASTAeffective <- factor(data$IASTAeffective, labels = c("Up", "Down", "Same"))


#contrasts(data$effective)=contr.sum(3)


#bad <- as.vector(deltadata$patient[deltadata$FRDelta < 0 & deltadata$hypnosis == 1])
 
#data <- droplevels(data[!data$patient %in% bad,]) # To eliminate subjects based on outliers
data <- droplevels(data[!data$patient %in% c("101"),]) # To eliminate subjects based on outliers
deltadata <- droplevels(deltadata[!deltadata$patient %in% c("101"),]) # To eliminate subjects based on outliers

#Our data is squeaky clean, we are ready to save!
save(data, file='BPCOData.R')
save(deltadata, file='BPCODeltaData.R')
#write.xlsx(data, "/home/subcomarc/Google_Drive/R/BPCOexcel_Clean.xlsx") #For sharing clean data with medical doctors, who need an excel file

###################### ##########################
#Quick stats - descriptive and some interactions#
##################### ###########################

#FR
evalq(plot(patient, FR, col = "red", type = "p"), data[data$hypnosis %in% c(1),]) 
lines(data$patient[data$hypnosis %in% c(0)], data$FR[data$hypnosis %in% c(0)], col = "blue", type = "p") #No big diffs
l <- lm(FR ~ patient, data[data$hypnosis == 1,])
abline(l, col = "red")
l <- lm(FR ~ patient, data[data$hypnosis == 0,])
abline(l, col = "blue")

#SO
evalq(plot(patient, SO, col = "red", type = "p"), data[data$hypnosis %in% c(1),]) 
lines(data$patient[data$hypnosis %in% c(0)], data$SO[data$hypnosis %in% c(0)], col = "blue", type = "p") #No big diffs
l <- lm(SO ~ patient, data[data$hypnosis == 1,])
abline(l, col = "red")
l <- lm(SO ~ patient, data[data$hypnosis == 0,])
abline(l, col = "blue")


#respqual
evalq(plot(patient, respqual, col = "red", type = "p"), data[data$hypnosis %in% c(1),]) 
lines(data$patient[data$hypnosis %in% c(0)], data$respqual[data$hypnosis %in% c(0)], col = "blue", type = "p") #No big diffs
l <- lm(respqual ~ patient, data[data$hypnosis == 1,])
abline(l, col = "red")
l <- lm(respqual ~ patient, data[data$hypnosis == 0,])
abline(l, col = "blue")

#borg
evalq(plot(patient, borg, col = "red", type = "p"), data[data$hypnosis %in% c(1),]) 
lines(data$patient[data$hypnosis %in% c(0)], data$borg[data$hypnosis %in% c(0)], col = "blue", type = "p") #No big diffs
l <- lm(borg ~ patient, data[data$hypnosis == 1,])
abline(l, col = "red")
l <- lm(borg ~ patient, data[data$hypnosis == 0,])
abline(l, col = "blue")

#IASTA
evalq(plot(patient, IASTAscore, col = "red", type = "p"), data[data$hypnosis %in% c(1),]) 
lines(data$patient[data$hypnosis %in% c(0)], data$IASTAscore[data$hypnosis %in% c(0)], col = "blue", type = "p") #No big diffs
l <- lm(IASTAscore ~ patient, data[data$hypnosis == 1,])
abline(l, col = "red")
l <- lm(IASTAscore ~ patient, data[data$hypnosis == 0,])
abline(l, col = "blue")

###################### #########################
#All contrasts and dataframes for linear models#
##################### ##########################
datatemp<-data
datatemp[datatemp$hypnosis==1,]$hypnosis<-"Hypnosis"  #for plots
datatemp[datatemp$hypnosis==0,]$hypnosis<-"Control"  #for plots
datatemp[datatemp$whenexp==1,]$whenexp<-"Before"  #for plots
datatemp[datatemp$whenexp==2,]$whenexp<-"After"  #for plots

deltadatatemp<-deltadata
deltadatatemp[deltadatatemp$hypnosis==1,]$hypnosis<-"Hypnosis"  #for plots
deltadatatemp[deltadatatemp$hypnosis==0,]$hypnosis<-"Control"  #for plots

data$hypnosis <- factor(data$hypnosis, labels = c("Control","Hypnosis"))
data$whenexp <- factor(data$whenexp, labels = c("Before","After"))
datatemp$hypnosis <- factor(data$hypnosis, labels = c("Control", "Hypnosis")) #for plots
datatemp$whenexp <- factor(data$whenexp, labels = c("Before","After")) #for plots

contrasts(data$hypnosis)=contr.sum(2)
contrasts(data$whenexp)=contr.sum(2)
contrasts(datatemp$hypnosis)=contr.sum(2)
contrasts(datatemp$whenexp)=contr.sum(2)

deltadata$hypnosis <- factor(deltadata$hypnosis, labels = c("Control","Hypnosis"))
deltadatatemp$hypnosis <- factor(deltadatatemp$hypnosis,labels = c("Control","Hypnosis")) #for plots

contrasts(deltadata$hypnosis)=contr.sum(2)
contrasts(deltadatatemp$hypnosis)=contr.sum(2)


control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))

###################### ###################### ######################
#STATISTICAL ANALYSES#
###################### ###################### ######################

###HYPNOSIS ACROSS INSTANCE (FR)

#Can we apply linear models?? Let's check the conditions (Linearity, absence of colinearity, 
# homoskedacity, Normality of residuals, Independence)

HYPS <- lmer(FR ~ hypnosis*whenexp + (1|patient), data)

HYPS <- lmer(FR ~ whenexp + (1|patient), droplevels(data[data$hypnosis %in% c("Hypnosis")
                                                         ,]))
HYPSnull <- lmer(SO ~ 1 + (1|patient), droplevels(data[data$hypnosis %in% c("Control")
                                                       ,]))

HYPS <- lmer(FRDeltaPercent ~ hypnosis + (1|patient), deltadata
                                                         )

Anova(HYPS)
anova(HYPS,HYPSnull)
exp((152.19 - 152.29)/2) #exp((BICnull - BICfull)/2)

HYPS <- lmer(borgDelta ~ hypnosis + (1|patient), droplevels(deltadata[deltadata$patient %in% 
                                                                      deltadata[deltadata$FReffective %in% c("Down","Up","Same")
                                                                    & deltadata$hypnosis %in% c("Hypnosis"),]$patient
                                                          ,]))

HYPSnull <- lmer(borgDelta ~ 1 + (1|patient), droplevels(deltadata[deltadata$patient %in% 
                                                                      deltadata[deltadata$FReffective %in% c("Down", "Up", "Same")
                                                                                & deltadata$hypnosis %in% c("Hypnosis"),]$patient
                                                                    ,]))


Anova(HYPS)

anova(HYPS,HYPSnull)
exp((112.3 - 115.7)/2) #exp((BICnull - BICfull)/2)

hist(residuals(HYPS)) # Residuals are normally distributed= "Tukey")))
qqnorm(residuals(HYPS)) # The quantiles progress sort of linearly too
plot(fitted(HYPS),residuals(HYPS)) ##

Anova(HYPS)
lsm.options(pbkrtest.limit = 9557)
lsmeans(HYPS, pairwise ~ hypnosis*whenexp, adjust="tukey")
lsmeans(HYPS, pairwise ~ hypnosis, adjust="tukey") #for diff
lsmeans(HYPS, pairwise ~ whenexp, adjust="tukey") #for diff


###HYPNOSIS ACROSS INSTANCE (SO)

#Can we apply linear models?? Let's check the conditions (Linearity, absence of colinearity, 
# homoskedacity, Normality of residuals, Independence)

HYPS <- lmer(SO ~ hypnosis*whenexp + (1|patient), data)
HYPS <- lmer(SODelta ~ hypnosis + (1|patient), deltadata) #for the diff

hist(residuals(HYPS)) # Residuals are normally distributed= "Tukey")))

qqnorm(residuals(HYPS)) # The quantiles progress sort of linearly too
plot(fitted(HYPS),residuals(HYPS)) ##There's some trouble with some of the answers... my money is on subject 14!

Anova(HYPS)
lsm.options(pbkrtest.limit = 9557)
lsmeans(HYPS, pairwise ~ hypnosis*whenexp, adjust="tukey")
lsmeans(HYPS, pairwise ~ hypnosis, adjust="tukey") #for diff

###HYPNOSIS ACROSS INSTANCE (respqual)

#Can we apply linear models?? Let's check the conditions (Linearity, absence of colinearity, 
# homoskedacity, Normality of residuals, Independence)

HYPS <- lmer(respqual ~ hypnosis*whenexp + (1|patient), data)

hist(residuals(HYPS)) # Residuals are normally distributed= "Tukey")))

qqnorm(residuals(HYPS)) # The quantiles progress sort of linearly too
plot(fitted(HYPS),residuals(HYPS)) ##There's some trouble with some of the answers... my money is on subject 14!

Anova(HYPS)
lsm.options(pbkrtest.limit = 9557)
lsmeans(HYPS, pairwise ~ hypnosis*whenexp, adjust="tukey")

###HYPNOSIS ACROSS INSTANCE (borg)

#Can we apply linear models?? Let's check the conditions (Linearity, absence of colinearity, 
# homoskedacity, Normality of residuals, Independence)

HYPS <- lmer(borg ~ hypnosis*whenexp + (1|patient), data)
HYPS <- lmer(borgDelta ~ hypnosis + (1|patient), deltadata) #for the diff


hist(residuals(HYPS)) # Residuals are normally distributed= "Tukey")))

qqnorm(residuals(HYPS)) # The quantiles progress sort of linearly too
plot(fitted(HYPS),residuals(HYPS)) ##There's some trouble with some of the answers... my money is on subject 14!

Anova(HYPS)
lsm.options(pbkrtest.limit = 9557)
lsmeans(HYPS, pairwise ~ hypnosis*whenexp, adjust="tukey")
lsmeans(HYPS, pairwise ~ hypnosis, adjust="tukey") #for diff

###HYPNOSIS ACROSS INSTANCE (IASTAscore)

#Can we apply linear models?? Let's check the conditions (Linearity, absence of colinearity, 
# homoskedacity, Normality of residuals, Independence)

HYPS <- lmer(IASTAscore ~ hypnosis*whenexp + (1|patient), data)
HYPS <- lmer(IASTADelta ~ hypnosis + (1|patient), deltadata) #for the diff

hist(residuals(HYPS)) # Residuals are normally distributed= "Tukey")))

qqnorm(residuals(HYPS)) # The quantiles progress sort of linearly too
plot(fitted(HYPS),residuals(HYPS)) ##There's some trouble with some of the answers... my money is on subject 14!

Anova(HYPS)
lsm.options(pbkrtest.limit = 9557)
lsmeans(HYPS, pairwise ~ hypnosis*whenexp, adjust="tukey")


###################### ###################### ######################
#NICER PLOTS#
###################### ###################### ######################

#ggplot(data[data$hypnosis==1,], aes(x=whenexp, y=IASTAscore, fill=whenexp)) + #Differences in IASTA scores before and after hyp
#ggplot(data[data$hypnosis==0,], aes(x=whenexp, y=IASTAscore, fill=whenexp)) + #Differences in IASTA scores before and after ctrl
#ggplot(data[data$hypnosis==1,], aes(x=whenexp, y=borg, fill=whenexp)) +  #Differences in borg scores before and after hypnosis
#ggplot(data[data$hypnosis==0,], aes(x=whenexp, y=borg, fill=whenexp)) +  #Differences in borg scores before and after control
#ggplot(data[data$whenexp==2,], aes(x=hypnosis, y=borg, fill=hypnosis)) + #Differences in borg scores between hypnosis and control (after)
#ggplot(data, aes(x=hypnosis, y=borg, fill=whenexp)) +  #Differences in borg scores bef ctrl vs. aft hyp. With geom_col
#ggplot(data[data$hypnosis==1,], aes(x=whenexp, y=FR, fill=whenexp)) + #Differences in respqual scores before and after hypnosis
#ggplot(data, aes(x=hypnosis, y=respqual, fill=whenexp)) +  #Differences in respqual scores bef ctrl vs. aft hyp. With geom_col
#ggplot(data[data$hypnosis==0,], aes(x=whenexp, y=respqual, fill=whenexp)) + #Differences in respqual scores before and after control
#ggplot(data[data$whenexp==2,], aes(x=hypnosis, y=respqual, fill=hypnosis)) + #Differences in respqual scores between hypnosis and control (after)
#ggplot(data[data$hypnosis==1,], aes(x=whenexp, y=FR, fill=hypnosis)) + #Differences in SO and FR scores separatedly before and after hypnosis
#ggplot(data, aes(x=whenexp, y=borg, fill=whenexp)) + #Differences in SO before and after intervention
#ggplot(data[data$whenexp==2,], aes(x=hypnosis, y=FR, fill=hypnosis)) + #Differences in SO and FR scores separatedly before and after hypnosis
#ggplot(deltadata[deltadata$effective == 1,], aes(x=hypnosis, y=FRDelta, fill=hypnosis)) + #Differences in SO and FR with delta and effective
ggplot(deltadata, aes(x=hypnosis, y=borg, fill=hypnosis)) + #Differences in SO and FR with delta and effective
#ggplot(data[data$hypnosis=="Control",], aes(x=whenexp, y=IASTAscore)) + #Differences in SO and FR with delta and effective  
  # scale_shape_discrete(solid=T, legend=F)
  #geom_point(aes(size=3, color=IASTAeffective)) +
  #geom_line(aes(group=patient, x=whenexp, y=IASTAscore, color=IASTAeffective)) +
  #geom_line(size=1) +
  geom_boxplot(size=1) +
  #geom_text() + #remember to add label variable in then aes on top if you use this
  # annotate("text", label=mean(deltadata[deltadata$hypnosis == 1,]$SODelta), x=2, y=86, size=8, colour="mediumslateblue") +
  # annotate("text", label=mean(deltadata[deltadata$hypnosis == 0,]$SODelta), x=1, y=86, size=8, colour="lightblue1") +
  #annotate("text", label=ceiling(mean(deltadata[deltadata$hypnosis == 0,]$IASTADelta)), x=1, y=-12, size=8, colour="lightblue1") +
  #annotate("text", label=ceiling(mean(deltadata[deltadata$hypnosis == 1,]$IASTADelta)), x=2, y=-3, size=8, colour="mediumslateblue") +
  # geom_col(position = "dodge") +
  #scale_fill_discrete(labels=c("Control", "Hypnosis")) +
  scale_fill_manual(labels=c("Control", "Hypnosis"), values = c("lightblue1","mediumslateblue")) +
  #scale_color_manual(labels=c("Control", "Hypnosis"), values = c("red","lightblue1", "black")) +
  #geom_text(size=14, color="black", aes(x=1.5, y=26, label="***")) +
  #geom_segment(size=2, color="black", aes(x=1,xend=2,y=25,yend=25)) +
  #geom_segment(size=2, color="black", aes(x=0.75,xend=0.75,y=9,yend=15.1)) + #for barplot
  #geom_segment(size=2, color="black", aes(x=2.25,xend=2.25,y=6,yend=15.1)) + #for barplot
  #stat_summary(geom = "errorbar", position = "dodge",  ymin=c(7)-2, ymax=7+2) + #for barplot
  labs(title = "Change in Oxygen Saturation score\n", 
       x = "\n", 
       y = "% Change") +
  # geom_line() +
  # scale_fill_manual(values = c("red", "blue")) +
  #ggtitle("Differences in respqual scores before and after hypnosis") +
  #theme(axis.text.x = element_text(angle=0, face="bold", colour="black")) +
  theme_bw() +
  theme(legend.text = element_text(angle=0, face="bold", colour="black", size="14"),
        legend.title = element_text(angle=0, face="bold", colour="white"),
        panel.grid.major = element_blank(), # panel.grid.minor = element_blank(), 
        axis.line=element_line(size=1, color="black"), panel.border=element_blank(), 
        axis.text.x=element_text(size=14), axis.title.x=element_text(size=16),
        axis.text.y=element_text(size=14), axis.title.y=element_text(size=16),
        plot.title=element_text(size=20, face="bold", color="black"))
        #xlab(NULL)
 

#+ stat_smooth(method='lm', se=FALSE)
#geom_abline
#dev.off()
#qplot(condition, RT , data=e)
#lines(e$condition, e$RT)
