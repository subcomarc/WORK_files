##Analysis Script for the Stroop Pilot Experiment
##By subcomarc, with the resources of Jérôme Sackur
##Data provided by Pauline Demory and Aurore Lemonnier
## Needs file "StroopData.R"

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
#options(xtable.timestamp = "")


##################### ###################################################
#Load preprocessed data (if this is the first time, DO PREPPING INSTEAD)#
#################### ####################################################

load("/home/subcomarc/Google_Drive/R/StroopData.R")
load("/home/subcomarc/Google_Drive/R/RawStroopData.R")

####################### ########################################
#Prepping the data for stats and adding some convenient columns#
###################### ######################################### ##################

# Load Raw data from Excel file and prepare a "data" dataframe to work with without jeopardizing the original #
dataRAW <- read.xlsx(file="/home/subcomarc/Google_Drive/R/Pps1_64_Hernan.xlsx", 1)
#dataRAW_2 <- read.xlsx(file="/home/subcomarc/Google_Drive/R/Pps1_20_Hernan.xlsx", 1)
data <- dataRAW
#data <- dataRAW_2

# Change number of color for color name  #
for (i in 1:ncol(data)) {
  data[[i]] <- gsub ("88,41,0", "marron", data[[i]])
  data[[i]] <- gsub ("205,0,0", "rouge", data[[i]])
  data[[i]] <- gsub ("237,137,10", "orange", data[[i]])
  data[[i]] <- gsub ("0,205,0", "vert", data[[i]])
}

# Change NULLs for NAs #
for (i in 1:ncol(data)) {
  data[[i]] <- gsub ("NULL", NA, data[[i]])
}

## Change NAs for 9999 # (should you want to)
#datamatrix <- as.matrix(data)
#datamatrix[is.na(data) == TRUE] <- 9999
#data <- as.data.frame(datamatrix)

# Rename columns into intelligible info (and get rid of NAs when convenient) #
data$subject <- data$Subject
data$session <- data$Session
data$age <- data$Age
data$handedness <- data$Handedness
data$gender <- data$Sex
data$language <- data$langue.RESP
data$language[data$language == 2] <- 0 
data$condition <- data$condition.Block. 
data$condition[is.na(data$condition == TRUE)] <- data$condition.LogLevel5.[is.na(data$condition == TRUE)]
data$training <- 1 ## if at some point we wanna contrast training vs experiment
data$training[is.na(data$condition.Block. == TRUE)] <- 0 ## if at some point we wanna contrast training vs experiment
data$color <- data$color.Block.
data$color[is.na(data$color == TRUE)] <- data$color.LogLevel5.[is.na(data$color == TRUE)]
data$target <- data$target.Block.
data$target[is.na(data$target == TRUE)] <- data$target.LogLevel5.[is.na(data$target == TRUE)]
data$RT <- data$Target.RT.Block.
data$RT[is.na(data$RT == TRUE)] <- data$Target.RT.LogLevel5.[is.na(data$RT == TRUE)]
data$RT <- as.numeric(data$RT)
data$list <- data$Running.LogLevel5.
data$list[is.na(data$list == TRUE)] <- 0
data$accuracy <- data$Accuracy
data$accuracy[is.na(data$accuracy == TRUE)] <- 2
data$congruence[is.na(data$congruence == TRUE)] <- 2
data$congruence <- gsub("Congruent", 1, data$congruence) #we'll turn the categories into numbers for ease of analysis
data$congruence <- gsub("Incongruent", 0, data$congruence) #we'll turn the categories into numbers for ease of analysis
data$congruency <- data$congruence
data$condition <- gsub("Training", 0, data$condition) #we'll turn the categories into numbers for ease of analysis
data$condition <- gsub("Neutral", 1, data$condition) #we'll turn the categories into numbers for ease of analysis
data$condition <- gsub("Emotionnal", 2, data$condition) #we'll turn the categories into numbers for ease of analysis
data$condition <- gsub("Semantic", 3, data$condition) #we'll turn the categories into numbers for ease of analysis
data$condition <- gsub("Color", 4, data$condition) #we'll turn the categories into numbers for ease of analysis


# Recreate the dataframe only with clean columns, create one with the training separatedly #
# data = subject session age handedness gender language condition color target  RT list accuracy congruency #
# For condition: 0-Training 1-Neutral 2-Emotionnal 3-Semantic 4-Color#
# For congruency: 1-Congruent 0-Incongruent 2-Nonsensical#
allcleandata <- as.data.frame(data[21:34])
trainingdata <- allcleandata[allcleandata$training == 1,]
trainingdata <- trainingdata[,-8]
data <- allcleandata[allcleandata$training == 0,]
data <- data[,-8]

# Let's have a look at the training sessions #

hist(trainingdata$RT)
hist(data$RT)
plot(density(data$RT)) #in black
lines(density(trainingdata$RT), col = 2) #in red


# Check for outliers #

unique(data$accuracy) 
table(data$accuracy, data$subject) # it appears something is wrong with subject number 5, I'll kill him - what's going on??
data <- droplevels(data[!data$subject %in% c("5", "21", "32", "37"),]) # To eliminate subjects based on experiment behavior

unique(data$RT) 
plot(data$subject, data$RT) # subjects 2 and 14 have 0 RT on 1 trial each - what's going on??

data$accuracy <- as.numeric(data$accuracy)# make sure they are numbers!
data$RT <- as.numeric(data$RT)
data$condition <- as.numeric(data$condition)
data$subject <- as.numeric(data$subject)

data <- data[data$RT > 40,]##less that 40 ms and more than 3 seconds are OUT
data <- data[data$RT < 3000,]##less that 40 ms and more than 3 seconds are OUT
#data <- data[data$accuracy == 1]##wrong answers are OUT!

e <- evalq(aggregate(list(accuracy=accuracy, RT=RT), list(subject=subject, condition=condition), mean), data) #build a subset to have a look
evalq(plot(subject, accuracy, col=condition, type='p'), e)
evalq(plot(subject, RT, col=condition, type='p'), e)
l <- lm(RT~subject, e)
abline(l, col = "red") 


hist(e$RT); mean(e$RT); sd(e$RT)

limRT <- mean(e$RT)+3*sd(e$RT) #check mean RT
bad <- as.vector(e$subject[e$RT > limRT]) 
data <- droplevels(data[!data$subject %in% bad,]) # To eliminate subjects based on outliers

#build a subset to have a look with median
e <- evalq(aggregate(list(RT=RT), list(subject=subject, condition=condition), median), data) 
limRT <- median(e$RT)+3*sd(e$RT) #with mean
bad <- as.vector(e$subject[e$RT > limRT]) 
data <- droplevels(data[!data$subject %in% bad,]) # To eliminate subjects based on outliers


#Our data is squeaky clean, we are ready to save!
save(data, file='StroopData.R')
save(dataRAW, file='RawStroopData.R')

###################### ##########################
#Quick stats - descriptive and some interactions#
##################### ###########################

e <- evalq(aggregate(list(accuracy=accuracy, RT=RT), list(subject=subject, congruency=congruency), mean), 
           data[data$condition %in% c(3, 4),])#Any terrible differences related to congruency?

evalq(plot(subject, RT, col = "red", type = "p"), e[e$congruency == 0,]) 
lines(e$subject[e$congruency == 1], e$RT[e$congruency == 1], col = "blue", type = "p") #No big diffs

e <- evalq(aggregate(list(RT=RT), list(subject=subject, congruency=congruency, condition=condition), mean), 
           droplevels(data[data$condition %in% c(3, 4) & data$congruency %in% c(0, 1),]))
evalq(interaction.plot(condition, congruency, RT), e) ## No semantic stroop? Clear color stroop

e <- evalq(aggregate(list(RT=RT), list(subject=subject, congruency=congruency, condition=condition), mean), 
           droplevels(data[data$condition %in% c(2, 4) & data$congruency %in% c(0, 2),]))
evalq(plot(condition, RT), e) ## Emotional vs color stroop


e <- evalq(aggregate(list(RT=RT), list(subject=subject, condition=condition, congruency=congruency), mean), 
           droplevels(data[data$condition %in% c(1, 3) & data$congruency %in% c(1, 2),]))
evalq(plot(condition, RT), e) ## No  semantic stroop compared to normal?

e <- evalq(aggregate(list(RT=RT), list(subject=subject, condition=condition), mean), 
           droplevels(data[data$condition %in% c(1, 2) & data$congruency %in% c(2),]))
evalq(plot(condition, RT), e) ## Very slight emotional stroop? Emotional Vs Normal


e <- evalq(aggregate(list(RT=RT), list(subject=subject, condition=condition), mean), 
           droplevels(data[data$condition %in% c(2, 3) & data$congruency %in% c(1, 2),]))
evalq(plot(condition, RT), e) ## Emotional Vs Semantic



###################### #########################
#All contrasts and dataframes for linear models#
##################### ##########################
data$condition <- factor(data$condition, labels = c("Neutral","Emotional","Semantic","Color"))
data$congruency <- factor(data$congruency)
data$color <- factor(data$color)
data$gender <- factor(data$gender)
data$list <- factor(data$list) #for checking if there's a list effect in the semantic condition

contrasts(data$condition)=contr.sum(4)
contrasts(data$congruency)=contr.sum(3)
contrasts(data$color)=contr.sum(6)
contrasts(data$gender)=contr.sum(2)
contrasts(data$list)=contr.sum(8) #for checking if there's a list effect in the semantic condition

control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000))

###################### ###################### ######################
#STATISTICAL ANALYSES#
###################### ###################### ######################
# data = subject session age handedness gender language condition color target  RT list accuracy congruency #
# For condition: 0-Training 1-Neutral 2-Emotionnal 3-Semantic 4-Color#
# For congruency: 1-Congruent 0-Incongruent 2-Nonsensical#

###SEMANTIC VS COLOR STROOP

#Can we apply linear models?? Let's check the conditions (Linearity, absence of colinearity, 
# homoskedacity, Normality of residuals, Independence)

SemVsCol <- lmer(RT ~ condition*congruency + (1|subject), data=droplevels(data[data$condition %in% c(3, 4),]))
SemVsCol2 <- lmer(RT ~ condition+congruency + (1|subject), data=droplevels(data[data$condition %in% c(3, 4),]))
SemVsCol3 <- lmer(RT ~ condition + (1|subject), data=droplevels(data[data$condition %in% c(3, 4),]))

#SemVsCol <- lmer(RT ~ congruency + (1|subject) + (1|condition), data=droplevels(data[data$condition %in% c(3, 4),]))
#SemVsCol <- lmer(RT ~ congruency * condition + (1+condition|subject), data=droplevels(data[data$condition %in% c(3, 4),]))

hist(residuals(SemVsCol)) # Residuals are normally distributed= "Tukey")))

qqnorm(residuals(SemVsCol)) # The quantiles progress sort of linearly too
plot(fitted(SemVsCol),residuals(SemVsCol)) ##There's some trouble with some of the answers... my money is on subject 14!


SemVsCol <- lmer(RT ~ condition * congruency + (1|subject), data=droplevels(data[data$condition %in% c(3, 4),]))
anova(SemVsCol,SemVsCol2,SemVsCol3)
Anova(SemVsCol)
lsm.options(pbkrtest.limit = 9557)
lsmeans(SemVsCol, pairwise ~ interaction(condition,congruency), adjust="tukey")
#lsmeans(SemVsCol, pairwise ~ congruency, adjust="tukey")
#lsmeans(SemVsCol, pairwise ~ congruency*condition, adjust="tukey")


datasubs <- droplevels(data[data$condition %in% c(3, 4),])# redo with subscripts and multicomp
datasubs$dataint1 <- interaction(datasubs$condition, datasubs$congruency) # redo with subscripts and multicomp
SemVsCol <- lmer(RT ~ dataint1 + (1|subject), data=datasubs) # # redo with subscripts and multicomp
#SemVsCol <- lmer(RT ~ dataint1 + (1+condition|subject), data=datasubs) # # redo with subscripts and multicomp

Anova(SemVsCol)
summary(glht(SemVsCol, linfct=mcp(dataint1="Tukey")))


###SEMANTIC ALONE

datasubs <- droplevels(data[data$condition %in% c("Semantic"),])
SemConVsInc <- lmer(RT ~ congruency + (1|subject), data=datasubs)
Anova(SemConVsInc)
summary(glht(SemConVsInc, linfct=mcp(congruency="Tukey"))) # redo with multicomp

datasubs <- droplevels(data[data$condition %in% c(3),])#Semantic and list
datasubs$dataint1 <- interaction(datasubs$congruency, datasubs$list) # redo with subscripts and multicomp
SemConVsInc <- lmer(RT ~ dataint1 + (1|subject), data=datasubs)#Semantic and list
Anova(SemConVsInc)
summary(glht(SemConVsInc, linfct=mcp(dataint1="Tukey")))


###COLOR ALONE

datasubs <- droplevels(data[data$condition %in% c(4),])
ColConVsInc <- lmer(RT ~ congruency + (1|subject), data=datasubs)
Anova(ColConVsInc)
summary(glht(ColConVsInc, linfct=mcp(congruency="Tukey"))) # redo with multicomp

datasubs <- droplevels(data[data$condition %in% c(4),])#Color and list
datasubs$dataint1 <- interaction(datasubs$congruency, datasubs$list) # redo with subscripts and multicomp
ColConVsInc <- lmer(RT ~ dataint1 + (1|subject), data=datasubs)#Color and list
Anova(ColConVsInc)
summary(glht(ColConVsInc, linfct=mcp(dataint1="Tukey")))

###EMOTIONAL VS NEUTRAL STROOP


#Can we apply linear models?? Let's check the conditions (Linearity, absence of colinearity, 
# homoskedacity, Normality of residuals, Independence)

EmVsNor <- lmer(RT ~ condition + (1|subject), data=droplevels(data[data$condition %in% c(1, 2),]))

hist(residuals(EmVsNor)) # Residuals are normally distributed= "Tukey")))

qqnorm(residuals(EmVsNor)) # The quantiles progress sort of linearly too
plot(fitted(EmVsNor),residuals(EmVsNor)) ##There's some trouble with some of the answers... my money is on subject 14!

Anova(EmVsNor)
lsm.options(pbkrtest.limit = 9557)
lsmeans(EmVsNor, pairwise ~ condition, adjust="tukey")

summary(glht(EmVsNor, linfct=mcp(condition="Tukey"))) # redo with multicomp


###SEMANTIC VS NEUTRAL STROOP

datasubs <- droplevels(data[data$condition %in% c("Semantic", "Neutral"),])#subseting Normal trials as false congruent trials to enable comparison
SemVsNor <- lmer(RT ~ condition + (1|subject), data=droplevels(datasubs[datasubs$congruency %in% c(1, 2),])) #change 0 for 1 to test cong
Anova(SemVsNor)
summary(glht(SemVsNor, linfct=mcp(condition="Tukey"))) # redo with multicomp

###COLOR VS NEUTRAL STROOP

datasubs <- droplevels(data[data$condition %in% c("Neutral", "Color"),])#subseting Normal trials as false congruent trials to enable comparison,
ColVsNor <- lmer(RT ~ condition + (1|subject), droplevels(datasubs[datasubs$congruency %in% c(1, 2),])) #change 0 for 1 to test cong
Anova(ColVsNor)
summary(glht(ColVsNor, linfct=mcp(condition="Tukey"))) # redo with multicomp

###SEMANTIC VS EMOTIONAL STROOP

datasubs <- droplevels(data[data$condition %in% c("Semantic", "Emotional") & data$congruency %in% c(1,2),])#
datasubs <- droplevels(data[data$condition %in% c("Emotional", "Semantic"),])#subseting Neutral trials as false congruent trials to enable comparison

datasubs$condition <- factor(datasubs$condition, labels = c("Emotional","Semantic")) #Relabel for plot

data$condition <- factor(data$condition, labels = c("Neutral","Emotional","Semantic","Color")) #Relabel for plot
contrasts(datasubs$condition)=contr.sum(2)
contrasts(data$condition)=contr.sum(4)

SemVsEmr <- lmer(RT ~ condition + (1|subject), data=droplevels(datasubs[datasubs$congruency %in% c(0,2),])) #change 0 for 1 to test cong
SemVsEmr <- lmer(RT ~ condition + (1|subject), datasubs) 
Anova(SemVsEmr)
summary(SemVsEmr)
summary(glht(SemVsEmr, linfct=mcp(condition="Tukey"))) # redo with multicomp


###################### ###################### ######################
#NICER PLOTS#
###################### ###################### ######################
# data = subject session age handedness gender language condition color target  RT list accuracy congruency #
# For condition: 0-Training 1-Neutral 2-Emotionnal 3-Semantic 4-Color#
# For congruency: 1-Congruent 0-Incongruent 2-Nonsensical#

##Colors## (normal hues that you can get from the internet or any palette generator, like http://colorschemedesigner.com)

SemanticCongruent <- c("#DC0055")
SemanticInc <- c("#EE6B9E")
Emotional <- c("#9C02A7")
Neutral <- c("#00AF64")
ColorCongruent <- c("#006363")
ColorInc <- c("#33CCCC")

##Semantic vs. Color stroop (congruence as a factor)##

s <- evalq(aggregate(list(RT=RT), list(subject=subject), mean), 
           droplevels(data[data$condition %in% c(3, 4),]))
grandmean <- mean(s$RT)
e <- evalq(aggregate(list(RT=RT), list(condition=condition, congruency=congruency), mean), 
           droplevels(data[data$condition %in% c(3, 4),]))
f <- evalq(aggregate(list(RT=RT), list(condition=condition, congruency=congruency), 
                        function(x){sd(x)/sqrt(nrow(s))}),
           #function(x){(sd(x)/sqrt(nrow(s)) * sqrt(4/3))}), #with morey correction for a 2*2 design
           droplevels(data[data$condition %in% c(3, 4),])) 
e$se <- f$RT

e$condition <- factor(e$condition, labels = c("Semantic","Color")) #Relabel for plot
e$congruency <- factor(e$congruency, labels = c("Incongruent","Congruent") ) #Relabel for plot
ggplot(e, aes(x=condition, y=RT, color=congruency, group=congruency)) +  #Interaction condition:congruency SemVsCol
geom_line(size=2) +
geom_point(aes(shape=congruency), size=5) +
geom_errorbar(data=e, mapping=aes(x=condition, ymax=RT+se, ymin=RT-se), width=.05) +
scale_color_manual(values = c(SemanticCongruent,ColorCongruent)) +
geom_text(size=14, color="black", aes(x=1, y=720, label="n.s.")) +
geom_text(size=14, color="black", aes(x=2, y=760, label="***")) +
#geom_segment(size=2, color="black", aes(x=0.55,xend=1.05,y=700,yend=700)) +
labs(title = "Interaction between Congruency and Condition\n (Semantic Vs. Color stroop)", 
       x = "Stroop Type\n", 
       y = "RT (ms)") +  
theme_bw() +
  theme(legend.text = element_text(angle=0, face="bold", colour="black", size="14"),
        legend.title = element_text(angle=0, face="bold", colour="white"),
        panel.grid.major = element_blank(), # panel.grid.minor = element_blank(), 
        axis.line=element_line(size=1, color="black"), panel.border=element_blank(), 
        axis.text.x=element_text(size=14), axis.title.x=element_text(size=16),
        axis.text.y=element_text(size=14), axis.title.y=element_text(size=16),
        plot.title=element_text(size=20, face="bold", color="black"))

##Color stroop (congruence as a factor)##

s <- evalq(aggregate(list(RT=RT), list(subject=subject), mean), 
           droplevels(data[data$condition %in% c(4),]))
grandmean <- mean(s$RT)
e <- evalq(aggregate(list(RT=RT), list(congruency=congruency), mean), 
           droplevels(data[data$condition %in% c(4),]))
f <- evalq(aggregate(list(RT=RT), list(congruency=congruency), 
                     function(x){sd(x)/sqrt(nrow(s))}),
           #function(x){(sd(x)/sqrt(nrow(s)) * sqrt(4/3))}), #with morey correction for a 2*2 design
           droplevels(data[data$condition %in% c(4),])) 
e$se <- f$RT


e$congruency <- factor(e$congruency, labels = c("Incongruent","Congruent") ) #Relabel for plot
ggplot(e, aes(x=congruency, y=RT, color=congruency, group=congruency)) +  #Interaction condition:congruency SemVsCol
  geom_line(size=2) +
  geom_point(aes(shape=congruency), size=5) +
  geom_errorbar(data=e, mapping=aes(x=congruency, ymax=RT+se, ymin=RT-se), width=.05) +
  scale_color_manual(values = c(ColorInc, ColorCongruent)) +
  #geom_text(size=14, color="black", aes(x=1, y=720, label="n.s.")) +
  geom_text(size=14, color="black", aes(x=1.5, y=765, label="***")) +
  geom_segment(size=2, color="black", aes(x=1,xend=2,y=760,yend=760)) +
  labs(title = "Color stroop (Congruent vs. Incongruent)", 
       x = "Congruence\n", 
       y = "RT (ms)") +  
  theme_bw() +
  theme(legend.text = element_text(angle=0, face="bold", colour="black", size="14"),
        legend.title = element_text(angle=0, face="bold", colour="white"),
        panel.grid.major = element_blank(), # panel.grid.minor = element_blank(), 
        axis.line=element_line(size=1, color="black"), panel.border=element_blank(), 
        axis.text.x=element_text(size=14), axis.title.x=element_text(size=16),
        axis.text.y=element_text(size=14), axis.title.y=element_text(size=16),
        plot.title=element_text(size=20, face="bold", color="black"))


##Emotional Stroop vs. Neutral ##

s <- evalq(aggregate(list(RT=RT), list(subject=subject), mean), 
           droplevels(data[data$condition %in% c(1, 2),]))
grandmean <- mean(s$RT)
e <- evalq(aggregate(list(RT=RT), list(subject=subject, condition=condition), mean), 
           droplevels(data[data$condition %in% c(1, 2),]))
f <- evalq(aggregate(list(RT=RT), list(subject=subject, condition=condition), 
                     function(x){sd(x)/sqrt(nrow(s))}), droplevels(data[data$condition %in% c(1, 2),])) 
e$se <- f$RT

e$condition <- factor(e$condition, labels = c("Neutral (no stroop)","Emotional")) #Relabel for plot
ggplot(e, aes(x=condition, y=RT, color=condition, group=condition)) +  #Interaction condition:congruency SemVsCol
  #geom_line(size=2) +
  #geom_point(aes(shape=congruency), size=5) +
  geom_boxplot(size=1) +
  #geom_errorbar(data=e, mapping=aes(x=condition, ymax=RT+se, ymin=RT-se), width=.05) +
  scale_color_manual(values = c(Neutral,Emotional)) +
  geom_text(size=14, color="black", aes(x=1.5, y=900, label="***")) +
# geom_text(size=14, color="black", aes(x=2, y=760, label="**")) +
  geom_segment(size=2, color="black", aes(x=1,xend=2,y=880,yend=880)) +
  labs(title = "Emotional Stroop Vs. Neutral control\n", 
       x = "Stroop Type\n", 
       y = "RT (ms)") +  
  theme_bw() +
  theme(legend.text = element_text(angle=0, face="bold", colour="black", size="14"),
        legend.title = element_text(angle=0, face="bold", colour="white"),
        panel.grid.major = element_blank(), # panel.grid.minor = element_blank(), 
        axis.line=element_line(size=1, color="black"), panel.border=element_blank(), 
        axis.text.x=element_text(size=14), axis.title.x=element_text(size=16),
        axis.text.y=element_text(size=14), axis.title.y=element_text(size=16),
        plot.title=element_text(size=20, face="bold", color="black"))


##Semantic vs. Neutral##

s <- evalq(aggregate(list(RT=RT), list(subject=subject), mean), 
           droplevels(data[data$condition %in% c(1, 3),]))
grandmean <- mean(s$RT)
e <- evalq(aggregate(list(RT=RT), list(subject=subject, condition=condition), mean), 
           droplevels(data[data$condition %in% c(1, 3)
                           & data$congruency %in% c(1, 2),]))
f <- evalq(aggregate(list(RT=RT), list(subject=subject, condition=condition), 
                     function(x){sd(x)/sqrt(nrow(s))}), droplevels(data[data$condition %in% c(1, 3)
                           & data$congruency %in% c(1, 2),])) 
e$se <- f$RT

e$condition <- factor(e$condition, labels = c("Neutral (no stroop)","Semantic")) #Relabel for plot
ggplot(e, aes(x=condition, y=RT, color=condition, group=condition)) +  #Interaction condition:congruency SemVsCol
  #geom_line(size=2) +
  #geom_point(aes(shape=congruency), size=5) +
  geom_boxplot(size=1) +
  #geom_errorbar(data=e, mapping=aes(x=condition, ymax=RT+se, ymin=RT-se), width=.05) +
  scale_color_manual(values = c(Neutral, SemanticInc)) +
  geom_text(size=14, color="black", aes(x=1.5, y=900, label="***")) +
  # geom_text(size=14, color="black", aes(x=2, y=760, label="**")) +
  geom_segment(size=2, color="black", aes(x=1,xend=2,y=880,yend=880)) +
  labs(title = "Semantic Stroop (congruent)\nVs. Neutral control\n", 
       x = "Condition\n", 
       y = "RT (ms)") +  
  theme_bw() +
  theme(legend.text = element_text(angle=0, face="bold", colour="black", size="14"),
        legend.title = element_text(angle=0, face="bold", colour="white"),
        panel.grid.major = element_blank(), # panel.grid.minor = element_blank(), 
        axis.line=element_line(size=1, color="black"), panel.border=element_blank(), 
        axis.text.x=element_text(size=14), axis.title.x=element_text(size=16),
        axis.text.y=element_text(size=14), axis.title.y=element_text(size=16),
        plot.title=element_text(size=20, face="bold", color="black"))

##Emotional vs. Semantic Stroop##


s <- evalq(aggregate(list(RT=RT), list(subject=subject), mean), 
           droplevels(data[data$condition %in% c(2, 3),]))
grandmean <- mean(s$RT)
e <- evalq(aggregate(list(RT=RT), list(subject=subject, condition=condition), mean), 
           droplevels(data[data$condition %in% c(2, 3)
                           & data$congruency %in% c(0, 2),]))
f <- evalq(aggregate(list(RT=RT), list(subject=subject, condition=condition), 
                     function(x){sd(x)/sqrt(nrow(s))}), droplevels(data[data$condition %in% c(2, 3)
                                                                        & data$congruency %in% c(0, 2),])) 
e$se <- f$RT

e$condition <- factor(e$condition, labels = c("Emotional","Semantic")) #Relabel for plot
ggplot(e, aes(x=condition, y=RT, color=condition, group=condition)) +  #Interaction condition:congruency SemVsCol
  #geom_line(size=2) +
  #geom_point(aes(shape=congruency), size=5) +
  geom_boxplot(size=1) +
  #geom_errorbar(data=e, mapping=aes(x=condition, ymax=RT+se, ymin=RT-se), width=.05) +
  scale_color_manual(values = c(Emotional,SemanticCongruent)) +
  geom_text(size=14, color="black", aes(x=1.5, y=940, label="n.s.")) +
  # geom_text(size=14, color="black", aes(x=2, y=760, label="**")) +
  geom_segment(size=2, color="black", aes(x=1,xend=2,y=880,yend=880)) +
#  labs(title = "Semantic Stroop (incongruent)\nVs. Emotional Stroop\n", 
  labs(title = "Semantic Stroop (congruent)\nVs. Emotional Stroop\n", 
       x = "Condition\n", 
       y = "RT (ms)") +  
  theme_bw() +
  theme(legend.text = element_text(angle=0, face="bold", colour="black", size="14"),
        legend.title = element_text(angle=0, face="bold", colour="white"),
        panel.grid.major = element_blank(), # panel.grid.minor = element_blank(), 
        axis.line=element_line(size=1, color="black"), panel.border=element_blank(), 
        axis.text.x=element_text(size=14), axis.title.x=element_text(size=16),
        axis.text.y=element_text(size=14), axis.title.y=element_text(size=16),
        plot.title=element_text(size=20, face="bold", color="black"))


##Color Stroop vs. Neutral##

s <- evalq(aggregate(list(RT=RT), list(subject=subject), mean), 
           droplevels(data[data$condition %in% c(1, 4),]))
grandmean <- mean(s$RT)
e <- evalq(aggregate(list(RT=RT), list(subject=subject, condition=condition), mean), 
           droplevels(data[data$condition %in% c(1, 4)
                           & data$congruency %in% c(1, 2),]))
f <- evalq(aggregate(list(RT=RT), list(subject=subject, condition=condition), 
                     function(x){sd(x)/sqrt(nrow(s))}), droplevels(data[data$condition %in% c(1, 4)
                                                                        & data$congruency %in% c(1, 2),])) 
e$se <- f$RT

e$condition <- factor(e$condition, labels = c("Neutral","Color")) #Relabel for plot
ggplot(e, aes(x=condition, y=RT, color=condition, group=condition)) +  #Interaction condition:congruency SemVsCol
  #geom_line(size=2) +
  #geom_point(aes(shape=congruency), size=5) +
  geom_boxplot(size=1) +
  #geom_errorbar(data=e, mapping=aes(x=condition, ymax=RT+se, ymin=RT-se), width=.05) +
  scale_color_manual(values = c(Neutral,ColorInc)) +
  #geom_text(size=14, color="black", aes(x=1.5, y=955, label="***")) +
  geom_text(size=14, color="black", aes(x=1.5, y=890, label="n.s.")) +
  geom_segment(size=2, color="black", aes(x=1,xend=2,y=855,yend=855)) +
#  labs(title = "Color Stroop (incongruent)\nVs. Neutral\n",
  labs(title = "Color Stroop (congruent)\nVs. Neutral\n",
       x = "Condition\n", 
       y = "RT (ms)") +  
  theme_bw() +
  theme(legend.text = element_text(angle=0, face="bold", colour="black", size="14"),
        legend.title = element_text(angle=0, face="bold", colour="white"),
        panel.grid.major = element_blank(), # panel.grid.minor = element_blank(), 
        axis.line=element_line(size=1, color="black"), panel.border=element_blank(), 
        axis.text.x=element_text(size=14), axis.title.x=element_text(size=16),
        axis.text.y=element_text(size=14), axis.title.y=element_text(size=16),
        plot.title=element_text(size=20, face="bold", color="black"))

  




