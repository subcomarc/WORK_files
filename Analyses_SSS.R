##Analysis Script for the SSS French validation
##By subcomarc, with the resources of Jérôme Sackur
## Needs file "Resultats_SSS.xlsx"

############### ################################
#Preloading settings, directories and libraries#
############### ################################
source("R_Lib_Sackur.r")

library("nlme")
library('xtable')
library('knitr')
library('lme4')
library('car')
library('multcomp')
library('lsmeans')
library('emmeans')
library("Hmisc")
library("openxlsx")
library("xlsx")
library("ggplot2")


options(xtable.floating = FALSE)
options(xtable.timestamp = "")

##################### ###################################################
#Load preprocessed data (if this is the first time, DO PREPPING INSTEAD)#
#################### ####################################################

load("/home/subcomarc/WORK/R/SSS_Data.R")

####################### ########################################
#Prepping the data for stats and adding some convenient columns#
###################### ######################################### ##################

# Q1 : Je suis facilement influencé-e par l'opinion des autres.	
# Q2 : Je trouve beaucoup de conseils pratiques à la télévision ou dans des magazines.	
# Q3: Quand quelqu'un tousse ou éternue, j'ai généralement envie de faire de même.	
# Q4: Imaginer une boisson rafraîchissante peut me donner soif.	
# Q5 : Un bon vendeur peut vraiment me donner envie de son produit.	
# Q6 : J'ai repris de nombreuses habitudes de mes ami-es	
# Q7 : Pour moi, il est important de bien s'intégrer dans le groupe dont on fait partie.	
# Q8 : Quand je vois quelqu'un frissonner, souvent cela me donne froid.	
# Q9 : Mon style s'inspire de celui de certaines célébrités.	
# Q10 : Quand les gens me racontent ce qu'ils ressentent, je remarque souvent que je ressens la même chose.	
# Q11 : Lorsque je prends une décision, je suis souvent les conseils des autres.	
# Q12 : Lire des descriptions de plats savoureux me donne l'eau à la bouche.	
# Q13 : Je puise beaucoup de bonnes idées chez les autres.	
# Q14 : Il m'arrive de me laisser influencer par une bonne publicité.	
# Q15 : Après avoir vu une publicité pour une crème hydratante, parfois ma peau me semble sèche.	
# Q16 : J'ai découvert beaucoup des choses que je préfère grâces à mes ami-es.	
# Q17 : Quand un produit est joliment présenté, je veux généralement l'acheter.	
# Q18 : Quand je pense à quelque chose d'effrayant, je sens mon cœur se mettre à battre plus fort.	
# Q19 : Je change souvent d'opinion après avoir discuté avec d'autres personnes.	
# Q20 : Si on me dit que je n'ai pas l'air bien, je commence à me sentir malade.	
# Q21 : Je suis les tendances de la mode.	
# 	
# 	
# Sessions : 	1 = Mardi 10h
# 	2 = Jeudi 18h30
# 	3 = Vendredi 10h


#Load from Excel

dataRAW <- read.xlsx("/home/subcomarc/WORK/R/Resultats_SSS.xlsx", 1) #This alone should be enough to import everything as numeric,
                                                                    #from Windows when in linux, BUT it is known to fail 
                                                                    #because of java, so we have to manually
                                                                    #make sure that NaNs are NaNs and numbers are numbers. 
                                                                    #Hence:
data <- dataRAW
for (i in 1:ncol(data)) {
  data[[i]] <- gsub ("NA", NA, data[[i]])
}
data[, c(2:23)] <- matrix(as.numeric(unlist(data[, c(2:23)])))


# Add a column with the actual score for each participant - The score is just the sum of all items (21-105) #
# For now I'm leaving everyone in, including those who missed a question or two. I'd remove them in the future tho...

data$Score <- rowSums(data[, c(3:23)], na.rm = 1)
#data$Score <- rowSums(dataRAW[, c(3:23)], na.rm = 0)


#Save as R file to play later
save(data, file='SSS_Data.R')

###################### ##########################
#Quick stats - descriptive and some interactions#
##################### ###########################

hist(data$Score) # Nice and rather normal dist, consider I did not do any cleaning of the data at this point, 
                  #we may need the outliers as Highs and Lows

plot(density(data$Score))
#lines(density(data$Score2), col = 2) #in red

# Differences between sessions?
data$Session <- as.factor(data$Session)
Anova(lm(Score ~ Session, data)) #Nope, but surprisingly not that far. Did something happen in one of the sessions?
                                  #By the way, why lm and not lmer? Because I only have one observation per participant
