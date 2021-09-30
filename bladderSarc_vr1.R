
######################################################################
###  Analysis of impact of bladder sarcopenia
###  Alan McWilliam
###  29th Sept 2021
###
###
######################################################################

library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(summarytools)
library(survival)
library(KMsurv)
library(survminer)

######################################################################
### read data and join
### filter out incorrect segmentations

bladderClin <- read.csv("C:/Users/alan_/Desktop/bladderSarcopenia/bladderClinData.csv")
bladderStats <- read.csv("C:/Users/alan_/Desktop/bladderSarcopenia/bladderSarc_stats_20210929.csv")
bladderJoin <- read.csv("C:/Users/alan_/Desktop/bladderSarcopenia/NewBigList.csv")
bladderRemove <- read.csv("C:/Users/alan_/Desktop/bladderSarcopenia/remove.csv")

bladderClin <- bladderClin %>%
  rename(ID = "Christie.no.")
bladderJoin <- bladderJoin %>%
  distinct(Number, .keep_all = TRUE)
bladderSarc <- merge(bladderJoin, bladderClin, by = "ID")

bladderStats$ID <- str_remove(bladderStats$ID, '.xdr')
bladderStats <- bladderStats %>%
  rename(joinID = "ID")
bladderSarc$joinID <- paste(bladderSarc$anonymisation, bladderSarc$Number, sep = "")
bladderSarc <- merge(bladderSarc, bladderStats, by = 'joinID')

bladderSarc <- bladderSarc %>%
  filter(!Number %in% bladderRemove$ptNumber)
View(bladderSarc)

######################################################################
### select out variables for analysis and clean

bladderClean <- bladderSarc %>%
  select(RIP., survival.months, SM.Density, SM.Area, height, weight, age, gender, PS, T.stage, ACE.27) %>%
  filter(T.stage != '1',
         T.stage != '2/3')


stview(dfSummary(bladderClean))

summary(factor(bladderClean$T.stage))

######################################################################
### 

uniCox <- coxph(Surv(time = survival.months, event = RIP.)~SM.Area, data = bladderClean)
summary(uniCox)

multiCox <- coxph(Surv(time = survival.months, event = RIP.)~SM.Area + weight + age + factor(gender) + factor(T.stage), data = bladderClean)
summary(multiCox)

summary(factor(bladderClean$RIP.))
tapply(bladderClean$survival.months, bladderClean$RIP., summary)

summary(bladderSarc$Muscle.area)
summary(bladderSarc$SM.Area)        
summary(bladderSarc$muscle.density)
summary(bladderSarc$SM.Density)        
