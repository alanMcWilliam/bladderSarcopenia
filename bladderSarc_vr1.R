
######################################################################
###  Analysis of impact of bladder sarcopenia
###  Alan McWilliam
###  29th Sept 2021
###
######################################################################

library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(summarytools)


######################################################################
### read data and join
### filter out incorrect segmentations

bladderClin <- read.csv("C:/Users/alan_/Desktop/bladderSarcopenia/bladderClinData.csv")
bladderStats <- read.csv("C:/Users/alan_/Desktop/bladderSarcopenia/bladderSarc_stats_20210929.csv")
bladderJoin <- read.csv("C:/Users/alan_/Desktop/bladderSarcopenia/NewBigList.csv")
bladderRemove <- read.csv("C:/Users/alan_/Desktop/bladderSarcopenia/remove.csv")

bladderClin <- bladderClin %>%
  rename(ID = "Christie.no.")
bladderSarc <- merge(bladderClin, bladderJoin, by = "ID")

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
  select(RIP., survival.months, muscle.density, Muscle.area, height, weight, age, gender, PS, T.stage, ACE.27)
