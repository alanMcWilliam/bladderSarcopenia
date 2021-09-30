
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

bladderSarc$SMI <- bladderSarc$SM.Area/(bladderSarc$height^2)

######################################################################
### select out variables for analysis and clean

bladderClean <- bladderSarc %>%
  select(RIP., survival.months, SM.Density, SM.Area, SMI, height, weight, age, gender, PS, T.stage, ACE.27, cycles, pre.NAC.BMI, bladder.Ca.related., progression., PFS) %>%
  filter(T.stage != '1',
         T.stage != '2/3') %>%
  filter(ACE.27 != '3',
         ACE.27 != 'NK')


### dataframe with blood
bladderCleanBloods <- bladderSarc %>%
  select(SM.Density, SM.Area, SMI, albumin, gfr, Hb, neutrophil, lymphocyte, N.L.ratio)

stview(dfSummary(bladderClean))

summary(factor(bladderClean$T.stage))
summary(bladderClean$SMI)


######################################################################
### plot and visualise data

summary(bladderClean$SM.Density)
ggplot(data=bladderClean, aes(SM.Density)) + 
  geom_histogram(breaks=seq(-5,20, by = 1.5),
                 col = "skyblue", fill = "lightblue") +
  labs(title = "", x = "Muscle density" ) +
  theme(panel.background = element_blank())

summary(bladderClean$SMI)
ggplot(data=bladderClean, aes(SMI)) + 
  geom_histogram(breaks=seq(0,35, by = 1.5),
                 col = "skyblue", fill = "lightblue") +
  labs(title = "", x = "SMI" ) +
  theme(panel.background = element_blank())

ggplot(data=bladderClean, aes(x=SMI, fill=gender)) + 
  geom_histogram(alpha=0.6, position = 'identity',
                breaks=seq(0,35, by = 1)) +
  labs(title = "", x = "SMI" ) +
  theme(panel.background = element_blank())


summary(factor(bladderClean$PS))
tapply(bladderClean$SMI, factor(bladderClean$PS), summary)
PSaov <- aov(bladderClean$SMI~factor(bladderClean$PS))
summary(PSaov)
TukeyHSD(PSaov)

ggplot(data=bladderClean, aes(x=SMI, fill=factor(PS))) + 
  geom_histogram(alpha=0.4, position = 'identity',
                 breaks=seq(0,35, by = 1)) +
  labs(title = "", x = "SMI" ) +
  theme(panel.background = element_blank())

######################################################################
### analysis against overall survival

summary(bladderSarc$Muscle.area)
summary(bladderSarc$SM.Area)        
summary(bladderSarc$muscle.density)
summary(bladderSarc$SM.Density) 


uniCox <- coxph(Surv(time = survival.months, event = RIP.)~SMI, data = bladderClean)
summary(uniCox)
uniCox <- coxph(Surv(time = survival.months, event = bladder.Ca.related.)~SMI, data = bladderClean)
summary(uniCox)


uniCox <- coxph(Surv(time = survival.months, event = RIP.)~factor(PS), data = bladderClean)
summary(uniCox)

multiCox <- coxph(Surv(time = survival.months, event = RIP.)~SMI + age + factor(gender) + factor(T.stage) + factor(ACE.27) + pre.NAC.BMI + factor(PS), data = bladderClean)
summary(multiCox)
multiCox <- coxph(Surv(time = survival.months, event = RIP.)~SMI + age + factor(gender) + factor(T.stage) + pre.NAC.BMI, data = bladderClean)
summary(multiCox)

summary(factor(bladderClean$RIP.))
tapply(bladderClean$survival.months, bladderClean$RIP., summary)


### progression free survival
uniCox <- coxph(Surv(time = PFS, event = progression.)~SMI, data = bladderClean)
summary(uniCox)
multiCox <- coxph(Surv(time = PFS, event = progression.)~SMI + age + factor(gender) + factor(T.stage) + pre.NAC.BMI, data = bladderClean)
summary(multiCox)      

######################################################################
###





######################################################################
### correlations and plots against blood counts


plot(bladderCleanBloods$SM.Density, bladderCleanBloods$albumin)
plot(bladderCleanBloods$SM.Density, bladderCleanBloods$neutrophil)
plot(bladderCleanBloods$SM.Density, bladderCleanBloods$lymphocyte)
plot(bladderCleanBloods$SM.Density, bladderCleanBloods$N.L.ratio)

plot(bladderCleanBloods$SM.Area, bladderCleanBloods$albumin)
plot(bladderCleanBloods$SM.Area, bladderCleanBloods$neutrophil)
plot(bladderCleanBloods$SM.Area, bladderCleanBloods$lymphocyte)
plot(bladderCleanBloods$SM.Area, bladderCleanBloods$N.L.ratio)

plot(bladderCleanBloods$SMI, bladderCleanBloods$albumin)
plot(bladderCleanBloods$SMI, bladderCleanBloods$neutrophil)
plot(bladderCleanBloods$SMI, bladderCleanBloods$lymphocyte)
cor.test(bladderCleanBloods$SMI, bladderCleanBloods$lymphocyte)

plot(bladderCleanBloods$SMI, bladderCleanBloods$N.L.ratio)


#######################################################################
### go fishing

opt_cut_mean <- survminer::surv_cutpoint(bladderClean, time = "survival.months", event = "RIP.", "SMI",  minprop = 0.1, progressbar = TRUE) #smethod="logrank" set within)
summary(opt_cut_mean)
cat_mean <-survminer::surv_categorize(opt_cut_mean)
cat_split_mean <- survfit(Surv(time = survival.months, event = RIP.)~SMI, data = cat_mean)
ggsurvplot(cat_split_mean, risk.table = TRUE, conf.int = TRUE, surv.median.line = "hv", pval = TRUE, ncensor.plot = FALSE)

summary(bladderClean$SMI)
bladderClean$SMImedian <- bladderClean$SMI < median(bladderClean$SMI)
sarc <- survfit(Surv(time = survival.months, event = RIP.)~SMImedian, data = bladderClean)
ggsurvplot(sarc, risk.table = TRUE, conf.int = TRUE, surv.median.line = "hv", pval = TRUE, ncensor.plot = FALSE)

bladderClean <- within(bladderClean, quartileSMI <- as.integer(cut(SMI, quantile(SMI, probs=0:4/4, na.rm = TRUE)), include.lowest=TRUE))
tapply(bladderClean$SMI, bladderClean$quartileSMI, summary)
sarcQuart <- survfit(Surv(time = survival.months, event = RIP.)~quartileSMI, data = bladderClean)
ggsurvplot(sarcQuart, risk.table = TRUE, conf.int = TRUE, pval = FALSE, ncensor.plot = FALSE)

bladderClean$testQuart <- bladderClean$SMI < 15.8
sarcQuart2 <- survfit(Surv(time = survival.months, event = RIP.)~testQuart, data = bladderClean)
ggsurvplot(sarcQuart2, risk.table = TRUE, conf.int = TRUE, pval = TRUE, ncensor.plot = FALSE)


bladderCleanM <- bladderClean %>%
  filter(gender == 'M')
summary(bladderCleanM$SMI)

bladderCleanF <- bladderClean %>%
  filter(gender == 'F')
summary(bladderCleanF$SMI)

bladderCleanM$testQuart <- bladderCleanM$SMI < 17.5
sarcQuartM <- survfit(Surv(time = survival.months, event = RIP.)~testQuart, data = bladderCleanM)
ggsurvplot(sarcQuartM, risk.table = TRUE, conf.int = TRUE, pval = TRUE, ncensor.plot = FALSE)

bladderCleanF$testQuart <- bladderCleanF$SMI < 13.1
sarcQuartF <- survfit(Surv(time = survival.months, event = RIP.)~testQuart, data = bladderCleanF)
ggsurvplot(sarcQuartF, risk.table = TRUE, conf.int = TRUE, pval = TRUE, ncensor.plot = FALSE)
