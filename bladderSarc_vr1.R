
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
### Bootstrapping
### LASSO elastic net 
#https://www.r-bloggers.com/variable-selection-with-elastic-net/
#https://cran.r-project.org/web/packages/glmnet/glmnet.pdf

### set-up datafrme for bootstrapping 
## RIP., survival.months, bladder.Ca.related  --- for survival
## PFS, progression.  --- for disese progression
bladderCleanBootstrap <- bladderSarc %>%
  select(SM.Density, SM.Area, SMI, weight, age, gender, PS, T.stage, ACE.27, pre.NAC.BMI, progression.) %>%
  filter(T.stage != '1',
         T.stage != '2/3') %>%
  filter(ACE.27 != '3',
         ACE.27 != 'NK')
bladderCleanBootstrap$PS <- as.factor(bladderCleanBootstrap$PS)
stview(dfSummary(bladderCleanBootstrap))

pkgs <- list("glmnet", "doParallel", "foreach", "pROC")
lapply(pkgs, require, character.only = T)
registerDoParallel(cores = 4)

performElasticNet <- function(dataAll){
  #set coefficents for glm 
  mdlY <- as.factor(as.matrix(dataAll["progression."]))  ## need to update for event column
  mdlX <- model.matrix(~.-1, dataAll[-ncol(dataAll)])
  
  coefOut <- matrix(NA, nrow = 3, ncol = 22)  ## this needs updated for size of dataframe, match below
  
  #full LASSO
  cv1 <- cv.glmnet(mdlX, mdlY, family = "binomial", nfold = 10, type.measure = "deviance", parallel = TRUE, alpha = 1)
  md1 <- glmnet(mdlX, mdlY, family = "binomial", lambda = cv1$lambda.1se, alpha = 1)
  tmp_LASSO <- data.frame(coef.name = dimnames(coef(md1))[[1]], coef.value = matrix(coef(md1)))
  ##print(tmp_LASSO)
  coefOut[1,] <- tmp_LASSO[,2]
  
  #rigid
  md2 <- glmnet(mdlX, mdlY, family = "binomial", lambda = cv1$lambda.1se, alpha = 0) 
  tmp_rigid <- data.frame(coef.name = dimnames(coef(md2))[[1]], coef.value = matrix(coef(md2)))
  coefOut[2,] <- tmp_rigid[,2]
  
  #elastic net LASSO
  a <- seq(0.1, 0.9, 0.05) #change final number for fine tuning to be faster was 0.01
  search <- foreach(i = a, .combine = rbind) %dopar% {
    cv <- cv.glmnet(mdlX, mdlY, family = "binomial", nfold = 10, type.measure = "deviance", paralle = TRUE, alpha = i)
    data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
  }
  cv3 <- search[search$cvm == min(search$cvm), ]
  md3 <- glmnet(mdlX, mdlY, family = "binomial", lambda = cv3$lambda.1se, alpha = cv3$alpha)
  tmp_elasticNet <- data.frame(coef.name = dimnames(coef(md3))[[1]], coef.value = matrix(coef(md3)))
  coefOut[3,] <- tmp_elasticNet[,2]
  
  return (coefOut)
}

# B number of bootstraps
bootstrap_r <- function(ds, B) {
  #ds <- all2_EN
  # Preallocate storage for statistics
  boot_stat_LASSO <- matrix(NA, nrow = B, ncol = 22) #number for ncol - number coefficients plus one for intercept
  boot_stat_rigid <- matrix(NA, nrow = B, ncol = 22) #number for ncol 
  boot_stat_elasticNet <- matrix(NA, nrow = B, ncol = 22) #number for ncol 
  
  # Number of observations
  n <- nrow(ds)
  
  # Perform bootstrap
  for(i in seq_len(B)) {
    print(i)
    # Sample initial data
    gen_data <- ds[ sample(nrow(ds), replace=TRUE), ]
    # Calculate sample data mean and SD
    coefOut2 <- performElasticNet(gen_data)
    #print(coefOut2)
    
    
    boot_stat_LASSO[i,] <- coefOut2[1,]
    boot_stat_rigid[i,] <- coefOut2[2,]
    boot_stat_elasticNet[i,] <- coefOut2[3,]
  }
  
  boot_stat <- rbind(boot_stat_LASSO, boot_stat_rigid, boot_stat_elasticNet)
  
  # Return bootstrap result
  return(boot_stat)
}


#set a seed for bootstrapping
set.seed(883)
b = 10
resultsAll <- bootstrap_r(bladderCleanBootstrap, b)
View(resultsAll)

#split into results for each model, add colum names and convert to data frames
boot_LASSO <- resultsAll[1:b,]
boot_rigid <- resultsAll[(b+1):(2*b),]
boot_elastic <- resultsAll[(2*b+1):nrow(resultsAll),]

colnames(boot_rigid) <- c("Intercept", "Intercept", "SM.Density", "SM.Area", "SMI", "weight", "age", "genderF", "GenderM", "PS1", "PS2", "T.stage1", "T.stage2", "T.stage2/3", "T.stage3", "T.stage4", "ACE.270", "ACE.271", "ACE.272", "ACE.273", "ACE.27NK", "pre.NAC.BMI")
colnames(boot_LASSO) <- c("Intercept", "Intercept", "SM.Density", "SM.Area", "SMI", "weight", "age", "genderF", "GenderM", "PS1", "PS2", "T.stage1", "T.stage2", "T.stage2/3", "T.stage3", "T.stage4", "ACE.270", "ACE.271", "ACE.272", "ACE.273", "ACE.27NK", "pre.NAC.BMI")
colnames(boot_elastic) <- c("Intercept", "Intercept", "SM.Density", "SM.Area", "SMI", "weight", "age", "genderF", "GenderM", "PS1", "PS2", "T.stage1", "T.stage2", "T.stage2/3", "T.stage3", "T.stage4", "ACE.270", "ACE.271", "ACE.272", "ACE.273", "ACE.27NK", "pre.NAC.BMI")

boot_rigid <- as.data.frame(boot_rigid)
boot_LASSO <- as.data.frame(boot_LASSO)
boot_elastic <- as.data.frame(boot_elastic)

View(boot_rigid)
View(boot_LASSO)
View(boot_elastic)


generateStats <- function(df){
  for(i in colnames(df)){
    print(i)
    print(mean(df[[i]]))
    print(summary(df[[i]]))
    print(quantile(df[[i]], probs = c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)))
    print(sum(df[[i]] > 0))
    
    #plot histogram
    ##LASSO y 0, 150, seq -0.05, 0.05
    #ggplot(data=df, aes(df[[i]])) + 
    #geom_histogram(breaks=seq(-0.05,0.05, by = 0.001),
    #               col = "skyblue", fill = "lightblue") +
    #labs(title = i, x = "coefficent" ) +
    #ylim(0,500) +
    #theme(panel.background = element_blank())
    
    # ggsave(paste("C:\\Users\\alan_\\Desktop\\templateHeart\\results3\\", i, ".jpg", sep=""))
  }
}


#read data in from bootstarpping
boot_rigid <- read.csv("C:\\Users\\alan_\\Desktop\\templateHeart\\bootstrapped\\bootRigid500.csv")
boot_LASSO <- read.csv("C:\\Users\\alan_\\Desktop\\templateHeart\\bootstrapped\\bootLASSO500.csv")
boot_elastic <- read.csv("C:\\Users\\alan_\\Desktop\\templateHeart\\bootstrapped\\bootElastic500.csv")


#how to select variables from bootstrapping
generateStats(boot_elastic)




summary(bladderSarc$Muscle.area)
summary(bladderSarc$SM.Area)        
summary(bladderSarc$muscle.density)
summary(bladderSarc$SM.Density)       

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
