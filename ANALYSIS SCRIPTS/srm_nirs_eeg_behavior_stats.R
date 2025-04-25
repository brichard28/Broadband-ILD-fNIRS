# Author: Benjamin Richardson
# Uses information from srm_nirs_eeg_analyze_behavior.m


library(tidyverse)
library(ggpubr)
library(ggplot2)
library(rstatix)
library(afex)
library(dplyr)

# SummarySE Function ####
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


####################################################
##    Hit Rates    ##
####################################################

hit_rates <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\RESULTS DATA\\SRM-NIRS-EEG-1_Hit_Rates.csv")

# Remove unneeded columns, put in long format
hit_rates$OriginalVariableNames <- array(0:29)
colnames(hit_rates) <- c("S","ITD50_Noise","ITD500_Noise","ILD70n_Noise","ILD10_Noise","ITD50_Speech","ITD500_Speech","ILD70n_Speech","ILD10_Speech")
hit_rates <- pivot_longer(hit_rates, cols=c("ITD50_Noise","ITD500_Noise","ILD70n_Noise","ILD10_Noise","ITD50_Speech","ITD500_Speech","ILD70n_Speech","ILD10_Speech"),
                          names_to = c("Spatialization","Masker"), names_sep = "_", values_to = "HitRate")

# Organize Factors
to.factor <- c('S','Masker','Spatialization')
hit_rates[, to.factor] <- lapply(hit_rates[, to.factor], as.factor)

# Summary Statistics
# hit_rates %>% group_by(Condition, Masker) %>% get_summary_stats(HitRate, type = "mean_sd")
# 
# # Boxplot
# bxp <- ggboxplot(hit_rates, x = "Condition", y = "HitRate", color = "Masker", palette = "jco")
# bxp

# Check for normality, remove outliers
# hit_rates %>% group_by(Condition, Masker) %>% shapiro_test(HitRate)

#hit_rates_no_outliers <- hit_rates %>%
#  group_by(Condition, Masker) %>%
#  identify_outliers("HitRate") %>%
#  filter(!is.outlier)



# Create a QQ plot
# ggqqplot(hit_rates, "HitRate", ggtheme = theme_bw()) + facet_grid(Condition ~ Masker, labeller = "label_both")


model_hitrate <- mixed(HitRate ~ Spatialization*Masker + (1|S),
                data= hit_rates, 
                control = lmerControl(optimizer = "bobyqa"), method = 'LRT')

model_hitrate

# Compare Speech vs. Noise

hit_rates$Masker <- relevel(hit_rates$Masker, "Speech")

posthoc_hitrate_speech_v_noise <- lmer(HitRate ~ Masker + (1|S),
                                data= hit_rates, 
                                control = lmerControl(optimizer = "bobyqa"))

summary(posthoc_hitrate_speech_v_noise)

# Compare All spatializations to each other within speech masker
# ITD50 as reference
hit_rates$Spatialization <- relevel(hit_rates$Spatialization, "ITD50")
posthoc_hitrate_itd50_speech <- lmer(HitRate ~ Spatialization + (1|S),
                              data= subset(hit_rates, Masker == "Speech"), 
                              control = lmerControl(optimizer = "bobyqa"))
summary(posthoc_hitrate_itd50_speech)

# ITD500 as reference
hit_rates$Spatialization <- relevel(hit_rates$Spatialization, "ITD500")
posthoc_hitrate_itd500_speech <- lmer(HitRate ~ Spatialization + (1|S),
                              data= subset(hit_rates, Masker == "Speech"), 
                              control = lmerControl(optimizer = "bobyqa"))
summary(posthoc_hitrate_itd500_speech)

# ILD70n as reference
hit_rates$Spatialization <- relevel(hit_rates$Spatialization, "ILD70n")
posthoc_hitrate_ild70n_speech <- lmer(HitRate ~ Spatialization + (1|S),
                              data= subset(hit_rates, Masker == "Speech"), 
                              control = lmerControl(optimizer = "bobyqa"))
summary(posthoc_hitrate_ild70n_speech)

# ILD10 as reference
hit_rates$Spatialization <- relevel(hit_rates$Spatialization, "ILD10")
posthoc_hitrate_ild10_speech <- lmer(HitRate ~ Spatialization + (1|S),
                                      data= subset(hit_rates, Masker == "Speech"), 
                                      control = lmerControl(optimizer = "bobyqa"))
summary(posthoc_hitrate_ild10_speech)

# Compare all spatializations to each other within noise masker
# ITD50 as reference
hit_rates$Spatialization <- relevel(hit_rates$Spatialization, "ITD50")
posthoc_hitrate_itd50_noise <- lmer(HitRate ~ Spatialization + (1|S),
                                     data= subset(hit_rates, Masker == "Noise"), 
                                     control = lmerControl(optimizer = "bobyqa"))
summary(posthoc_hitrate_itd50_noise)

# ITD500 as reference
hit_rates$Spatialization <- relevel(hit_rates$Spatialization, "ITD500")
posthoc_hitrate_itd500_noise <- lmer(HitRate ~ Spatialization + (1|S),
                                      data= subset(hit_rates, Masker == "Noise"), 
                                      control = lmerControl(optimizer = "bobyqa"))
summary(posthoc_hitrate_itd500_noise)

# ILD70n as reference
hit_rates$Spatialization <- relevel(hit_rates$Spatialization, "ILD70n")
posthoc_hitrate_ild70n_noise <- lmer(HitRate ~ Spatialization + (1|S),
                                      data= subset(hit_rates, Masker == "Noise"), 
                                      control = lmerControl(optimizer = "bobyqa"))
summary(posthoc_hitrate_ild70n_noise)

# ILD10 as reference
hit_rates$Spatialization <- relevel(hit_rates$Spatialization, "ILD10")
posthoc_hitrate_ild10_noise <- lmer(HitRate ~ Spatialization + (1|S),
                                     data= subset(hit_rates, Masker == "Noise"), 
                                     control = lmerControl(optimizer = "bobyqa"))
summary(posthoc_hitrate_ild10_noise)





####################################################
##    FA Rates    ##
####################################################

FA_rates <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\RESULTS DATA\\SRM-NIRS-EEG-1_FA_Rates.csv")

# Remove unneeded columns, put in long format
FA_rates$OriginalVariableNames <- array(1:30)
colnames(FA_rates) <- c("S","ITD50","ITD500","ILD70n","ILD10")
FA_rates <- pivot_longer(FA_rates, cols=c("ITD50","ITD500","ILD70n","ILD10"),
                          names_to = c("Spatialization"), values_to = "FARate")

# Organize Factors
to.factor <- c('S','Spatialization')
FA_rates[, to.factor] <- lapply(FA_rates[, to.factor], as.factor)

# LMEM
model_farate <- mixed(FARate ~ Spatialization + (1|S),data= FA_rates,control = lmerControl(optimizer = "bobyqa"), method = 'LRT')

model_farate

# Post hocs
# ITD50 as reference
FA_rates$Spatialization <- relevel(FA_rates$Spatialization, "ITD50")
posthoc_farate_itd50 <- lmer(FARate ~ Spatialization + (1|S),
                        data= FA_rates, 
                        control = lmerControl(optimizer = "bobyqa"))

summary(posthoc_farate_itd50)

# ITD500 as reference
FA_rates$Spatialization <- relevel(FA_rates$Spatialization, "ITD500")
posthoc_farate_itd500 <- lmer(FARate ~ Spatialization + (1|S),
                             data= FA_rates, 
                             control = lmerControl(optimizer = "bobyqa"))

summary(posthoc_farate_itd500)

# ILD70n as reference
FA_rates$Spatialization <- relevel(FA_rates$Spatialization, "ILD70n")
posthoc_farate_ild70n <- lmer(FARate ~ Spatialization + (1|S),
                              data= FA_rates, 
                              control = lmerControl(optimizer = "bobyqa"))

summary(posthoc_farate_ild70n)

# ILD10 as reference
FA_rates$Spatialization <- relevel(FA_rates$Spatialization, "ILD10")
posthoc_farate_ild10 <- lmer(FARate ~ Spatialization + (1|S),
                              data= FA_rates, 
                              control = lmerControl(optimizer = "bobyqa"))

summary(posthoc_farate_ild10)


####################################################
##    D primes    ##
####################################################

d_primes <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\RESULTS DATA\\SRM-NIRS-EEG-1_d_primes.csv")

# Remove unneeded columns, put in long format
d_primes$OriginalVariableNames <- array(1:30)
colnames(d_primes) <- c("S","ITD50","ITD500","ILD70n","ILD10")
d_primes <- pivot_longer(d_primes, cols=c("ITD50","ITD500","ILD70n","ILD10"),
                         names_to = c("Spatialization"), values_to = "DPrime")

# Organize Factors
to.factor <- c('S','Spatialization')
d_primes[, to.factor] <- lapply(d_primes[, to.factor], as.factor)

# LMEM
model_dprime <- mixed(DPrime ~ Spatialization + (1|S),data= d_primes,control = lmerControl(optimizer = "bobyqa"))

model_dprime


# Post hocs
# ITD50 as reference
d_primes$Spatialization <- relevel(d_primes$Spatialization, "ITD50")
posthoc_dprime_itd50 <- lmer(DPrime ~ Spatialization + (1|S),
                             data= d_primes, 
                             control = lmerControl(optimizer = "bobyqa"))

summary(posthoc_dprime_itd50)

# ITD500 as reference
d_primes$Spatialization <- relevel(d_primes$Spatialization, "ITD500")
posthoc_dprime_itd500 <- lmer(DPrime ~ Spatialization + (1|S),
                             data= d_primes, 
                             control = lmerControl(optimizer = "bobyqa"))

summary(posthoc_dprime_itd500)

# ILD70n as reference
d_primes$Spatialization <- relevel(d_primes$Spatialization, "ILD70n")
posthoc_dprime_ild70n <- lmer(DPrime ~ Spatialization + (1|S),
                             data= d_primes, 
                             control = lmerControl(optimizer = "bobyqa"))

summary(posthoc_dprime_ild70n)

# ILD10 as reference
d_primes$Spatialization <- relevel(d_primes$Spatialization, "ILD10")
posthoc_dprime_ild10 <- lmer(DPrime ~ Spatialization + (1|S),
                              data= d_primes, 
                              control = lmerControl(optimizer = "bobyqa"))

summary(posthoc_dprime_ild10)












####################################################
##    Object Rates    ##
####################################################
# 
object_rates <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\RESULTS DATA\\SRM-NIRS-EEG-1_object_Rates.csv")

# Remove unneeded columns, put in long format
object_rates$OriginalVariableNames <- array(1:30)
colnames(object_rates) <- c("S","ITD50_Noise","ITD500_Noise","ILD70n_Noise","ILD10_Noise","ITD50_Speech","ITD500_Speech","ILD70n_Speech","ILD10_Speech")
object_rates <- pivot_longer(object_rates, cols=c("ITD50_Noise","ITD500_Noise","ILD70n_Noise","ILD10_Noise","ITD50_Speech","ITD500_Speech","ILD70n_Speech","ILD10_Speech"),
                          names_to = c("Spatialization","Masker"), names_sep = "_", values_to = "ObjectRate")

# Organize Factors
to.factor <- c('S','Masker','Spatialization')
object_rates[, to.factor] <- lapply(object_rates[, to.factor], as.factor)

# Summary Statistics
object_rates %>% group_by(Spatialization, Masker) %>% get_summary_stats(ObjectRate, type = "mean_sd")


# LMEM
model_objectrate <- mixed(ObjectRate ~ Spatialization*Masker + (1|S),
                      data= object_rates, 
                      control = lmerControl(optimizer = "bobyqa"), method = 'LRT')

model_objectrate

# Post hocs

# ITD50 as reference
object_rates$Spatialization <- relevel(object_rates$Spatialization, "ITD50")
posthoc_objectrate_itd50_speech <- lmer(ObjectRate ~ Spatialization + (1|S),
                             data= subset(object_rates, Masker == "Speech"), 
                             control = lmerControl(optimizer = "bobyqa"))

summary(posthoc_objectrate_itd50_speech)

# ITD500 as reference
object_rates$Spatialization <- relevel(object_rates$Spatialization, "ITD500")
posthoc_objectrate_itd500_speech <- lmer(ObjectRate ~ Spatialization + (1|S),
                              data= subset(object_rates, Masker == "Speech"), 
                              control = lmerControl(optimizer = "bobyqa"))

summary(posthoc_objectrate_itd500_speech)

# ILD70n as reference
object_rates$Spatialization <- relevel(object_rates$Spatialization, "ILD70n")
posthoc_objectrate_ild70n_speech <- lmer(ObjectRate ~ Spatialization + (1|S),
                              data= subset(object_rates, Masker == "Speech"), 
                              control = lmerControl(optimizer = "bobyqa"))

summary(posthoc_objectrate_ild70n_speech)

# ILD10 as reference
object_rates$Spatialization <- relevel(object_rates$Spatialization, "ILD10")
posthoc_objectrate_ild10_speech <- lmer(ObjectRate ~ Spatialization + (1|S),
                             data= subset(object_rates, Masker == "Speech"), 
                             control = lmerControl(optimizer = "bobyqa"))

summary(posthoc_objectrate_ild10_speech)


# ITD50 as reference
object_rates$Spatialization <- relevel(object_rates$Spatialization, "ITD50")
posthoc_objectrate_itd50_noise <- lmer(ObjectRate ~ Spatialization + (1|S),
                                        data= subset(object_rates, Masker == "Noise"), 
                                        control = lmerControl(optimizer = "bobyqa"))

summary(posthoc_objectrate_itd50_noise)

# ITD500 as reference
object_rates$Spatialization <- relevel(object_rates$Spatialization, "ITD500")
posthoc_objectrate_itd500_noise <- lmer(ObjectRate ~ Spatialization + (1|S),
                                         data= subset(object_rates, Masker == "Noise"), 
                                         control = lmerControl(optimizer = "bobyqa"))

summary(posthoc_objectrate_itd500_noise)

# ILD70n as reference
object_rates$Spatialization <- relevel(object_rates$Spatialization, "ILD70n")
posthoc_objectrate_ild70n_noise <- lmer(ObjectRate ~ Spatialization + (1|S),
                                         data= subset(object_rates, Masker == "Noise"), 
                                         control = lmerControl(optimizer = "bobyqa"))

summary(posthoc_objectrate_ild70n_noise)

# ILD10 as reference
object_rates$Spatialization <- relevel(object_rates$Spatialization, "ILD10")
posthoc_objectrate_ild10_noise <- lmer(ObjectRate ~ Spatialization + (1|S),
                                        data= subset(object_rates, Masker == "Noise"), 
                                        control = lmerControl(optimizer = "bobyqa"))

summary(posthoc_objectrate_ild10_noise)