# Author: Benjamin Richardson
# Uses information from srm_nirs_eeg_analyze_behavior.m


library(tidyverse)
library(ggpubr)
library(ggplot2)
library(rstatix)
library(afex)
library(dplyr)

####################################################
##    Hit Rates    ##
####################################################

hit_rates <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\RESULTS DATA\\SRM-NIRS-EEG-1_Hit_Rates.csv")

# Remove unneeded columns, put in long format
hit_rates$OriginalVariableNames <- array(1:30)
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
model_farate <- mixed(FARate ~ Spatialization + (1|S),
                       data= FA_rates, 
                       control = lmerControl(optimizer = "bobyqa"), method = 'LRT')

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

# d_primes <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\RESULTS DATA\\SRM-NIRS-EEG-1_d_primes.csv")
# 
# # Remove unneeded columns, put in long format
# d_primes$OriginalVariableNames <- array(1:30)
# colnames(d_primes) <- c("S","ITD50","ITD500","ILD70n","ILD10")
# d_primes <- pivot_longer(d_primes, cols=c("ITD50","ITD500","ILD70n","ILD10"),
#                          names_to = c("Condition"), values_to = "DPrime")
# 
# # Organize Factors
# to.factor <- c('S','Condition')
# d_primes[, to.factor] <- lapply(d_primes[, to.factor], as.factor)
# 
# # Summary Statistics
# d_primes %>% group_by(Condition) %>% get_summary_stats(DPrime, type = "mean_sd")
# 
# # Boxplot
# bxp <- ggboxplot(d_primes, x = "Condition", y = "DPrime", palette = "jco")
# bxp
# 
# # Check for normality, remove outliers
# d_primes %>% group_by(Condition) %>% shapiro_test(DPrime)
# 
# # Create a QQ plot
# ggqqplot(d_primes, "DPrime", ggtheme = theme_bw()) # + facet_grid(Condition) #, labeller = "label_both")
# 
# # Run ANOVA
# res.aov <- anova_test(data = d_primes, dv = DPrime, wid = S, within = c(Condition))
# get_anova_table(res.aov)
# 
# # Pairwise comparisons between conditions
# pwc_masker <- d_primes  %>% pairwise_t_test(DPrime ~ Condition, paired = TRUE, p.adjust.method = "bonferroni")
# print(pwc_masker)
# 
# 
# 
# 
















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