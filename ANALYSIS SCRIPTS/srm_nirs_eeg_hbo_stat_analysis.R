
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(rstatix)
library(afex)
library(dplyr)
library(car)

# Load in Data
speech_masker_data <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\all_subjects_mean_during_stim_speech_masker.csv")
noise_masker_data <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\all_subjects_mean_during_stim_noise_masker.csv")

colnames(speech_masker_data) <- c("S","Channel","ITD50","ITD500","ILD70n","ILD10")
colnames(noise_masker_data) <- c("S","Channel","ITD50","ITD500","ILD70n","ILD10")
speech_masker_data <- pivot_longer(speech_masker_data, cols=c("ITD50","ITD500","ILD70n","ILD10"), names_to = "Spatialization", values_to = "MeanHbO")
noise_masker_data <- pivot_longer(noise_masker_data, cols=c("ITD50","ITD500","ILD70n","ILD10"), names_to = "Spatialization", values_to = "MeanHbO")


# Combine speech and noise
speech_masker_data$Masker<-"speech"
noise_masker_data$Masker<-"noise"
all_data<-rbind(speech_masker_data,noise_masker_data)

# Add ROI information
all_data$Roi<-NA
pfc_channels <- c(0,1,2,3,4,5)
stg_channels <- c(6,7,8,9,10,11,12,13)
all_data$Roi[which(all_data$Channel %in% pfc_channels)] <- "pfc"
all_data$Roi[which(all_data$Channel %in% stg_channels)] <- "stg"

# Change ROI information

#all_data$Roi[all_data$Roi == "0"] <- "pfc"
#all_data$Roi[all_data$Roi == "1"] <- "stg"

# Organize Factors
to.factor <- c('S', 'Roi','Spatialization')
all_data[, to.factor] <- lapply(all_data[, to.factor], as.factor)
#noise_masker_data[, to.factor] <- lapply(noise_masker_data[, to.factor], as.factor)




all_data$Masker <- as.factor(all_data$Masker)
all_data_cleaned <- na.omit(all_data)

# Take mean across ROI
r1<-with(all_data_cleaned, tapply(MeanHbO, Roi, mean))


# Boxplot PFC
#bxp_pfc <- ggboxplot(subset(all_data_cleaned,Roi=="pfc"), x = "Masker", y = "MeanHbO", color = "Spatialization", palette = "jco")
#bxp_pfc

# Boxplot STG
#bxp_stg <- ggboxplot(subset(all_data_cleaned,Roi=="stg"), x = "Masker", y = "MeanHbO", color = "Spatialization", palette = "jco")
#bxp_stg

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
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



# Scatterplot Speech vs. Noise
#ggplot(all_data_cleaned) +
#  geom_point(aes(x = Masker, y= MeanHbO)) + ggtitle("Speech vs. Noise Masker")

# Scatterplot PFC vs. STG
#ggplot(all_data_cleaned) +
#  geom_point(aes(x = Roi, y= MeanHbO)) + ggtitle("PFC vs. STG")

# Scatterplot across Spatialization in speech masker, PFC
#ggplot(subset(all_data_cleaned, Masker == "speech" & Roi == "pfc")) +
#  geom_point(aes(x = Spatialization, y= MeanHbO)) + ggtitle("Speech Masker, PFC")

# Scatterplot across Spatialization in speech masker, STG
#ggplot(subset(all_data_cleaned, Masker == "speech" & Roi == "stg")) +
#  geom_point(aes(x = Spatialization, y= MeanHbO)) + ggtitle("Speech Masker, STG")

# Scatterplot across Spatialization in noise masker, PFC
#ggplot(subset(all_data_cleaned, Masker == "noise" & Roi == "pfc")) +
#  geom_point(aes(x = Spatialization, y= MeanHbO)) + ggtitle("Noise Masker, PFC")

# Scatterplot across Spatialization in noise masker, STG
#ggplot(subset(all_data_cleaned, Masker == "noise" & Roi == "stg")) +
#  geom_point(aes(x = Spatialization, y= MeanHbO)) + ggtitle("Noise Masker, STG")

# Check for normality, remove outliers
all_data_cleaned %>% group_by(Spatialization, Masker,Roi) %>% shapiro_test(MeanHbO)


## ANOVA with factors of Masker Type, ROI, and spatial Spatialization

# Run ANOVA
#res.aov <- anova_test(data = all_data_cleaned, dv = MeanHbO, wid = S, within = c(Spatialization, Masker, Roi))
#get_anova_table(res.aov)

#perform three-way ANOVA
#model <- aov(MeanHbO ~ Spatialization * Masker * Roi, data=all_data_cleaned)

#view summary of three-way ANOVA
#summary(model)


# Pairwise comparisons between spatializations, gpiven maskers
#pwc_speech_masker <- subset(all_data_cleaned, Masker == "speech")  %>% pairwise_t_test(MeanHbO ~ Spatialization, paired = TRUE, p.adjust.method = "bonferroni")
#print(pwc_speech_masker)

#pwc_noise_masker <- subset(all_data_cleaned, Masker == "noise")  %>% pairwise_t_test(MeanHbO ~ Spatialization, paired = TRUE, p.adjust.method = "bonferroni")
#print(pwc_noise_masker)


#########################
# Mixed Effects model #
#########################

# R relevel factor will help reorder Spatializations if needed for contrast coding

# Simplest Model

# all_data_cleaned_collapsed <- all_data_cleaned %>% dplyr::group_by(S, Roi, Masker, Spatialization) %>%
#   dplyr::summarise(MeanHbO = mean(MeanHbO))
# 
# z1 <- mixed(MeanHbO ~ Spatialization*Roi*Masker + (1|S),
#       data= all_data_cleaned_collapsed, 
#       control = lmerControl(optimizer = "bobyqa"), 
#       method = 'LRT')
# 
# z1
# summary(z1)$AIC

# Model which includes individual channel information
all_data_cleaned_pfc <-  subset(all_data_cleaned, Roi == "pfc")
z2_pfc <- mixed(MeanHbO ~ Spatialization*Masker + (1|S) + (1|Channel),
            data= all_data_cleaned_pfc, 
            control = lmerControl(optimizer = "bobyqa"), method = 'LRT')

z2_pfc
#summary(z2_pfc)$AIC

all_data_cleaned_stg <- subset(all_data_cleaned, Roi == "stg")
z2_stg <- mixed(MeanHbO ~ Spatialization*Masker + (1|S)  + (1|Channel),
                data= all_data_cleaned_stg, 
                control = lmerControl(optimizer = "bobyqa"), method = 'LRT')

z2_stg
#summary(z2_stg)$AIC

# AIC is lower for z2 so let's use that

# ITD50 as reference

#all_data_cleaned_pfc$Masker <- relevel(all_data_cleaned_pfc$Masker, "speech")
all_data_cleaned_pfc$Spatialization <- relevel(all_data_cleaned_pfc$Spatialization, "ITD50")

pfc_lmer_itd50 <- lmer(MeanHbO ~ Spatialization + (1|S) + (1|Channel),
                data= all_data_cleaned_pfc,
                control = lmerControl(optimizer = "bobyqa"))

summary(pfc_lmer_itd50)

# ITD500 as reference

all_data_cleaned_pfc$Spatialization <- relevel(all_data_cleaned_pfc$Spatialization, "ITD500")

pfc_lmer_itd500 <- lmer(MeanHbO ~ Spatialization + (1|S) + (1|Channel),
                 data= all_data_cleaned_pfc,
                 control = lmerControl(optimizer = "bobyqa"))

summary(pfc_lmer_itd500)


# ILD70n as reference

all_data_cleaned_pfc$Spatialization <- relevel(all_data_cleaned_pfc$Spatialization, "ILD70n")

pfc_lmer_ild70n <- lmer(MeanHbO ~ Spatialization + (1|S) + (1|Channel),
                 data= all_data_cleaned_pfc,
                 control = lmerControl(optimizer = "bobyqa"))

summary(pfc_lmer_ild70n)

# ILD10 as reference

all_data_cleaned_pfc$Spatialization <- relevel(all_data_cleaned_pfc$Spatialization, "ILD10")

pfc_lmer_ild10 <- lmer(MeanHbO ~ Spatialization + (1|S) + (1|Channel),
                 data= all_data_cleaned_pfc,
                 control = lmerControl(optimizer = "bobyqa"))

summary(pfc_lmer_ild10)


pfc_se_data <- summarySE(subset(all_data_cleaned, Roi == "pfc"), measurevar="MeanHbO", groupvars=c("S","Masker","Spatialization"))
pfc_se_data <- summarySE(pfc_se_data, measurevar="MeanHbO", groupvars=c("Masker","Spatialization"))
ggplot(pfc_se_data, aes(x=factor(Spatialization, level=c('ITD50', 'ITD500', 'ILD70n', 'ILD10')), y=MeanHbO, color=Masker)) + 
  geom_errorbar(aes(ymin=MeanHbO-se, ymax=MeanHbO+se), width=.1, position=position_dodge(width=0.5)) +
  geom_line() +
  geom_point(position=position_dodge(width=0.5)) + 
  ggtitle("PFC") +
  labs(x="",y="Mean \u0394HbO") +
  ylim(0,0.12) +
  geom_signif(comparisons = list(c("ITD50","ITD500")), y_position = 0.1, tip_length = 0, color="black", annotation = c("**")) +
  geom_signif(comparisons = list(c("ITD50","ILD70n")), y_position = 0.095, tip_length = 0, color="black", annotation = c("**")) +
  geom_signif(comparisons = list(c("ITD50","ILD10")), y_position = 0.090, tip_length = 0, color="black", annotation = c("***")) +
  geom_signif(comparisons = list(c("ITD500","ILD70n")), y_position = 0.085, tip_length = 0, color="black", annotation = c("ns")) +
  geom_signif(comparisons = list(c("ITD500","ILD10")), y_position = 0.080, tip_length = 0, color="black", annotation =c("ns")) +
  geom_signif(comparisons = list(c("ILD70n","ILD10")), y_position = 0.075, tip_length = 0, color="black", annotation =c("ns")) + 
  scale_x_discrete(labels=c("ITD50" = "Small ITD", "ITD500" = "Large ITD","ILD70n" = "Natural ILD","ILD10" = "Broadband ILD"))


# ITD50 as reference

#all_data_cleaned_pfc$Masker <- relevel(all_data_cleaned_pfc$Masker, "speech")
all_data_cleaned_stg$Spatialization <- relevel(all_data_cleaned_stg$Spatialization, "ITD50")

stg_lmer_itd50 <- lmer(MeanHbO ~ Spatialization + (1|S) + (1|Channel),
                       data= all_data_cleaned_stg,
                       control = lmerControl(optimizer = "bobyqa"))

summary(stg_lmer_itd50)

# ITD500 as reference

all_data_cleaned_stg$Spatialization <- relevel(all_data_cleaned_stg$Spatialization, "ITD500")

stg_lmer_itd500 <- lmer(MeanHbO ~ Spatialization + (1|S) + (1|Channel),
                        data= all_data_cleaned_stg,
                        control = lmerControl(optimizer = "bobyqa"))

summary(stg_lmer_itd500)


# ILD70n as reference

all_data_cleaned_stg$Spatialization <- relevel(all_data_cleaned_stg$Spatialization, "ILD70n")

stg_lmer_ild70n <- lmer(MeanHbO ~ Spatialization + (1|S) + (1|Channel),
                        data= all_data_cleaned_stg,
                        control = lmerControl(optimizer = "bobyqa"))

summary(stg_lmer_ild70n)

# ILD10 as reference

all_data_cleaned_stg$Spatialization <- relevel(all_data_cleaned_stg$Spatialization, "ILD10")

stg_lmer_ild10 <- lmer(MeanHbO ~ Spatialization + (1|S) + (1|Channel),
                       data= all_data_cleaned_stg,
                       control = lmerControl(optimizer = "bobyqa"))

summary(stg_lmer_ild10)




stg_se_data <- summarySE(subset(all_data_cleaned, Roi == "stg"), measurevar="MeanHbO", groupvars=c("S","Masker","Spatialization"))
stg_se_data <- summarySE(stg_se_data, measurevar="MeanHbO", groupvars=c("Masker","Spatialization"))
ggplot(stg_se_data, aes(x=factor(Spatialization, level=c('ITD50', 'ITD500', 'ILD70n', 'ILD10')), y=MeanHbO, color=Masker)) + 
  geom_errorbar(aes(ymin=MeanHbO-se, ymax=MeanHbO+se), width=.1, position=position_dodge(width=0.5)) +
  geom_line() +
  geom_point(position=position_dodge(width=0.5)) + 
  ggtitle("STG") +
  labs(x="",y="Mean \u0394HbO", parse=TRUE) +
  ylim(0,0.12) +
  geom_signif(comparisons = list(c("ITD50","ITD500")), y_position = 0.1, tip_length = 0, color="black", annotation = c("ns")) +
  geom_signif(comparisons = list(c("ITD50","ILD70n")), y_position = 0.095, tip_length = 0, color="black", annotation = c("ns")) +
  geom_signif(comparisons = list(c("ITD50","ILD10")), y_position = 0.090, tip_length = 0, color="black", annotation = c("*")) +
  geom_signif(comparisons = list(c("ITD500","ILD70n")), y_position = 0.085, tip_length = 0, color="black", annotation = c("ns")) +
  geom_signif(comparisons = list(c("ITD500","ILD10")), y_position = 0.080, tip_length = 0, color="black", annotation =c("ns")) +
  geom_signif(comparisons = list(c("ILD70n","ILD10")), y_position = 0.075, tip_length = 0, color="black", annotation =c("**")) +
  scale_x_discrete(labels=c("ITD50" = "Small ITD", "ITD500" = "Large ITD","ILD70n" = "Natural ILD","ILD10" = "Broadband ILD"))



