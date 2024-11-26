library(tidyverse)
library(ggpubr)
library(ggplot2)
library(rstatix)
library(afex)
library(dplyr)
library(car)
require(gridExtra)

#### Data Preparation ###########
# Load in Data
speech_masker_data <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\all_subjects_uncorr_GLM_speech_masker.csv")
noise_masker_data <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\all_subjects_uncorr_GLM_noise_masker.csv")

colnames(speech_masker_data) <- c("S","Channel","ITD50","ITD500","ILD70n","ILD10")
colnames(noise_masker_data) <- c("S","Channel","ITD50","ITD500","ILD70n","ILD10")
speech_masker_data <- pivot_longer(speech_masker_data, cols=c("ITD50","ITD500","ILD70n","ILD10"), names_to = "Spatialization", values_to = "Beta")
noise_masker_data <- pivot_longer(noise_masker_data, cols=c("ITD50","ITD500","ILD70n","ILD10"), names_to = "Spatialization", values_to = "Beta")


# Combine speech and noise
speech_masker_data$Masker<-"speech"
noise_masker_data$Masker<-"noise"
all_data<-rbind(speech_masker_data,noise_masker_data)

# Add ROI information
all_data$Roi<- "NA"
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
r1<-with(all_data_cleaned, tapply(Beta, Roi, mean))

### Summary SE function ########
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  #library(plyr)
  
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


# Check for normality, remove outliers
all_data_cleaned %>% group_by(Spatialization, Masker,Roi) %>% shapiro_test(Beta)




##### PFC, Speech Masker Model #############
all_data_cleaned_pfc_speech <- subset(all_data_cleaned, Roi == "pfc" & Masker == "speech")
z2_pfc_speech <- mixed(Beta ~ Spatialization + (1|S) + (1|Channel),
                       data= all_data_cleaned_pfc_speech, 
                       control = lmerControl(optimizer = "bobyqa"), method = 'LRT')
z2_pfc_speech

# Pairwise Comparisons (treatment coding)

# ITD50 as reference
all_data_cleaned_pfc_speech$Spatialization <- relevel(all_data_cleaned_pfc_speech$Spatialization, "ITD50")
pfc_speech_lmer_itd50 <- lmer(Beta ~ Spatialization + (1|S) + (1|Channel),
                              data= all_data_cleaned_pfc_speech,
                              control = lmerControl(optimizer = "bobyqa"))#
#summary(pfc_speech_lmer_itd50)

# ITD500 as reference
all_data_cleaned_pfc_speech$Spatialization <- relevel(all_data_cleaned_pfc_speech$Spatialization, "ITD500")
pfc_speech_lmer_itd500 <- lmer(Beta ~ Spatialization + (1|S) + (1|Channel),
                               data= all_data_cleaned_pfc_speech,
                               control = lmerControl(optimizer = "bobyqa"))#
#summary(pfc_speech_lmer_itd500)

# ILD70n as reference
all_data_cleaned_pfc_speech$Spatialization <- relevel(all_data_cleaned_pfc_speech$Spatialization, "ILD70n")
pfc_speech_lmer_ild70n <- lmer(Beta ~ Spatialization + (1|S) + (1|Channel),
                               data= all_data_cleaned_pfc_speech,
                               control = lmerControl(optimizer = "bobyqa"))#
#summary(pfc_speech_lmer_ild70n)

# ILD10 as reference
all_data_cleaned_pfc_speech$Spatialization <- relevel(all_data_cleaned_pfc_speech$Spatialization, "ILD10")
pfc_speech_lmer_ild10 <- lmer(Beta ~ Spatialization + (1|S) + (1|Channel),
                              data= all_data_cleaned_pfc_speech,
                              control = lmerControl(optimizer = "bobyqa"))#
#summary(pfc_speech_lmer_ild10)



##### PFC, Noise Masker Model ########
all_data_cleaned_pfc_noise <- subset(all_data_cleaned, Roi == "pfc" & Masker == "noise")
z2_pfc_noise <- mixed(Beta ~ Spatialization + (1|S) + (1|Channel),
                      data= all_data_cleaned_pfc_noise, 
                      control = lmerControl(optimizer = "bobyqa"), method = 'LRT')
z2_pfc_noise

# Pairwise Comparisons (treatment coding)

# ITD50 as reference
all_data_cleaned_pfc_noise$Spatialization <- relevel(all_data_cleaned_pfc_noise$Spatialization, "ITD50")
pfc_noise_lmer_itd50 <- lmer(Beta ~ Spatialization + (1|S) + (1|Channel),
                             data= all_data_cleaned_pfc_noise,
                             control = lmerControl(optimizer = "bobyqa"))#
#summary(pfc_noise_lmer_itd50)

# ITD500 as reference
all_data_cleaned_pfc_noise$Spatialization <- relevel(all_data_cleaned_pfc_noise$Spatialization, "ITD500")
pfc_noise_lmer_itd500 <- lmer(Beta ~ Spatialization + (1|S) + (1|Channel),
                              data= all_data_cleaned_pfc_noise,
                              control = lmerControl(optimizer = "bobyqa"))#
#summary(pfc_noise_lmer_itd500)

# ILD70n as reference
all_data_cleaned_pfc_noise$Spatialization <- relevel(all_data_cleaned_pfc_noise$Spatialization, "ILD70n")
pfc_noise_lmer_ild70n <- lmer(Beta ~ Spatialization + (1|S) + (1|Channel),
                              data= all_data_cleaned_pfc_noise,
                              control = lmerControl(optimizer = "bobyqa"))#
#summary(pfc_noise_lmer_ild70n)

# ILD10 as reference
all_data_cleaned_pfc_noise$Spatialization <- relevel(all_data_cleaned_pfc_noise$Spatialization, "ILD10")
pfc_noise_lmer_ild10 <- lmer(Beta ~ Spatialization + (1|S) + (1|Channel),
                             data= all_data_cleaned_pfc_noise,
                             control = lmerControl(optimizer = "bobyqa"))#
#summary(pfc_lmer_noise_ild10)


##### stg, Speech Masker Model #############
all_data_cleaned_stg_speech <- subset(all_data_cleaned, Roi == "stg" & Masker == "speech")
z2_stg_speech <- mixed(Beta ~ Spatialization + (1|S) + (1|Channel),
                       data= all_data_cleaned_stg_speech, 
                       control = lmerControl(optimizer = "bobyqa"), method = 'LRT')
z2_stg_speech

# Pairwise Comparisons (treatment coding)

# ITD50 as reference
all_data_cleaned_stg_speech$Spatialization <- relevel(all_data_cleaned_stg_speech$Spatialization, "ITD50")
stg_speech_lmer_itd50 <- lmer(Beta ~ Spatialization + (1|S) + (1|Channel),
                              data= all_data_cleaned_stg_speech,
                              control = lmerControl(optimizer = "bobyqa"))#
#summary(stg_speech_lmer_itd50)

# ITD500 as reference
all_data_cleaned_stg_speech$Spatialization <- relevel(all_data_cleaned_stg_speech$Spatialization, "ITD500")
stg_speech_lmer_itd500 <- lmer(Beta ~ Spatialization + (1|S) + (1|Channel),
                               data= all_data_cleaned_stg_speech,
                               control = lmerControl(optimizer = "bobyqa"))#
#summary(stg_speech_lmer_itd500)

# ILD70n as reference
all_data_cleaned_stg_speech$Spatialization <- relevel(all_data_cleaned_stg_speech$Spatialization, "ILD70n")
stg_speech_lmer_ild70n <- lmer(Beta ~ Spatialization + (1|S) + (1|Channel),
                               data= all_data_cleaned_stg_speech,
                               control = lmerControl(optimizer = "bobyqa"))#
#summary(stg_speech_lmer_ild70n)

# ILD10 as reference
all_data_cleaned_stg_speech$Spatialization <- relevel(all_data_cleaned_stg_speech$Spatialization, "ILD10")
stg_speech_lmer_ild10 <- lmer(Beta ~ Spatialization + (1|S) + (1|Channel),
                              data= all_data_cleaned_stg_speech,
                              control = lmerControl(optimizer = "bobyqa"))#
#summary(stg_speech_lmer_ild10)


##### stg, Noise Masker Model #####
all_data_cleaned_stg_noise <- subset(all_data_cleaned, Roi == "stg" & Masker == "noise")
z2_stg_noise <- mixed(Beta ~ Spatialization + (1|S) + (1|Channel),
                      data= all_data_cleaned_stg_noise, 
                      control = lmerControl(optimizer = "bobyqa"), method = 'LRT')
z2_stg_noise

# Pairwise Comparisons (treatment coding)

# ITD50 as reference
all_data_cleaned_stg_noise$Spatialization <- relevel(all_data_cleaned_stg_noise$Spatialization, "ITD50")
stg_noise_lmer_itd50 <- lmer(Beta ~ Spatialization + (1|S) + (1|Channel),
                             data= all_data_cleaned_stg_noise,
                             control = lmerControl(optimizer = "bobyqa"))#
#summary(stg_noise_lmer_itd50)

# ITD500 as reference
all_data_cleaned_stg_noise$Spatialization <- relevel(all_data_cleaned_stg_noise$Spatialization, "ITD500")
stg_noise_lmer_itd500 <- lmer(Beta ~ Spatialization + (1|S) + (1|Channel),
                              data= all_data_cleaned_stg_noise,
                              control = lmerControl(optimizer = "bobyqa"))#
#summary(stg_noise_lmer_itd500)

# ILD70n as reference
all_data_cleaned_stg_noise$Spatialization <- relevel(all_data_cleaned_stg_noise$Spatialization, "ILD70n")
stg_noise_lmer_ild70n <- lmer(Beta ~ Spatialization + (1|S) + (1|Channel),
                              data= all_data_cleaned_stg_noise,
                              control = lmerControl(optimizer = "bobyqa"))#
#summary(stg_noise_lmer_ild70n)

# ILD10 as reference
all_data_cleaned_stg_noise$Spatialization <- relevel(all_data_cleaned_stg_noise$Spatialization, "ILD10")
stg_noise_lmer_ild10 <- lmer(Beta ~ Spatialization + (1|S) + (1|Channel),
                             data= all_data_cleaned_stg_noise,
                             control = lmerControl(optimizer = "bobyqa"))#
#summary(stg_lmer_noise_ild10)




##### PFC Plot ##########
pfc_se_data_speech <- summarySE(all_data_cleaned_pfc_speech, measurevar="Beta", groupvars=c("S","Masker","Spatialization"), na.rm = TRUE)
pfc_se_data_speech <- summarySE(pfc_se_data_speech, measurevar="Beta", groupvars=c("Masker","Spatialization"), na.rm = TRUE)
plotspeech <- ggplot(pfc_se_data_speech, aes(x=factor(Spatialization, level=c('ITD50', 'ITD500', 'ILD70n', 'ILD10')), y=Beta, color = Spatialization)) + 
  scale_color_manual(values = c("ITD50" = "red","ITD500" =  "blue","ILD70n" = "magenta","ILD10" = "green")) +
  geom_errorbar(aes(ymin=Beta-se, ymax=Beta+se, color=Spatialization), width=.1, position=position_dodge(width=0.5)) +
  geom_point(aes(color = Spatialization),size = 4, position=position_dodge(width=0.5)) + 
  ggtitle("Speech Masker PFC") +
  labs(x="",y="Mean \u0394HbO (\u03BCM)") +
  ylim(-0.25,1.0) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_text(size=12)) +
  scale_x_discrete(labels=c("ITD50" = "Small\nITD", "ITD500" = "Large\nITD","ILD70n" = "Natural\nILD","ILD10" = "Broadband\nILD")) +
  theme(legend.position="none")
#geom_signif(comparisons = list(c("ITD50","ITD500")), y_position = 0.09, tip_length = 0, color="black", annotation = c("***"), textsize = 5) +
#geom_signif(comparisons = list(c("ITD50","ILD70n")), y_position = 0.10, tip_length = 0, color="black", annotation = c("***"), textsize = 5) +
#geom_signif(comparisons = list(c("ITD50","ILD10")), y_position = 0.11, tip_length = 0, color="black", annotation = c("***"), textsize = 5)


pfc_se_data_noise <- summarySE(all_data_cleaned_pfc_noise, measurevar="Beta", groupvars=c("S","Masker","Spatialization"), na.rm = TRUE)
pfc_se_data_noise <- summarySE(pfc_se_data_noise, measurevar="Beta", groupvars=c("Masker","Spatialization"), na.rm = TRUE)
plotnoise <- ggplot(pfc_se_data_noise, aes(x=factor(Spatialization, level=c('ITD50', 'ITD500', 'ILD70n', 'ILD10')), y=Beta, color = Spatialization)) + 
  scale_color_manual(values = c("ITD50" = "red","ITD500" =  "blue","ILD70n" = "magenta","ILD10" = "green")) +
  geom_errorbar(aes(ymin=Beta-se, ymax=Beta+se, color=Spatialization), width=.1, position=position_dodge(width=0.5)) +
  geom_point(aes(color = Spatialization),size = 4, position=position_dodge(width=0.5)) + 
  ggtitle("Noise Masker PFC") +
  labs(x="",y="") +
  ylim(-0.25,1.0) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_blank()) +
  scale_x_discrete(labels=c("ITD50" = "Small\nITD", "ITD500" = "Large\nITD","ILD70n" = "Natural\nILD","ILD10" = "Broadband\nILD")) +
  theme(legend.position="none")

grid.arrange(plotspeech,plotnoise, ncol=2, widths = c(1,0.9))




######### STG Plot ###############

stg_se_data_speech <- summarySE(all_data_cleaned_stg_speech, measurevar="Beta", groupvars=c("S","Spatialization"), na.rm=TRUE)
stg_se_data_speech <- summarySE(stg_se_data_speech, measurevar="Beta", groupvars=c("Spatialization"),  na.rm=TRUE)
plotspeech <- ggplot(stg_se_data_speech, aes(x=factor(Spatialization, level=c('ITD50', 'ITD500', 'ILD70n', 'ILD10')), y=Beta, color = Spatialization)) + 
  scale_color_manual(values = c("ITD50" = "red","ITD500" =  "blue","ILD70n" = "magenta","ILD10" = "green")) +
  geom_errorbar(aes(ymin=Beta-se, ymax=Beta+se, color=Spatialization), width=.1, position=position_dodge(width=0.5)) +
  geom_point(aes(color = Spatialization),size = 4, position=position_dodge(width=0.5)) +  
  ggtitle("Speech Masker STG") +
  labs(x="",y="Mean \u0394HbO (\u03BCM)", parse=TRUE) +
  ylim(-0.25,1.0) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_text(size=12)) +
  scale_x_discrete(labels=c("ITD50" = "Small\nITD", "ITD500" = "Large\nITD","ILD70n" = "Natural\nILD","ILD10" = "Broadband\nILD")) +
  theme(legend.position="none")
  #geom_signif(comparisons = list(c("ITD50","ILD10")), y_position = 0.090, tip_length = 0, color="black", annotation = c("**"), textsize = 5) +
  #geom_signif(comparisons = list(c("ITD500","ILD10")), y_position = 0.080, tip_length = 0, color="black", annotation =c("*"), textsize = 5) +
  #geom_signif(comparisons = list(c("ILD70n","ILD10")), y_position = 0.070, tip_length = 0, color="black", annotation =c("**"), textsize = 5) 


stg_se_data_noise <- summarySE(all_data_cleaned_stg_noise, measurevar="Beta", groupvars=c("S","Spatialization"), na.rm=TRUE)
stg_se_data_noise <- summarySE(stg_se_data_noise, measurevar="Beta", groupvars=c("Spatialization"),  na.rm=TRUE)
plotnoise <- ggplot(stg_se_data_noise, aes(x=factor(Spatialization, level=c('ITD50', 'ITD500', 'ILD70n', 'ILD10')), y=Beta, color = Spatialization)) + 
  scale_color_manual(values = c("ITD50" = "red","ITD500" =  "blue","ILD70n" = "magenta","ILD10" = "green")) +
  geom_errorbar(aes(ymin=Beta-se, ymax=Beta+se, color=Spatialization), width=.1, position=position_dodge(width=0.5)) +
  geom_point(aes(color = Spatialization),size = 4, position=position_dodge(width=0.5)) +  
  ggtitle("Noise Masker STG") +
  labs(x="", y="") +
  ylim(-0.25,1.0) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_blank()) +
  scale_x_discrete(labels=c("ITD50" = "Small\nITD", "ITD500" = "Large\nITD","ILD70n" = "Natural\nILD","ILD10" = "Broadband\nILD")) +
  theme(legend.position="none") 

grid.arrange(plotspeech,plotnoise, ncol=2, widths = c(1,0.9))


