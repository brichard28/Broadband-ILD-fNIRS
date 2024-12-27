

library(tidyverse)
library(ggpubr)
library(ggplot2)
library(rstatix)
library(afex)
library(dplyr)
library(car)
require(gridExtra)

#### Data Preparation HBO ###########
# Load in Data
speech_masker_data_hbo <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\all_subjects_mean_during_stim_speech_masker.csv")
noise_masker_data_hbo <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\all_subjects_mean_during_stim_noise_masker.csv")

colnames(speech_masker_data_hbo) <- c("S","Channel","ITD50","ITD500","ILD70n","ILD10")
colnames(noise_masker_data_hbo) <- c("S","Channel","ITD50","ITD500","ILD70n","ILD10")
speech_masker_data_hbo <- pivot_longer(speech_masker_data_hbo, cols=c("ITD50","ITD500","ILD70n","ILD10"), names_to = "Spatialization", values_to = "MeanHb")
noise_masker_data_hbo <- pivot_longer(noise_masker_data_hbo, cols=c("ITD50","ITD500","ILD70n","ILD10"), names_to = "Spatialization", values_to = "MeanHb")


# Combine speech and noise
speech_masker_data_hbo$Masker<-"speech"
noise_masker_data_hbo$Masker<-"noise"
all_data_hbo<-rbind(speech_masker_data_hbo,noise_masker_data_hbo)

# Add ROI information
all_data_hbo$Roi<- "NA"
pfc_channels <- c(0,1,2,3,4,5)
stg_channels <- c(6,7,8,9,10,11,12,13)
all_data_hbo$Roi[which(all_data_hbo$Channel %in% pfc_channels)] <- "pfc"
all_data_hbo$Roi[which(all_data_hbo$Channel %in% stg_channels)] <- "stg"

# Organize Factors
to.factor <- c('S', 'Roi','Spatialization')
all_data_hbo[, to.factor] <- lapply(all_data_hbo[, to.factor], as.factor)

all_data_hbo$Masker <- as.factor(all_data_hbo$Masker)
all_data_cleaned_hbo <- na.omit(all_data_hbo)

all_data_cleaned_hbo %>% group_by(Spatialization, Masker,Roi) %>% shapiro_test(MeanHb)
all_data_cleaned_pfc_speech_hbo <- subset(all_data_cleaned_hbo, Roi == "pfc" & Masker == "speech")
all_data_cleaned_pfc_noise_hbo <- subset(all_data_cleaned_hbo, Roi == "pfc" & Masker == "noise")
all_data_cleaned_stg_speech_hbo <- subset(all_data_cleaned_hbo, Roi == "stg" & Masker == "speech")
all_data_cleaned_stg_noise_hbo <- subset(all_data_cleaned_hbo, Roi == "stg" & Masker == "noise")

#### Data Preparaion HBR #####
# Load in Data
speech_masker_data_hbr <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\all_subjects_mean_during_stim_speech_masker_hbr.csv")
noise_masker_data_hbr <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\all_subjects_mean_during_stim_noise_masker_hbr.csv")

colnames(speech_masker_data_hbr) <- c("S","Channel","ITD50","ITD500","ILD70n","ILD10")
colnames(noise_masker_data_hbr) <- c("S","Channel","ITD50","ITD500","ILD70n","ILD10")
speech_masker_data_hbr <- pivot_longer(speech_masker_data_hbr, cols=c("ITD50","ITD500","ILD70n","ILD10"), names_to = "Spatialization", values_to = "MeanHb")
noise_masker_data_hbr <- pivot_longer(noise_masker_data_hbr, cols=c("ITD50","ITD500","ILD70n","ILD10"), names_to = "Spatialization", values_to = "MeanHb")


# Combine speech and noise
speech_masker_data_hbr$Masker<-"speech"
noise_masker_data_hbr$Masker<-"noise"
all_data_hbr<-rbind(speech_masker_data_hbr,noise_masker_data_hbr)

# Add ROI information
all_data_hbr$Roi<-NA
pfc_channels <- c(0,1,2,3,4,5)
stg_channels <- c(6,7,8,9,10,11,12,13)
all_data_hbr$Roi[which(all_data_hbr$Channel %in% pfc_channels)] <- "pfc"
all_data_hbr$Roi[which(all_data_hbr$Channel %in% stg_channels)] <- "stg"

# Change ROI information

#all_data$Roi[all_data$Roi == "0"] <- "pfc"
#all_data$Roi[all_data$Roi == "1"] <- "stg"

# Organize Factors
to.factor <- c('S', 'Roi','Spatialization')
all_data_hbr[, to.factor] <- lapply(all_data_hbr[, to.factor], as.factor)

all_data_hbr$Masker <- as.factor(all_data_hbr$Masker)
all_data_cleaned_hbr <- na.omit(all_data_hbr)

all_data_cleaned_hbr %>% group_by(Spatialization, Masker,Roi) %>% shapiro_test(MeanHb)
all_data_cleaned_pfc_speech_hbr <- subset(all_data_cleaned_hbr, Roi == "pfc" & Masker == "speech")
all_data_cleaned_pfc_noise_hbr <- subset(all_data_cleaned_hbr, Roi == "pfc" & Masker == "noise")
all_data_cleaned_stg_speech_hbr <- subset(all_data_cleaned_hbr, Roi == "stg" & Masker == "speech")
all_data_cleaned_stg_noise_hbr <- subset(all_data_cleaned_hbr, Roi == "stg" & Masker == "noise")


##### Combining Data #####
all_data_cleaned_hbo$chromophore <- "HbO"
all_data_cleaned_hbr$chromophore <- "HbR"
all_data <- rbind(all_data_cleaned_hbo,all_data_cleaned_hbr)
all_data_pfc_speech <- subset(all_data, Roi == "pfc" & Masker == "speech")
all_data_pfc_noise <- subset(all_data, Roi == "pfc" & Masker == "noise")
all_data_stg_speech <- subset(all_data, Roi == "stg" & Masker == "speech")
all_data_stg_noise <- subset(all_data, Roi == "stg" & Masker == "noise")

### Summary SE function ########
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


# Check for normality, remove outliers
all_data_cleaned_hbo %>% group_by(Spatialization, Masker,Roi) %>% shapiro_test(MeanHb)



##### PFC Plot HbO ##########

# Speech
pfc_se_data_speech_all <- summarySE(all_data_pfc_speech, measurevar="MeanHb", groupvars=c("S","Masker","Spatialization","chromophore"), na.rm = TRUE)
pfc_se_data_speech_all <- summarySE(pfc_se_data_speech_all, measurevar="MeanHb", groupvars=c("Masker","Spatialization","chromophore"), na.rm = TRUE)

plotspeech <- ggplot(subset(pfc_se_data_speech_all, chromophore == "HbO"), aes(x=factor(Spatialization, level=c('ITD50', 'ITD500', 'ILD70n', 'ILD10')), y=MeanHb, shape = Spatialization, fill = "black")) + 
  scale_shape_manual(values = c("ITD50"=21, "ITD500" = 24, "ILD70n" = 22, "ILD10" = 25)) +
  scale_fill_manual(values = "black") +
  scale_color_manual(values = "black") +
  geom_errorbar(aes(ymin=MeanHb-se, ymax=MeanHb+se), width=.1, position=position_dodge(width=0.5)) +
  geom_point(aes(),size = 4, position=position_dodge(width=0.5)) +
  ggtitle("Speech Masker PFC") +
  labs(x="",y="Mean \u0394HbO (\u03BCM)") +
  ylim(0,0.12) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_text(size=12)) +
  scale_x_discrete(labels=c("ITD50" = "Small\nITD", "ITD500" = "Large\nITD","ILD70n" = "Natural\nILD","ILD10" = "Broadband\nILD")) +
  theme(legend.position="none") +
  geom_signif(comparisons = list(c("ITD50","ITD500")), y_position = 0.09, tip_length = 0, color="black", annotation = c("***"), textsize = 5) +
  geom_signif(comparisons = list(c("ITD50","ILD70n")), y_position = 0.10, tip_length = 0, color="black", annotation = c("***"), textsize = 5) +
  geom_signif(comparisons = list(c("ITD50","ILD10")), y_position = 0.11, tip_length = 0, color="black", annotation = c("***"), textsize = 5)

# Noise HbO and HbR

pfc_se_data_noise_all <- summarySE(all_data_pfc_noise, measurevar="MeanHb", groupvars=c("S","Masker","Spatialization","chromophore"), na.rm = TRUE)
pfc_se_data_noise_all <- summarySE(pfc_se_data_noise_all, measurevar="MeanHb", groupvars=c("Masker","Spatialization","chromophore"), na.rm = TRUE)

plotnoise <- ggplot(subset(pfc_se_data_noise_all, chromophore == "HbO"), aes(x=factor(Spatialization, level=c('ITD50', 'ITD500', 'ILD70n', 'ILD10')), y=MeanHb, shape = Spatialization, fill = "black")) + 
  scale_shape_manual(values = c("ITD50"=21, "ITD500" = 24, "ILD70n" = 22, "ILD10" = 25)) +
  scale_fill_manual(values = "black") +
  scale_color_manual(values = "black") +
  geom_errorbar(aes(ymin=MeanHb-se, ymax=MeanHb+se), width=.1, position=position_dodge(width=0.5)) +
  geom_point(aes(),size = 4, position=position_dodge(width=0.5)) +
  ggtitle("Noise Masker PFC") +
  labs(x="",y="") +
  ylim(0,0.12) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_blank()) +
  scale_x_discrete(labels=c("ITD50" = "Small\nITD", "ITD500" = "Large\nITD","ILD70n" = "Natural\nILD","ILD10" = "Broadband\nILD")) +
  theme(legend.position="none")

grid.arrange(plotspeech,plotnoise, ncol=2, widths = c(1,0.9))



######### STG Plot ###############

stg_se_data_speech_all <- summarySE(all_data_stg_speech, measurevar="MeanHb", groupvars=c("S","Spatialization"), na.rm=TRUE)
stg_se_data_speech_all <- summarySE(stg_se_data_speech_all, measurevar="MeanHb", groupvars=c("Spatialization"),  na.rm=TRUE)
plotspeech <- ggplot(subset(stg_se_data_speech, chromophore == "HbO"), aes(x=factor(Spatialization, level=c('ITD50', 'ITD500', 'ILD70n', 'ILD10')), y=MeanHb, color = Spatialization)) + 
  scale_shape_manual(values = c("ITD50"=21, "ITD500" = 24, "ILD70n" = 22, "ILD10" = 25)) +
  scale_fill_manual(values = "black") +
  scale_color_manual(values = "black") +
  geom_errorbar(aes(ymin=MeanHb-se, ymax=MeanHb+se, color=Spatialization), width=.1, position=position_dodge(width=0.5)) +
  geom_point(aes(color = Spatialization),size = 4, position=position_dodge(width=0.5)) +  
  ggtitle("Speech Masker STG") +
  labs(x="",y="Mean \u0394HbO (\u03BCM)", parse=TRUE) +
  ylim(0,0.12) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_text(size=12)) +
  scale_x_discrete(labels=c("ITD50" = "Small\nITD", "ITD500" = "Large\nITD","ILD70n" = "Natural\nILD","ILD10" = "Broadband\nILD")) +
  theme(legend.position="none")
  #geom_signif(comparisons = list(c("ITD50","ILD10")), y_position = 0.090, tip_length = 0, color="black", annotation = c("**"), textsize = 5) +
  #geom_signif(comparisons = list(c("ITD500","ILD10")), y_position = 0.080, tip_length = 0, color="black", annotation =c("*"), textsize = 5) +
  #geom_signif(comparisons = list(c("ILD70n","ILD10")), y_position = 0.070, tip_length = 0, color="black", annotation =c("**"), textsize = 5) 


stg_se_data_noise_all <- summarySE(all_data_stg_noise, measurevar="MeanHb", groupvars=c("S","Spatialization"), na.rm=TRUE)
stg_se_data_noise_all <- summarySE(stg_se_data_noise_all, measurevar="MeanHb", groupvars=c("Spatialization"),  na.rm=TRUE)
plotnoise <- ggplot(subset(stg_se_data_noise_all, chromophore == "HbO"), aes(x=factor(Spatialization, level=c('ITD50', 'ITD500', 'ILD70n', 'ILD10')), y=MeanHb, color = Spatialization)) + 
  scale_shape_manual(values = c("ITD50"=21, "ITD500" = 24, "ILD70n" = 22, "ILD10" = 25)) +
  scale_fill_manual(values = "black") +
  scale_color_manual(values = "black") +
  geom_errorbar(aes(ymin=MeanHb-se, ymax=MeanHb+se, color=Spatialization), width=.1, position=position_dodge(width=0.5)) +
  geom_point(aes(color = Spatialization),size = 4, position=position_dodge(width=0.5)) +  
  ggtitle("Noise Masker STG") +
  labs(x="", y="") +
  ylim(0,0.12) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_blank()) +
  scale_x_discrete(labels=c("ITD50" = "Small\nITD", "ITD500" = "Large\nITD","ILD70n" = "Natural\nILD","ILD10" = "Broadband\nILD")) +
  theme(legend.position="none") 

grid.arrange(plotspeech,plotnoise, ncol=2, widths = c(1,0.9))




##### STAT MODELS HBO #####

# PFC, Speech Masker Model #
model_pfc_speech_hbo <- mixed(MeanHb ~ Spatialization + (1|S) + (1|Channel),
                              data= all_data_cleaned_pfc_speech_hbo, 
                              control = lmerControl(optimizer = "bobyqa"), method = 'LRT')
model_pfc_speech_hbo 

# Pairwise Comparisons (treatment coding)

# ITD50 as reference
all_data_cleaned_pfc_speech_hbo$Spatialization <- relevel(all_data_cleaned_pfc_speech_hbo$Spatialization, "ITD50")
pfc_speech_lmer_itd50 <- lmer(MeanHb ~ Spatialization + (1|S) + (1|Channel),
                              data= all_data_cleaned_pfc_speech_hbo,
                              control = lmerControl(optimizer = "bobyqa"))#
#summary(pfc_speech_lmer_itd50)

# ITD500 as reference
all_data_cleaned_pfc_speech_hbo$Spatialization <- relevel(all_data_cleaned_pfc_speech_hbo$Spatialization, "ITD500")
pfc_speech_lmer_itd500 <- lmer(MeanHb ~ Spatialization + (1|S) + (1|Channel),
                               data= all_data_cleaned_pfc_speech_hbo,
                               control = lmerControl(optimizer = "bobyqa"))#
#summary(pfc_speech_lmer_itd500)

# ILD70n as reference
all_data_cleaned_pfc_speech_hbo$Spatialization <- relevel(all_data_cleaned_pfc_speech_hbo$Spatialization, "ILD70n")
pfc_speech_lmer_ild70n <- lmer(MeanHb ~ Spatialization + (1|S) + (1|Channel),
                               data= all_data_cleaned_pfc_speech_hbo,
                               control = lmerControl(optimizer = "bobyqa"))#
#summary(pfc_speech_lmer_ild70n)

# ILD10 as reference
all_data_cleaned_pfc_speech_hbo$Spatialization <- relevel(all_data_cleaned_pfc_speech_hbo$Spatialization, "ILD10")
pfc_speech_lmer_ild10 <- lmer(MeanHb ~ Spatialization + (1|S) + (1|Channel),
                              data= all_data_cleaned_pfc_speech_hbo,
                              control = lmerControl(optimizer = "bobyqa"))#
#summary(pfc_speech_lmer_ild10)



#PFC, Noise Masker Model #
model_pfc_noise_hbo <- mixed(MeanHb ~ Spatialization + (1|S) + (1|Channel),
                             data= all_data_cleaned_pfc_noise_hbo, 
                             control = lmerControl(optimizer = "bobyqa"), method = 'LRT')
model_pfc_noise_hbo

# Pairwise Comparisons (treatment coding)

# ITD50 as reference
all_data_cleaned_pfc_noise_hbo$Spatialization <- relevel(all_data_cleaned_pfc_noise_hbo$Spatialization, "ITD50")
pfc_noise_lmer_itd50 <- lmer(MeanHb ~ Spatialization + (1|S) + (1|Channel),
                             data= all_data_cleaned_pfc_noise_hbo,
                             control = lmerControl(optimizer = "bobyqa"))#
#summary(pfc_noise_lmer_itd50)

# ITD500 as reference
all_data_cleaned_pfc_noise_hbo$Spatialization <- relevel(all_data_cleaned_pfc_noise_hbo$Spatialization, "ITD500")
pfc_noise_lmer_itd500 <- lmer(MeanHb ~ Spatialization + (1|S) + (1|Channel),
                              data= all_data_cleaned_pfc_noise_hbo,
                              control = lmerControl(optimizer = "bobyqa"))#
#summary(pfc_noise_lmer_itd500)

# ILD70n as reference
all_data_cleaned_pfc_noise_hbo$Spatialization <- relevel(all_data_cleaned_pfc_noise_hbo$Spatialization, "ILD70n")
pfc_noise_lmer_ild70n <- lmer(MeanHb ~ Spatialization + (1|S) + (1|Channel),
                              data= all_data_cleaned_pfc_noise_hbo,
                              control = lmerControl(optimizer = "bobyqa"))#
#summary(pfc_noise_lmer_ild70n)

# ILD10 as reference
all_data_cleaned_pfc_noise_hbo$Spatialization <- relevel(all_data_cleaned_pfc_noise_hbo$Spatialization, "ILD10")
pfc_noise_lmer_ild10 <- lmer(MeanHb ~ Spatialization + (1|S) + (1|Channel),
                             data= all_data_cleaned_pfc_noise_hbo,
                             control = lmerControl(optimizer = "bobyqa"))#
#summary(pfc_lmer_noise_ild10)


# STG, Speech Masker Model #
model_stg_speech_hbo <- mixed(MeanHb ~ Spatialization + (1|S) + (1|Channel),
                              data= all_data_cleaned_stg_speech_hbo, 
                              control = lmerControl(optimizer = "bobyqa"), method = 'LRT')
model_stg_speech_hbo

# Pairwise Comparisons (treatment coding)

# ITD50 as reference
all_data_cleaned_stg_speech_hbo$Spatialization <- relevel(all_data_cleaned_stg_speech_hbo$Spatialization, "ITD50")
stg_speech_lmer_itd50 <- lmer(MeanHb ~ Spatialization + (1|S) + (1|Channel),
                              data= all_data_cleaned_stg_speech_hbo,
                              control = lmerControl(optimizer = "bobyqa"))#
#summary(stg_speech_lmer_itd50)

# ITD500 as reference
all_data_cleaned_stg_speech_hbo$Spatialization <- relevel(all_data_cleaned_stg_speech_hbo$Spatialization, "ITD500")
stg_speech_lmer_itd500 <- lmer(MeanHb ~ Spatialization + (1|S) + (1|Channel),
                               data= all_data_cleaned_stg_speech_hbo,
                               control = lmerControl(optimizer = "bobyqa"))#
#summary(stg_speech_lmer_itd500)

# ILD70n as reference
all_data_cleaned_stg_speech_hbo$Spatialization <- relevel(all_data_cleaned_stg_speech_hbo$Spatialization, "ILD70n")
stg_speech_lmer_ild70n <- lmer(MeanHb ~ Spatialization + (1|S) + (1|Channel),
                               data= all_data_cleaned_stg_speech_hbo,
                               control = lmerControl(optimizer = "bobyqa"))#
#summary(stg_speech_lmer_ild70n)

# ILD10 as reference
all_data_cleaned_stg_speech_hbo$Spatialization <- relevel(all_data_cleaned_stg_speech_hbo$Spatialization, "ILD10")
stg_speech_lmer_ild10 <- lmer(MeanHb ~ Spatialization + (1|S) + (1|Channel),
                              data= all_data_cleaned_stg_speech_hbo,
                              control = lmerControl(optimizer = "bobyqa"))#
#summary(stg_speech_lmer_ild10)


# STG, Noise Masker Model #
model_stg_noise_hbo <- mixed(MeanHb ~ Spatialization + (1|S) + (1|Channel),
                             data= all_data_cleaned_stg_noise_hbo, 
                             control = lmerControl(optimizer = "bobyqa"), method = 'LRT')
model_stg_noise_hbo



#### STAT MODELS HBR ####

# PFC Speech
model_pfc_speech_hbr <- mixed(MeanHb ~ Spatialization + (1|S) + (1|Channel),
                              data= all_data_cleaned_pfc_speech_hbr, 
                              control = lmerControl(optimizer = "bobyqa"), method = 'LRT')
model_pfc_speech_hbr


# PFC Noise
model_pfc_noise_hbr <- mixed(MeanHb ~ Spatialization + (1|S) + (1|Channel),
                             data= all_data_cleaned_pfc_noise_hbr, 
                             control = lmerControl(optimizer = "bobyqa"), method = 'LRT')
model_pfc_noise_hbr

# STG Speech
model_stg_speech_hbr <- mixed(MeanHb ~ Spatialization + (1|S) + (1|Channel),
                              data= all_data_cleaned_stg_speech_hbr, 
                              control = lmerControl(optimizer = "bobyqa"), method = 'LRT')
model_stg_speech_hbr

# STG Noise
model_stg_noise_hbr <- mixed(MeanHb ~ Spatialization + (1|S) + (1|Channel),
                             data= all_data_cleaned_stg_noise_hbr, 
                             control = lmerControl(optimizer = "bobyqa"), method = 'LRT')
model_stg_noise_hbr

