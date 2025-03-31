
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(rstatix)
library(afex)
library(dplyr)
library(car)

# Data Preparation ####
# Load in Data
attend_right_speech_data_hbo <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\all_subjects_mean_during_stim_lateralization_target_right_speech_masker.csv")
attend_left_speech_data_hbo <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\all_subjects_mean_during_stim_lateralization_target_left_speech_masker.csv")
attend_right_noise_data_hbo <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\all_subjects_mean_during_stim_lateralization_target_right_noise_masker.csv")
attend_left_noise_data_hbo <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\all_subjects_mean_during_stim_lateralization_target_left_noise_masker.csv")

attend_right_speech_data_hbr <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\all_subjects_mean_during_stim_lateralization_target_right_speech_masker_hbr.csv")
attend_left_speech_data_hbr <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\all_subjects_mean_during_stim_lateralization_target_left_speech_masker_hbr.csv")
attend_right_noise_data_hbr <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\all_subjects_mean_during_stim_lateralization_target_right_noise_masker_hbr.csv")
attend_left_noise_data_hbr <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\all_subjects_mean_during_stim_lateralization_target_left_noise_masker_hbr.csv")

colnames(attend_right_speech_data_hbo) <- c("S","Channel","ITD50","ITD500","ILD70n","ILD10")
colnames(attend_right_noise_data_hbo) <- c("S","Channel","ITD50","ITD500","ILD70n","ILD10")
colnames(attend_left_speech_data_hbo) <- c("S","Channel","ITD50","ITD500","ILD70n","ILD10")
colnames(attend_left_noise_data_hbo) <- c("S","Channel","ITD50","ITD500","ILD70n","ILD10")

colnames(attend_right_speech_data_hbr) <- c("S","Channel","ITD50","ITD500","ILD70n","ILD10")
colnames(attend_right_noise_data_hbr) <- c("S","Channel","ITD50","ITD500","ILD70n","ILD10")
colnames(attend_left_speech_data_hbr) <- c("S","Channel","ITD50","ITD500","ILD70n","ILD10")
colnames(attend_left_noise_data_hbr) <- c("S","Channel","ITD50","ITD500","ILD70n","ILD10")

attend_right_speech_data_hbo$chromophore <- "HbO"
attend_right_speech_data_hbr$chromophore <- "HbR"
attend_left_speech_data_hbo$chromophore <- "HbO"
attend_left_speech_data_hbr$chromophore <- "HbR"

attend_right_noise_data_hbo$chromophore <- "HbO"
attend_right_noise_data_hbr$chromophore <- "HbR"
attend_left_noise_data_hbo$chromophore <- "HbO"
attend_left_noise_data_hbr$chromophore <- "HbR"


attend_right_speech_data <- rbind(attend_right_speech_data_hbo,attend_right_speech_data_hbr)
attend_left_speech_data <- rbind(attend_left_speech_data_hbo,attend_left_speech_data_hbr)
attend_right_noise_data <- rbind(attend_right_noise_data_hbo,attend_right_noise_data_hbr)
attend_left_noise_data <- rbind(attend_left_noise_data_hbo,attend_left_noise_data_hbr)

attend_right_speech_data <- pivot_longer(attend_right_speech_data, cols=c("ITD50","ITD500","ILD70n","ILD10"), names_to = "Spatialization", values_to = "MeanHb")
attend_right_noise_data <- pivot_longer(attend_right_noise_data, cols=c("ITD50","ITD500","ILD70n","ILD10"), names_to = "Spatialization", values_to = "MeanHb")
attend_left_speech_data <- pivot_longer(attend_left_speech_data, cols=c("ITD50","ITD500","ILD70n","ILD10"), names_to = "Spatialization", values_to = "MeanHb")
attend_left_noise_data <- pivot_longer(attend_left_noise_data, cols=c("ITD50","ITD500","ILD70n","ILD10"), names_to = "Spatialization", values_to = "MeanHb")

attend_right_speech_data <- attend_right_speech_data[!attend_right_speech_data$S == '2',]
attend_right_speech_data <- attend_right_speech_data[!attend_right_speech_data$S == '1',]
attend_right_noise_data <- attend_right_noise_data[!attend_right_noise_data$S == '2',]
attend_right_noise_data <- attend_right_noise_data[!attend_right_noise_data$S == '1',]

attend_left_speech_data <- attend_left_speech_data[!attend_left_speech_data$S == '2',]
attend_left_speech_data <- attend_left_speech_data[!attend_left_speech_data$S == '1',]
attend_left_noise_data <- attend_left_noise_data[!attend_left_noise_data$S == '2',]
attend_left_noise_data <- attend_left_noise_data[!attend_left_noise_data$S == '1',]



attend_right_speech_data$Masker <- "speech"
attend_left_speech_data$Masker <- "speech"
attend_right_noise_data$Masker <- "noise"
attend_left_noise_data$Masker <- "noise"

attend_right_data <-rbind(attend_right_speech_data,attend_right_noise_data)
attend_left_data <-rbind(attend_left_speech_data,attend_left_noise_data)


# Add Ipsilateral/Contralateral Information
left_hemisphere_channels <- c(0,1,2,3,10,11,12,13)
right_hemisphere_channels <- c(4,5,6,7,8,9)

attend_right_data$Hemisphere <- NA
attend_left_data$Hemisphere <- NA

attend_right_data$Hemisphere[which(attend_right_data$Channel %in% right_hemisphere_channels)] <- "Ipsilateral"
attend_right_data$Hemisphere[which(attend_right_data$Channel %in% left_hemisphere_channels)] <- "Contralateral"

attend_left_data$Hemisphere[which(attend_left_data$Channel %in% right_hemisphere_channels)] <- "Contralateral"
attend_left_data$Hemisphere[which(attend_left_data$Channel %in% left_hemisphere_channels)] <- "Ipsilateral"

attend_right_data <- na.omit(attend_right_data)
attend_left_data <- na.omit(attend_left_data)

# Combine attend left and right
attend_right_data$Attend<-"right"
attend_left_data$Attend<-"left"
all_data<-rbind(attend_right_data,attend_left_data)

# Add ROI information
all_data$Roi<-NA
pfc_channels <- c(0,1,2,3,4,5)
stg_channels <- c(6,7,8,9,10,11,12,13)
all_data$Roi[which(all_data$Channel %in% pfc_channels)] <- "pfc"
all_data$Roi[which(all_data$Channel %in% stg_channels)] <- "stg"

# Organize Factors
to.factor <- c('S', 'Roi', 'Hemisphere', 'Masker', 'Attend', 'Spatialization')
all_data[, to.factor] <- lapply(all_data[, to.factor], as.factor)

all_data_cleaned <- na.omit(all_data)

# Check for normality, remove outliers
#shapiro.test(all_data_cleaned$MeanHb)


# Summary SE Function ####
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




# PFC, Speech Masker Model ####

all_data_cleaned_pfc_speechhbo <-  subset(all_data_cleaned, chromophore == "HbO" & Roi == "pfc" & Masker == "speech")
all_data_cleaned_pfc_speechhbr <-  subset(all_data_cleaned, chromophore == "HbR" & Roi == "pfc" & Masker == "speech")


z2_pfc_speech <- mixed(MeanHb ~ Spatialization*Hemisphere + (1|S) + (1|Channel),
                data= all_data_cleaned_pfc_speechhbo, 
                control = lmerControl(optimizer = "bobyqa"), method = 'LRT')

z2_pfc_speech

# Significant effect of spatialization, which we knew

# PFC, Noise Masker Model ####

all_data_cleaned_pfc_noisehbo <-  subset(all_data_cleaned,chromophore == "HbO" & Roi == "pfc" & Masker == "noise")
all_data_cleaned_pfc_noisehbo <-  subset(all_data_cleaned,chromophore == "HbR" & Roi == "pfc" & Masker == "noise")

z2_pfc_noise <- mixed(MeanHb ~ Spatialization*Hemisphere + (1|S) + (1|Channel),
                       data= all_data_cleaned_pfc_noisehbo, 
                       control = lmerControl(optimizer = "bobyqa"), method = 'LRT')

z2_pfc_noise

# No effects

# STG, Speech Masker Model ####

all_data_cleaned_stg_speechhbo <-  subset(all_data_cleaned, chromophore == "HbO" &  Roi == "stg" & Masker == "speech")
all_data_cleaned_stg_speechhbr <-  subset(all_data_cleaned, chromophore == "HbR" &  Roi == "stg" & Masker == "speech")

z2_stg_speech <- mixed(MeanHb ~ Spatialization*Hemisphere + (1|S) + (1|Channel),
                       data= all_data_cleaned_stg_speechhbo, 
                       control = lmerControl(optimizer = "bobyqa"), method = 'LRT')

z2_stg_speech

# Significant interaction between spatialization and hemisphere
# Pairwise Comparisons

EMM_stg_speech <- emmeans(z2_stg_speech, ~ Spatialization * Hemisphere)
pairs(EMM_stg_speech, simple = "Spatialization", adjust = "bonferroni")
pairs(EMM_stg_speech, simple = "Hemisphere", adjust = "bonferroni")


# STG, Noise Masker Model ####

all_data_cleaned_stg_noisehbo <-  subset(all_data_cleaned, chromophore == "HbO" & Roi == "stg" & Masker == "noise")
all_data_cleaned_stg_noisehbr <-  subset(all_data_cleaned, chromophore == "HbR" & Roi == "stg" & Masker == "noise")

z2_stg_noise <- mixed(MeanHb ~ Spatialization*Hemisphere + (1|S) + (1|Channel),
                       data= all_data_cleaned_stg_noisehbo, 
                       control = lmerControl(optimizer = "bobyqa"), method = 'LRT')

z2_stg_noise

# Significant interaction between spatialization and hemisphere
EMM_stg_noise <- emmeans(z2_stg_noise, ~ Spatialization * Hemisphere)
pairs(EMM_stg_noise, simple = "Spatialization", adjust = "bonferroni")
pairs(EMM_stg_noise, simple = "Hemisphere", adjust = "bonferroni")

# STG Plot ####
stg_se_data_speechhbo <- summarySE(all_data_cleaned_stg_speechhbo, measurevar="MeanHb", groupvars=c("S","Hemisphere","Spatialization"), na.rm = TRUE)
stg_se_data_speechhbo <- summarySE(stg_se_data_speechhbo, measurevar="MeanHb", groupvars=c("Hemisphere","Spatialization"), na.rm = TRUE)
plotspeechhbo <- ggplot(stg_se_data_speechhbo, aes(x=factor(Spatialization, level=c('ITD50', 'ITD500', 'ILD70n', 'ILD10')), y=MeanHb, shape = Hemisphere)) + 
  scale_shape_manual(values = c("Contralateral"=24, "Ipsilateral" = 22)) +
  geom_errorbar(aes(ymin=MeanHb-se, ymax=MeanHb+se), width=.1, position=position_dodge(width=0.5)) +
  geom_point(aes(),size = 4, position=position_dodge(width=0.5), color="red", fill="red") + 
  ggtitle("Speech Masker") +
  labs(x="",y="Mean \u0394HbO") +
  ylim(0,0.19) +
  scale_x_discrete(labels=c("ITD50" = "", "ITD500" = "","ILD70n" = "","ILD10" = "")) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_text(size=12)) +
  theme(legend.position="none") +
  geom_signif(xmin = c(0.9,1.9,2.9,3.9), xmax = c(1.1,2.1,3.1,4.1), y_position = c(0.11, 0.11, 0.11, 0.11), tip_length = 0.06, color="black", annotation = c("*","ns","ns","ns"), textsize = 5) +
  geom_signif(xmin = c(1,2,3), xmax = c(4,4,4), y_position = c(0.13, 0.15, 0.17), tip_length = 0.06, color="black", annotation = c("*","ns","p = 0.0624"), textsize = 5)

stg_se_data_noisehbo <- summarySE(all_data_cleaned_stg_noisehbo, measurevar="MeanHb", groupvars=c("S","Hemisphere","Spatialization"), na.rm = TRUE)
stg_se_data_noisehbo <- summarySE(stg_se_data_noisehbo, measurevar="MeanHb", groupvars=c("Hemisphere","Spatialization"), na.rm = TRUE)
plotnoisehbo <- ggplot(stg_se_data_noisehbo, aes(x=factor(Spatialization, level=c('ITD50', 'ITD500', 'ILD70n', 'ILD10')), y=MeanHb, shape  = Hemisphere)) + 
  scale_shape_manual(values = c("Contralateral"=24, "Ipsilateral" = 22)) +
  geom_errorbar(aes(ymin=MeanHb-se, ymax=MeanHb+se), width=.1, position=position_dodge(width=0.5)) +
  geom_point(aes(),size = 4, position=position_dodge(width=0.5), color="red", fill="red") + 
  ggtitle("Noise Masker") +
  labs(x="",y="") +
  ylim(0,0.19) +
  scale_x_discrete(labels=c("ITD50" = "", "ITD500" = "","ILD70n" = "","ILD10" = "")) + 
  theme_bw() + 
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_blank()) +
  geom_signif(xmin = c(0.9,1.9,2.9,3.9), xmax = c(1.1,2.1,3.1,4.1), y_position = c(0.11, 0.11, 0.11, 0.11), tip_length = 0.06, color="black", annotation = c("**","ns","ns","p = 0.077"), textsize = 5)


stg_se_data_speechhbr <- summarySE(all_data_cleaned_stg_speechhbr, measurevar="MeanHb", groupvars=c("S","Hemisphere","Spatialization"), na.rm = TRUE)
stg_se_data_speechhbr <- summarySE(stg_se_data_speechhbr, measurevar="MeanHb", groupvars=c("Hemisphere","Spatialization"), na.rm = TRUE)

plotspeechhbr <- ggplot(stg_se_data_speechhbr, aes(x=factor(Spatialization, level=c('ITD50', 'ITD500', 'ILD70n', 'ILD10')), y=MeanHb, shape = Hemisphere)) + 
  scale_shape_manual(values = c("Contralateral"=24, "Ipsilateral" = 22)) +
  geom_errorbar(aes(ymin=MeanHb-se, ymax=MeanHb+se), width=.1, position=position_dodge(width=0.5)) +
  geom_point(aes(),size = 4, position=position_dodge(width=0.5), color="blue", fill="blue") + 
  labs(x="",y="Mean \u0394HbR") +
  ylim(-0.03,0.005) +
  scale_x_discrete(labels=c("ITD50" = "Small\nITD", "ITD500" = "Large\nITD","ILD70n" = "Natural\nILD","ILD10" = "Broadband\nILD")) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_text(size=12)) +
  theme(legend.position="none")

stg_se_data_noisehbr <- summarySE(all_data_cleaned_stg_noisehbr, measurevar="MeanHb", groupvars=c("S","Hemisphere","Spatialization"), na.rm = TRUE)
stg_se_data_noisehbr <- summarySE(stg_se_data_noisehbr, measurevar="MeanHb", groupvars=c("Hemisphere","Spatialization"), na.rm = TRUE)


plotnoisehbr <- ggplot(stg_se_data_noisehbr, aes(x=factor(Spatialization, level=c('ITD50', 'ITD500', 'ILD70n', 'ILD10')), y=MeanHb, shape  = Hemisphere)) + 
  scale_shape_manual(values = c("Contralateral"=24, "Ipsilateral" = 22)) +
  geom_errorbar(aes(ymin=MeanHb-se, ymax=MeanHb+se), width=.1, position=position_dodge(width=0.5)) +
  geom_point(aes(),size = 4, position=position_dodge(width=0.5), color="blue", fill="blue") + 
  labs(x="",y="") +
  ylim(-0.03,0.005) +
  scale_x_discrete(labels=c("ITD50" = "Small\nITD", "ITD500" = "Large\nITD","ILD70n" = "Natural\nILD","ILD10" = "Broadband\nILD")) + 
  theme_bw() + 
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_blank())



grid.arrange(plotspeechhbo, plotnoisehbo, plotspeechhbr, plotnoisehbr, ncol=2,  widths = c(0.75,1), heights = c(4,2))

