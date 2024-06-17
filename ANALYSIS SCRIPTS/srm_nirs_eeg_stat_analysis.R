library(tidyr)
library(afex)

# Load in Data
speech_masker_data <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\all_subjects_mean_during_stim_speech_masker.csv")
noise_masker_data <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\all_subjects_mean_during_stim_noise_masker.csv")

colnames(speech_masker_data) <- c("S","Channel","ITD50","ITD500","ILD70n","ILD10")
colnames(noise_masker_data) <- c("S","Channel","ITD50","ITD500","ILD70n","ILD10")
speech_masker_data <- pivot_longer(speech_masker_data, cols=c("ITD50","ITD500","ILD70n","ILD10"), names_to = "Condition", values_to = "MeanHbO")
noise_masker_data <- pivot_longer(noise_masker_data, cols=c("ITD50","ITD500","ILD70n","ILD10"), names_to = "Condition", values_to = "MeanHbO")

# Organize Factors
to.factor <- c('S', 'Channel','Condition')
speech_masker_data[, to.factor] <- lapply(speech_masker_data[, to.factor], as.factor)
noise_masker_data[, to.factor] <- lapply(noise_masker_data[, to.factor], as.factor)

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

all_data$Roi <- as.factor(all_data$Roi)
all_data$Masker <- as.factor(all_data$Masker)

# R relevel factor will help reorder conditions if needed for contrast coding

# Simplest Model
z <- mixed(MeanHbO ~ Condition*Roi*Masker + (1|S),
      data= all_data, 
      control = lmerControl(optimizer = "bobyqa"), 
      method = 'LRT')

z

# Full Model
#z_full <- mixed(MeanHbO ~ Condition*Roi*Masker + (Condition*Roi*Masker|S),
#           data= all_data, 
#           control = lmerControl(optimizer = "bobyqa"), 
#           method = 'LRT')
#anova(z_full, z)

# What if we use Channel instead of ROI?
#z_channel <- mixed(MeanHbO ~ Condition*Channel*Masker + (1|S),
#           data= all_data, 
#           control = lmerControl(optimizer = "bobyqa"), 
#           method = 'LRT')