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

# Test for main effect of condition
z <- lmer(MeanHbO ~ Condition*Channel + (1|S),
      data= speech_masker_data, 
      control = lmerControl(optimizer = "bobyqa")) 
      #method = 'LRT')

summary(z)