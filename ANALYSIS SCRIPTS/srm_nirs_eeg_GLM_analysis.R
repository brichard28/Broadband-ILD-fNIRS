library(tidyr)
library(afex)

# Load in Data
uncorr_speech_masker_data <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\all_subjects_uncorr_GLM_speech_masker.csv")
colnames(uncorr_speech_masker_data) <- c("S","Channel","ITD50","ITD500","ILD70n","ILD10")
uncorr_speech_masker_data <- pivot_longer(uncorr_speech_masker_data, cols=c("ITD50","ITD500","ILD70n","ILD10"), names_to = "Condition", values_to = "Beta")

# Organize Factors
to.factor <- c('S', 'Channel','Condition')
uncorr_speech_masker_data[, to.factor] <- lapply(uncorr_speech_masker_data[, to.factor], as.factor)

# Test for main effect of condition
mixed(Beta ~ Condition*Channel + (1|S),
      data= uncorr_speech_masker_data, 
      control = lmerControl(optimizer = "bobyqa"), 
      method = 'LRT')