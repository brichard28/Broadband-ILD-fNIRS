
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(rstatix)
library(afex)
library(dplyr)
library(car)

# Data Preparation ####
# Load in Data
attend_right_speech_data <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\all_subjects_mean_during_stim_lateralization_target_right_speech_masker.csv")
attend_left_speech_data <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\all_subjects_mean_during_stim_lateralization_target_left_speech_masker.csv")
attend_right_noise_data <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\all_subjects_mean_during_stim_lateralization_target_right_noise_masker.csv")
attend_left_noise_data <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\all_subjects_mean_during_stim_lateralization_target_left_noise_masker.csv")

colnames(attend_right_speech_data) <- c("S","Channel","ITD50","ITD500","ILD70n","ILD10")
colnames(attend_right_noise_data) <- c("S","Channel","ITD50","ITD500","ILD70n","ILD10")
colnames(attend_left_speech_data) <- c("S","Channel","ITD50","ITD500","ILD70n","ILD10")
colnames(attend_left_noise_data) <- c("S","Channel","ITD50","ITD500","ILD70n","ILD10")

attend_right_speech_data <- pivot_longer(attend_right_speech_data, cols=c("ITD50","ITD500","ILD70n","ILD10"), names_to = "Spatialization", values_to = "MeanHbO")
attend_right_noise_data <- pivot_longer(attend_right_noise_data, cols=c("ITD50","ITD500","ILD70n","ILD10"), names_to = "Spatialization", values_to = "MeanHbO")
attend_left_speech_data <- pivot_longer(attend_left_speech_data, cols=c("ITD50","ITD500","ILD70n","ILD10"), names_to = "Spatialization", values_to = "MeanHbO")
attend_left_noise_data <- pivot_longer(attend_left_noise_data, cols=c("ITD50","ITD500","ILD70n","ILD10"), names_to = "Spatialization", values_to = "MeanHbO")

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
all_data_cleaned %>% group_by(Spatialization,Hemisphere,Roi) %>% shapiro_test(MeanHbO)


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

all_data_cleaned_pfc_speech <-  subset(all_data_cleaned, Roi == "pfc" & Masker == "speech")
z2_pfc_speech <- mixed(MeanHbO ~ Spatialization*Hemisphere + (1|S) + (1|Channel),
                data= all_data_cleaned_pfc_speech, 
                control = lmerControl(optimizer = "bobyqa"), method = 'LRT')

z2_pfc_speech

# Significant effect of spatialization, which we knew

# PFC, Noise Masker Model ####

all_data_cleaned_pfc_noise <-  subset(all_data_cleaned, Roi == "pfc" & Masker == "noise")
z2_pfc_noise <- mixed(MeanHbO ~ Spatialization*Hemisphere + (1|S) + (1|Channel),
                       data= all_data_cleaned_pfc_noise, 
                       control = lmerControl(optimizer = "bobyqa"), method = 'LRT')

z2_pfc_noise

# No effects

# STG, Speech Masker Model ####

all_data_cleaned_stg_speech <-  subset(all_data_cleaned, Roi == "stg" & Masker == "speech")
z2_stg_speech <- mixed(MeanHbO ~ Spatialization*Hemisphere + (1|S) + (1|Channel),
                       data= all_data_cleaned_stg_speech, 
                       control = lmerControl(optimizer = "bobyqa"), method = 'LRT')

z2_stg_speech

# Significant interaction between spatialization and hemisphere

# Pairwise Comparisons (treatment coding)
# Set contra always as reference
all_data_cleaned_stg_speech$Hemisphere <- relevel(all_data_cleaned_stg_speech$Hemisphere, "Contralateral")

# Ipsi vs. Contra within ITD50
# LANGUAGE FOR PAPER: for analyses where we examined effects of lateralization, we only included random by subject intercepts because models with random by channel intercepts failed to converge
stg_speech_lmer_itd50 <- lmer(MeanHbO ~ Hemisphere + (1|S),
                              data= subset(all_data_cleaned_stg_speech, Spatialization == "ITD50"),
                              control = lmerControl(optimizer = "bobyqa"))
#summary(stg_speech_lmer_itd50)

# Ipsi vs. Contra within ITD500
stg_speech_lmer_itd500 <- lmer(MeanHbO ~ Hemisphere + (1|S),
                              data= subset(all_data_cleaned_stg_speech, Spatialization == "ITD500"),
                              control = lmerControl(optimizer = "bobyqa"))
#summary(stg_speech_lmer_itd500)

# Ipsi vs. Contra within ILD70n
stg_speech_lmer_ild70n <- lmer(MeanHbO ~ Hemisphere + (1|S),
                               data= subset(all_data_cleaned_stg_speech, Spatialization == "ILD70n"),
                               control = lmerControl(optimizer = "bobyqa"))
#summary(stg_speech_lmer_ild70n)

# Ipsi vs. Contra within ITD500
stg_speech_lmer_ild10 <- lmer(MeanHbO ~ Hemisphere + (1|S),
                               data= subset(all_data_cleaned_stg_speech, Spatialization == "ILD10"),
                               control = lmerControl(optimizer = "bobyqa"))
#summary(stg_speech_lmer_ild10)



# STG, Noise Masker Model ####

all_data_cleaned_stg_noise <-  subset(all_data_cleaned, Roi == "stg" & Masker == "noise")
z2_stg_noise <- mixed(MeanHbO ~ Spatialization*Hemisphere + (1|S) + (1|Channel),
                       data= all_data_cleaned_stg_noise, 
                       control = lmerControl(optimizer = "bobyqa"), method = 'LRT')

z2_stg_noise

# Significant interaction between spatialization and hemisphere

# Pairwise Comparisons (treatment coding)
# Set contra always as reference
all_data_cleaned_stg_noise$Hemisphere <- relevel(all_data_cleaned_stg_noise$Hemisphere, "Contralateral")

# Ipsi vs. Contra within ITD50
stg_noise_lmer_itd50 <- lmer(MeanHbO ~ Hemisphere + (1|S),
                              data= subset(all_data_cleaned_stg_noise, Spatialization == "ITD50"),
                              control = lmerControl(optimizer = "bobyqa"))
#summary(stg_noise_lmer_itd50)

# Ipsi vs. Contra within ITD500
stg_noise_lmer_itd500 <- lmer(MeanHbO ~ Hemisphere + (1|S),
                               data= subset(all_data_cleaned_stg_noise, Spatialization == "ITD500"),
                               control = lmerControl(optimizer = "bobyqa"))
summary(stg_noise_lmer_itd500)

# Ipsi vs. Contra within ILD70n
stg_noise_lmer_ild70n <- lmer(MeanHbO ~ Hemisphere + (1|S),
                               data= subset(all_data_cleaned_stg_noise, Spatialization == "ILD70n"),
                               control = lmerControl(optimizer = "bobyqa"))
#summary(stg_noise_lmer_ild70n)

# Ipsi vs. Contra within ITD500
stg_noise_lmer_ild10 <- lmer(MeanHbO ~ Hemisphere + (1|S),
                              data= subset(all_data_cleaned_stg_noise, Spatialization == "ILD10"),
                              control = lmerControl(optimizer = "bobyqa"))
# summary(stg_noise_lmer_ild10)


# PFC Plot ####
pfc_se_data_speech <- summarySE(all_data_cleaned_pfc_speech, measurevar="MeanHbO", groupvars=c("S","Hemisphere","Spatialization"), na.rm = TRUE)
pfc_se_data_speech <- summarySE(pfc_se_data_speech, measurevar="MeanHbO", groupvars=c("Hemisphere","Spatialization"), na.rm = TRUE)
plotspeech <- ggplot(pfc_se_data_speech, aes(x=factor(Spatialization, level=c('ITD50', 'ITD500', 'ILD70n', 'ILD10')), y=MeanHbO, shape = Hemisphere, color = Spatialization)) + 
  scale_color_manual(values = c("ITD50" = "red","ITD500" =  "blue","ILD70n" = "magenta","ILD10" = "green")) +
  geom_errorbar(aes(ymin=MeanHbO-se, ymax=MeanHbO+se), width=.1, position=position_dodge(width=0.5)) +
  geom_point(aes(color = Spatialization),size = 4, position=position_dodge(width=0.5)) + 
  ggtitle("Speech Masker PFC") +
  labs(x="",y="Mean \u0394HbO") +
  ylim(-0.005,0.12) +
  scale_x_discrete(labels=c("ITD50" = "Small\nITD", "ITD500" = "Large\nITD","ILD70n" = "Natural\nILD","ILD10" = "Broadband\nILD")) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_text(size=12)) +
  theme(legend.position="none")

pfc_se_data_noise <- summarySE(all_data_cleaned_pfc_noise, measurevar="MeanHbO", groupvars=c("S","Hemisphere","Spatialization"), na.rm = TRUE)
pfc_se_data_noise <- summarySE(pfc_se_data_noise, measurevar="MeanHbO", groupvars=c("Hemisphere","Spatialization"), na.rm = TRUE)
plotnoise <- ggplot(pfc_se_data_noise, aes(x=factor(Spatialization, level=c('ITD50', 'ITD500', 'ILD70n', 'ILD10')), y=MeanHbO, shape = Hemisphere, color = Spatialization)) + 
  scale_color_manual(values = c("ITD50" = "red","ITD500" =  "blue","ILD70n" = "magenta","ILD10" = "green"), guide = 'none') +
  geom_errorbar(aes(ymin=MeanHbO-se, ymax=MeanHbO+se), width=.1, position=position_dodge(width=0.5)) +
  geom_point(aes(color = Spatialization),size = 4, position=position_dodge(width=0.5)) + 
  ggtitle("Noise Masker PFC") +
  labs(x="",y="") +
  ylim(-0.005,0.12) +
  scale_x_discrete(labels=c("ITD50" = "Small\nITD", "ITD500" = "Large\nITD","ILD70n" = "Natural\nILD","ILD10" = "Broadband\nILD")) + 
  theme_bw() + 
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_blank())

grid.arrange(plotspeech, plotnoise, ncol=2,  widths = c(0.75, 1))


# STG Plot ####
stg_se_data_speech <- summarySE(all_data_cleaned_stg_speech, measurevar="MeanHbO", groupvars=c("S","Hemisphere","Spatialization"), na.rm = TRUE)
stg_se_data_speech <- summarySE(stg_se_data_speech, measurevar="MeanHbO", groupvars=c("Hemisphere","Spatialization"), na.rm = TRUE)
plotspeech <- ggplot(stg_se_data_speech, aes(x=factor(Spatialization, level=c('ITD50', 'ITD500', 'ILD70n', 'ILD10')), y=MeanHbO, color = Hemisphere, fill = Hemisphere, shape = Spatialization)) + 
  scale_shape_manual(values = c("ITD50"=21, "ITD500" = 24, "ILD70n" = 22, "ILD10" = 25)) +
  scale_color_manual(values = c("Contralateral"="green","Ipsilateral" = "magenta")) +
  scale_fill_manual(values = c("Contralateral"="green","Ipsilateral" = "magenta")) +
  geom_errorbar(aes(ymin=MeanHbO-se, ymax=MeanHbO+se), width=.1, position=position_dodge(width=0.5)) +
  geom_point(aes(color = Hemisphere, shape = Spatialization),size = 4, position=position_dodge(width=0.5)) + 
  ggtitle("Speech Masker STG") +
  labs(x="",y="Mean \u0394HbO") +
  ylim(-0.005,0.12) +
  scale_x_discrete(labels=c("ITD50" = "Small\nITD", "ITD500" = "Large\nITD","ILD70n" = "Natural\nILD","ILD10" = "Broadband\nILD")) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_text(size=12)) +
  theme(legend.position="none") +
  geom_signif(xmin = c(0.9,1.9,2.9,3.9), xmax = c(1.1,2.1,3.1,4.1), y_position = c(0.11, 0.11, 0.11, 0.11), tip_length = 0.06, color="black", annotation = c("*","ns","ns","ns"), textsize = 5)
  

stg_se_data_noise <- summarySE(all_data_cleaned_stg_noise, measurevar="MeanHbO", groupvars=c("S","Hemisphere","Spatialization"), na.rm = TRUE)
stg_se_data_noise <- summarySE(stg_se_data_noise, measurevar="MeanHbO", groupvars=c("Hemisphere","Spatialization"), na.rm = TRUE)
plotnoise <- ggplot(stg_se_data_noise, aes(x=factor(Spatialization, level=c('ITD50', 'ITD500', 'ILD70n', 'ILD10')), y=MeanHbO, color = Hemisphere, fill = Hemisphere, shape  = Spatialization)) + 
  scale_shape_manual(values = c("ITD50"=21, "ITD500" = 24, "ILD70n" = 22, "ILD10" = 25)) +
  scale_color_manual(values = c("Contralateral"="green","Ipsilateral" = "magenta")) +
  scale_fill_manual(values = c("Contralateral"="green","Ipsilateral" = "magenta")) +
  geom_errorbar(aes(ymin=MeanHbO-se, ymax=MeanHbO+se), width=.1, position=position_dodge(width=0.5)) +
  geom_point(aes(color = Hemisphere, shape = Spatialization),size = 4, position=position_dodge(width=0.5)) + 
  ggtitle("Noise Masker STG") +
  labs(x="",y="") +
  ylim(-0.005,0.12) +
  scale_x_discrete(labels=c("ITD50" = "Small\nITD", "ITD500" = "Large\nITD","ILD70n" = "Natural\nILD","ILD10" = "Broadband\nILD")) + 
  theme_bw() + 
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_blank()) +
  geom_signif(xmin = c(0.9,1.9,2.9,3.9), xmax = c(1.1,2.1,3.1,4.1), y_position = c(0.11, 0.11, 0.11, 0.11), tip_length = 0.06, color="black", annotation = c("**","ns","ns","ns"), textsize = 5)


grid.arrange(plotspeech, plotnoise, ncol=2,  widths = c(0.75, 1))

