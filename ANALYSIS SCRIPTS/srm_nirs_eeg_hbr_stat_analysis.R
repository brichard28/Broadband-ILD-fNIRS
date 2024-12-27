library(tidyverse)
library(ggpubr)
library(ggplot2)
library(rstatix)
library(afex)
library(dplyr)
library(car)
require(gridExtra)

# Load in Data
speech_masker_data_hbr <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\all_subjects_mean_during_stim_speech_masker_hbr.csv")
noise_masker_data_hbr <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\ANALYSIS SCRIPTS\\Eli Analysis\\all_subjects_mean_during_stim_noise_masker_hbr.csv")

colnames(speech_masker_data_hbr) <- c("S","Channel","ITD50","ITD500","ILD70n","ILD10")
colnames(noise_masker_data_hbr) <- c("S","Channel","ITD50","ITD500","ILD70n","ILD10")
speech_masker_data_hbr <- pivot_longer(speech_masker_data_hbr, cols=c("ITD50","ITD500","ILD70n","ILD10"), names_to = "Spatialization", values_to = "MeanHbR")
noise_masker_data_hbr <- pivot_longer(noise_masker_data_hbr, cols=c("ITD50","ITD500","ILD70n","ILD10"), names_to = "Spatialization", values_to = "MeanHbR")


# Combine speech and noise
speech_masker_data_hbr$Masker<-"speech"
noise_masker_data_hbr$Masker<-"noise"
all_data_hbr<-rbind(speech_masker_data_hbr,noise_masker_data_hbr)

# Add ROI information
all_data$Roi<-NA
pfc_channels <- c(0,1,2,3,4,5)
stg_channels <- c(6,7,8,9,10,11,12,13)
all_data_hbr$Roi[which(all_data$Channel %in% pfc_channels)] <- "pfc"
all_data_hbr$Roi[which(all_data$Channel %in% stg_channels)] <- "stg"

# Change ROI information

#all_data$Roi[all_data$Roi == "0"] <- "pfc"
#all_data$Roi[all_data$Roi == "1"] <- "stg"

# Organize Factors
to.factor <- c('S', 'Roi','Spatialization')
all_data_hbr[, to.factor] <- lapply(all_data_hbr[, to.factor], as.factor)

all_data_hbr$Masker <- as.factor(all_data_hbr$Masker)
all_data_cleaned_hbr <- na.omit(all_data_hbr)

# Take mean across ROI
#r1<-with(all_data_cleaned, tapply(MeanHbR, Roi, mean))


# Boxplot PFC
#bxp_pfc <- ggboxplot(subset(all_data_cleaned,Roi=="pfc"), x = "Masker", y = "MeanHbR", color = "Spatialization", palette = "jco")
#bxp_pfc

# Boxplot STG
#bxp_stg <- ggboxplot(subset(all_data_cleaned,Roi=="stg"), x = "Masker", y = "MeanHbR", color = "Spatialization", palette = "jco")
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


# Check for normality, remove outliers
all_data_cleaned_hbr %>% group_by(Spatialization, Masker,Roi) %>% shapiro_test(MeanHbR)
all_data_cleaned_pfc_speech_hbr <- subset(all_data_cleaned_hbr, Roi == "pfc" & Masker == "speech")
all_data_cleaned_pfc_noise_hbr <- subset(all_data_cleaned_hbr, Roi == "pfc" & Masker == "noise")
all_data_cleaned_stg_speech_hbr <- subset(all_data_cleaned_hbr, Roi == "stg" & Masker == "speech")
all_data_cleaned_stg_noise_hbr <- subset(all_data_cleaned_hbr, Roi == "stg" & Masker == "noise")



#########################
# Mixed Effects models #
#########################

# PFC Speech
model_pfc_speech_hbr <- mixed(MeanHbR ~ Spatialization + (1|S) + (1|Channel),
                       data= all_data_cleaned_pfc_speech_hbr, 
                       control = lmerControl(optimizer = "bobyqa"), method = 'LRT')
model_pfc_speech_hbr


# PFC Noise
model_pfc_noise_hbr <- mixed(MeanHbR ~ Spatialization + (1|S) + (1|Channel),
                          data= all_data_cleaned_pfc_noise_hbr, 
                          control = lmerControl(optimizer = "bobyqa"), method = 'LRT')
model_pfc_noise_hbr

# STG Speech
model_stg_speech_hbr <- mixed(MeanHbR ~ Spatialization + (1|S) + (1|Channel),
                          data= all_data_cleaned_stg_speech_hbr, 
                          control = lmerControl(optimizer = "bobyqa"), method = 'LRT')
model_stg_speech_hbr

# STG Noise
model_stg_noise_hbr <- mixed(MeanHbR ~ Spatialization + (1|S) + (1|Channel),
                          data= all_data_cleaned_stg_noise_hbr, 
                          control = lmerControl(optimizer = "bobyqa"), method = 'LRT')
model_stg_noise_hbr


#########################
#        Plots          #
#########################


##### PFC Plot ##########
pfc_se_data_speech <- summarySE(all_data_cleaned_pfc_speech, measurevar="MeanHbR", groupvars=c("S","Masker","Spatialization"), na.rm = TRUE)
pfc_se_data_speech <- summarySE(pfc_se_data_speech, measurevar="MeanHbR", groupvars=c("Masker","Spatialization"), na.rm = TRUE)
plotspeech <- ggplot(pfc_se_data_speech, aes(x=factor(Spatialization, level=c('ITD50', 'ITD500', 'ILD70n', 'ILD10')), y=MeanHbR, color = Spatialization)) + 
  scale_color_manual(values = c("ITD50" = "red","ITD500" =  "blue","ILD70n" = "magenta","ILD10" = "green")) +
  geom_errorbar(aes(ymin=MeanHbR-se, ymax=MeanHbR+se, color=Spatialization), width=.1, position=position_dodge(width=0.5)) +
  geom_point(aes(color = Spatialization),size = 4, position=position_dodge(width=0.5)) + 
  ggtitle("Speech Masker PFC") +
  labs(x="",y="Mean \u0394HbR (\u03BCM)") +
  ylim(-0.03,0.03) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_text(size=12)) +
  scale_x_discrete(labels=c("ITD50" = "Small\nITD", "ITD500" = "Large\nITD","ILD70n" = "Natural\nILD","ILD10" = "Broadband\nILD")) +
  theme(legend.position="none")
#geom_signif(comparisons = list(c("ITD50","ITD500")), y_position = 0.09, tip_length = 0, color="black", annotation = c("***"), textsize = 5) +
#geom_signif(comparisons = list(c("ITD50","ILD70n")), y_position = 0.10, tip_length = 0, color="black", annotation = c("***"), textsize = 5) +
#geom_signif(comparisons = list(c("ITD50","ILD10")), y_position = 0.11, tip_length = 0, color="black", annotation = c("***"), textsize = 5)


pfc_se_data_noise <- summarySE(all_data_cleaned_pfc_noise, measurevar="MeanHbR", groupvars=c("S","Masker","Spatialization"), na.rm = TRUE)
pfc_se_data_noise <- summarySE(pfc_se_data_noise, measurevar="MeanHbR", groupvars=c("Masker","Spatialization"), na.rm = TRUE)
plotnoise <- ggplot(pfc_se_data_noise, aes(x=factor(Spatialization, level=c('ITD50', 'ITD500', 'ILD70n', 'ILD10')), y=MeanHbR, color = Spatialization)) + 
  scale_color_manual(values = c("ITD50" = "red","ITD500" =  "blue","ILD70n" = "magenta","ILD10" = "green")) +
  geom_errorbar(aes(ymin=MeanHbR-se, ymax=MeanHbR+se, color=Spatialization), width=.1, position=position_dodge(width=0.5)) +
  geom_point(aes(color = Spatialization),size = 4, position=position_dodge(width=0.5)) + 
  ggtitle("Noise Masker PFC") +
  labs(x="",y="") +
  ylim(-0.03,0.03) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_blank()) +
  scale_x_discrete(labels=c("ITD50" = "Small\nITD", "ITD500" = "Large\nITD","ILD70n" = "Natural\nILD","ILD10" = "Broadband\nILD")) +
  theme(legend.position="none")

grid.arrange(plotspeech,plotnoise, ncol=2, widths = c(1,0.9))




######### STG Plot ###############

stg_se_data_speech <- summarySE(all_data_cleaned_stg_speech, measurevar="MeanHbR", groupvars=c("S","Spatialization"), na.rm=TRUE)
stg_se_data_speech <- summarySE(stg_se_data_speech, measurevar="MeanHbR", groupvars=c("Spatialization"),  na.rm=TRUE)
plotspeech <- ggplot(stg_se_data_speech, aes(x=factor(Spatialization, level=c('ITD50', 'ITD500', 'ILD70n', 'ILD10')), y=MeanHbR, color = Spatialization)) + 
  scale_color_manual(values = c("ITD50" = "red","ITD500" =  "blue","ILD70n" = "magenta","ILD10" = "green")) +
  geom_errorbar(aes(ymin=MeanHbR-se, ymax=MeanHbR+se, color=Spatialization), width=.1, position=position_dodge(width=0.5)) +
  geom_point(aes(color = Spatialization),size = 4, position=position_dodge(width=0.5)) +  
  ggtitle("Speech Masker STG") +
  labs(x="",y="Mean \u0394HbR (\u03BCM)", parse=TRUE) +
  ylim(-0.03,0.03) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_text(size=12)) +
  scale_x_discrete(labels=c("ITD50" = "Small\nITD", "ITD500" = "Large\nITD","ILD70n" = "Natural\nILD","ILD10" = "Broadband\nILD")) +
  theme(legend.position="none")
#geom_signif(comparisons = list(c("ITD50","ILD10")), y_position = 0.090, tip_length = 0, color="black", annotation = c("**"), textsize = 5) +
#geom_signif(comparisons = list(c("ITD500","ILD10")), y_position = 0.080, tip_length = 0, color="black", annotation =c("*"), textsize = 5) +
#geom_signif(comparisons = list(c("ILD70n","ILD10")), y_position = 0.070, tip_length = 0, color="black", annotation =c("**"), textsize = 5) 


stg_se_data_noise <- summarySE(all_data_cleaned_stg_noise, measurevar="MeanHbR", groupvars=c("S","Spatialization"), na.rm=TRUE)
stg_se_data_noise <- summarySE(stg_se_data_noise, measurevar="MeanHbR", groupvars=c("Spatialization"),  na.rm=TRUE)
plotnoise <- ggplot(stg_se_data_noise, aes(x=factor(Spatialization, level=c('ITD50', 'ITD500', 'ILD70n', 'ILD10')), y=MeanHbR, color = Spatialization)) + 
  scale_color_manual(values = c("ITD50" = "red","ITD500" =  "blue","ILD70n" = "magenta","ILD10" = "green")) +
  geom_errorbar(aes(ymin=MeanHbR-se, ymax=MeanHbR+se, color=Spatialization), width=.1, position=position_dodge(width=0.5)) +
  geom_point(aes(color = Spatialization),size = 4, position=position_dodge(width=0.5)) +  
  ggtitle("Noise Masker STG") +
  labs(x="", y="") +
  ylim(-0.03,0.03) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_blank()) +
  scale_x_discrete(labels=c("ITD50" = "Small\nITD", "ITD500" = "Large\nITD","ILD70n" = "Natural\nILD","ILD10" = "Broadband\nILD")) +
  theme(legend.position="none") 

grid.arrange(plotspeech,plotnoise, ncol=2, widths = c(1,0.9))

