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
right_pfc_channels <- c(0,1,2,3)
left_pfc_channels <- c(4,5)
right_stg_channels <- c(6,7,8,9)
left_stg_channels <- c(10,11,12,13)
all_data_hbo$Roi[which(all_data_hbo$Channel %in% right_pfc_channels)] <- "right_pfc"
all_data_hbo$Roi[which(all_data_hbo$Channel %in% left_pfc_channels)] <- "left_pfc"
all_data_hbo$Roi[which(all_data_hbo$Channel %in% right_stg_channels)] <- "right_stg"
all_data_hbo$Roi[which(all_data_hbo$Channel %in% left_stg_channels)] <- "left_stg"

# Organize Factors
to.factor <- c('S', 'Roi','Spatialization')
all_data_hbo[, to.factor] <- lapply(all_data_hbo[, to.factor], as.factor)

all_data_hbo$Masker <- as.factor(all_data_hbo$Masker)
all_data_cleaned_hbo <- na.omit(all_data_hbo)

all_data_cleaned_hbo %>% group_by(Spatialization, Masker,Roi) %>% shapiro_test(MeanHb)
all_data_cleaned_pfc_speech_hbo <- subset(all_data_cleaned_hbo, Roi %in% c("right_pfc","left_pfc") & Masker == "speech")
all_data_cleaned_pfc_noise_hbo <- subset(all_data_cleaned_hbo, Roi %in% c("right_pfc","left_pfc") & Masker == "noise")
all_data_cleaned_stg_speech_hbo <- subset(all_data_cleaned_hbo, Roi %in% c("right_stg","left_stg") & Masker == "speech")
all_data_cleaned_stg_noise_hbo <- subset(all_data_cleaned_hbo, Roi %in% c("right_stg","left_stg") & Masker == "noise")

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
all_data_hbr$Roi[which(all_data_hbr$Channel %in% right_pfc_channels)] <- "right_pfc"
all_data_hbr$Roi[which(all_data_hbr$Channel %in% left_pfc_channels)] <- "left_pfc"
all_data_hbr$Roi[which(all_data_hbr$Channel %in% right_stg_channels)] <- "right_stg"
all_data_hbr$Roi[which(all_data_hbr$Channel %in% left_stg_channels)] <- "left_stg"
# Change ROI information

#all_data$Roi[all_data$Roi == "0"] <- "pfc"
#all_data$Roi[all_data$Roi == "1"] <- "stg"

# Organize Factors
to.factor <- c('S', 'Roi','Spatialization')
all_data_hbr[, to.factor] <- lapply(all_data_hbr[, to.factor], as.factor)

all_data_hbr$Masker <- as.factor(all_data_hbr$Masker)
all_data_cleaned_hbr <- na.omit(all_data_hbr)

all_data_cleaned_hbr %>% group_by(Spatialization, Masker,Roi) %>% shapiro_test(MeanHb)
all_data_cleaned_pfc_speech_hbr <- subset(all_data_cleaned_hbr, Roi %in% c("right_pfc","left_pfc") & Masker == "speech")
all_data_cleaned_pfc_noise_hbr <- subset(all_data_cleaned_hbr, Roi %in% c("right_pfc","left_pfc") & Masker == "noise")
all_data_cleaned_stg_speech_hbr <- subset(all_data_cleaned_hbr, Roi %in% c("right_stg","left_stg") & Masker == "speech")
all_data_cleaned_stg_noise_hbr <- subset(all_data_cleaned_hbr, Roi %in% c("right_stg","left_stg") & Masker == "noise")


##### Combining Data #####
all_data_cleaned_hbo$chromophore <- "HbO"
all_data_cleaned_hbr$chromophore <- "HbR"
all_data <- rbind(all_data_cleaned_hbo,all_data_cleaned_hbr)
all_data_pfc_speech <- subset(all_data, Roi %in% c("right_pfc","left_pfc") & Masker == "speech")
all_data_pfc_noise <- subset(all_data, Roi %in% c("right_pfc","left_pfc") & Masker == "noise")
all_data_stg_speech <- subset(all_data, Roi %in% c("right_stg","left_stg") & Masker == "speech")
all_data_stg_noise <- subset(all_data, Roi %in% c("right_stg","left_stg") & Masker == "noise")


# COMPARE WITH BEHAVIOR
hit_rates <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\RESULTS DATA\\SRM-NIRS-EEG-1_Hit_Rates.csv")
FA_rates <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\RESULTS DATA\\SRM-NIRS-EEG-1_FA_Rates.csv")
d_primes <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\RESULTS DATA\\SRM-NIRS-EEG-1_d_primes.csv")

# Remove unneeded columns, put in long format
hit_rates$OriginalVariableNames <- array(0:29)
colnames(hit_rates) <- c("S","ITD50_Noise","ITD500_Noise","ILD70n_Noise","ILD10_Noise","ITD50_Speech","ITD500_Speech","ILD70n_Speech","ILD10_Speech")
hit_rates <- pivot_longer(hit_rates, cols=c("ITD50_Noise","ITD500_Noise","ILD70n_Noise","ILD10_Noise","ITD50_Speech","ITD500_Speech","ILD70n_Speech","ILD10_Speech"),
                          names_to = c("Spatialization","Masker"), names_sep = "_", values_to = "HitRate")

FA_rates$OriginalVariableNames <- array(0:29)
colnames(FA_rates) <- c("S","ITD50","ITD500","ILD70n","ILD10")
FA_rates <- pivot_longer(FA_rates, cols=c("ITD50","ITD500","ILD70n","ILD10"),names_to = c("Spatialization"), values_to = "FARate")

d_primes$OriginalVariableNames <- array(0:29)
colnames(d_primes) <- c("S","ITD50","ITD500","ILD70n","ILD10")
d_primes <- pivot_longer(d_primes, cols=c("ITD50","ITD500","ILD70n","ILD10"),names_to = c("Spatialization"), values_to = "d_prime")

# Organize Factors
to.factor <- c('S','Masker','Spatialization')
hit_rates[, to.factor] <- lapply(hit_rates[, to.factor], as.factor)

to.factor <- c('S','Spatialization')
FA_rates[, to.factor] <- lapply(FA_rates[, to.factor], as.factor)
d_primes[, to.factor] <- lapply(d_primes[, to.factor], as.factor)

#### Overall PFC Mean Speech Masker Correlation Analysis ####
mean_hbo_to_compare <- subset(all_data_cleaned_hbo, Masker == "speech" & Roi %in% c("right_pfc","left_pfc"))
mean_hbo_to_compare <- mean_hbo_to_compare[!is.na(mean_hbo_to_compare$MeanHb), ]
mean_hbo_to_compare <- aggregate(MeanHb ~ S + Roi + Spatialization, 
                                 data = mean_hbo_to_compare, 
                                 FUN = mean)

df_list <- list(mean_hbo_to_compare, subset(hit_rates, Masker == "Speech"), FA_rates, d_primes)
all_to_compare_pfc_speech_masker <-  Reduce(function(x, y) merge(x, y, all=TRUE), df_list) 


# Function to compute regression stats per Spatialization level
get_lm_stats <- function(df, xvar, yvar = "MeanHb") {
  df %>%
    group_by(Spatialization) %>%
    do({
      model <- lm(as.formula(paste(yvar, "~", xvar)), data = .)
      tidy_model <- tidy(model)
      data.frame(
        Spatialization = unique(.$Spatialization),
        R2 = summary(model)$r.squared,
        pval = tidy_model$p.value[2]
      )
    })
}

# Function to generate a single ggplot
make_regression_plot <- function(df, xvar, x_label, title, xlim_vals, ylim_vals) {
  lm_stats <- get_lm_stats(df, xvar)
  lm_stats$label <- paste0("R² = ", round(lm_stats$R2, 3), 
                           ", p = ", round(lm_stats$pval, 3))
  
  # Predefine positions for annotation to avoid overlap
  lm_stats$label_x <- xlim_vals[1] + 0.05 * diff(xlim_vals)
  y_offset <- seq(ylim_vals[2] - 0.05, ylim_vals[2] - 0.25, length.out = nrow(lm_stats))
  lm_stats$label_y <- y_offset
  
  p <- ggplot(df, aes_string(x = xvar, y = "MeanHb", color = "Spatialization")) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, aes(group = Spatialization)) +
    geom_text(data = lm_stats, aes(x = label_x, y = label_y, label = label, color = Spatialization),
              inherit.aes = FALSE, hjust = 0, size = 4) +
    labs(title = title, x = x_label, y = "MeanHb") +
    theme_minimal() +
    xlim(xlim_vals) +
    ylim(ylim_vals)
  
  return(p)
}

# Generate and print the plots
p1 <- make_regression_plot(all_to_compare_pfc_speech_masker, "HitRate",
                           x_label = "HitRate",
                           title = "PFC HitRate vs MeanHbO (Speech Masker)",
                           xlim_vals = c(0, 1),
                           ylim_vals = c(-0.5, 0.5))
print(p1)

p2 <- make_regression_plot(all_to_compare_pfc_speech_masker, "FARate",
                           x_label = "FARate",
                           title = "PFC FARate vs MeanHbO (Speech Masker)",
                           xlim_vals = c(0, 1),
                           ylim_vals = c(-0.5, 0.5))
print(p2)

p3 <- make_regression_plot(all_to_compare_pfc_speech_masker, "d_prime",
                           x_label = "d_prime",
                           title = "PFC d_prime vs MeanHbO (Speech Masker)",
                           xlim_vals = c(-1, 7),
                           ylim_vals = c(-0.5, 0.5))
print(p3)



#### Overall STG Mean Speech Masker Correlation Analysis ####
mean_hbo_to_compare <- subset(all_data_cleaned_hbo, Masker == "speech" & Roi %in% c("left_stg","right_stg"))
mean_hbo_to_compare <- mean_hbo_to_compare[!is.na(mean_hbo_to_compare$MeanHb), ]
mean_hbo_to_compare <- aggregate(MeanHb ~ S + Roi + Spatialization, 
                                 data = mean_hbo_to_compare, 
                                 FUN = mean)

df_list <- list(mean_hbo_to_compare, subset(hit_rates, Masker == "Speech"), FA_rates, d_primes)
all_to_compare_stg_speech_masker <-  Reduce(function(x, y) merge(x, y, all=TRUE), df_list) 


# Generate and print the plots
p1 <- make_regression_plot(all_to_compare_stg_speech_masker, "HitRate",
                           x_label = "HitRate",
                           title = "STG HitRate vs MeanHbO (Speech Masker)",
                           xlim_vals = c(0, 1),
                           ylim_vals = c(-0.5, 0.5))
print(p1)

p2 <- make_regression_plot(all_to_compare_stg_speech_masker, "FARate",
                           x_label = "FARate",
                           title = "STG FARate vs MeanHbO (Speech Masker)",
                           xlim_vals = c(0, 1),
                           ylim_vals = c(-0.5, 0.5))
print(p2)

p3 <- make_regression_plot(all_to_compare_stg_speech_masker, "d_prime",
                           x_label = "d_prime",
                           title = "STG d_prime vs MeanHbO (Speech Masker)",
                           xlim_vals = c(-1, 7),
                           ylim_vals = c(-0.5, 0.5))
print(p3)







#### Overall PFC Mean Noise Masker Correlation Analysis ####
mean_hbo_to_compare <- subset(all_data_cleaned_hbo, Masker == "noise" & Roi %in% c("left_pfc","right_pfc"))
mean_hbo_to_compare <- mean_hbo_to_compare[!is.na(mean_hbo_to_compare$MeanHb), ]
mean_hbo_to_compare <- aggregate(MeanHb ~ S + Roi + Spatialization, 
                                 data = mean_hbo_to_compare, 
                                 FUN = mean)
df_list <- list(mean_hbo_to_compare, subset(hit_rates, Masker == "Noise"))
all_to_compare_pfc_noise_masker <-  Reduce(function(x, y) merge(x, y, all=TRUE), df_list) 


# Generate and print the plots
p1 <- make_regression_plot(all_to_compare_pfc_noise_masker, "HitRate",
                           x_label = "HitRate",
                           title = "PFC HitRate vs MeanHbO (Noise Masker)",
                           xlim_vals = c(0, 1),
                           ylim_vals = c(-0.5, 0.5))
print(p1)

#### Overall STG Mean Noise Masker Correlation Analysis ####
mean_hbo_to_compare <- subset(all_data_cleaned_hbo, Masker == "noise" & Roi %in% c("left_stg","right_stg"))
mean_hbo_to_compare <- mean_hbo_to_compare[!is.na(mean_hbo_to_compare$MeanHb), ]
mean_hbo_to_compare <- aggregate(MeanHb ~ S + Roi + Spatialization, 
                                 data = mean_hbo_to_compare, 
                                 FUN = mean)

df_list <- list(mean_hbo_to_compare, subset(hit_rates, Masker == "Noise"), FA_rates, d_primes)
all_to_compare_stg_noise_masker <-  Reduce(function(x, y) merge(x, y, all=TRUE), df_list) 


# Generate and print the plots
p1 <- make_regression_plot(all_to_compare_stg_noise_masker, "HitRate",
                           x_label = "HitRate",
                           title = "STG HitRate vs MeanHbO (Noise Masker)",
                           xlim_vals = c(0, 1),
                           ylim_vals = c(-0.5, 0.5))
print(p1)

