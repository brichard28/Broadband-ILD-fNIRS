# Author: Benjamin Richardson
# Uses information from srm_nirs_eeg_analyze_behavior.m


library(tidyverse)
library(ggpubr)
library(ggplot2)
library(rstatix)
library(afex)
library(dplyr)

hit_rates <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\RESULTS DATA\\SRM-NIRS-EEG-1_Hit_Rates.csv")

# Remove unneeded columns, put in long format
hit_rates$OriginalVariableNames <- array(1:26)
colnames(hit_rates) <- c("S","ITD50_Noise","ITD500_Noise","ILD70n_Noise","ILD10_Noise","ITD50_Speech","ITD500_Speech","ILD70n_Speech","ILD10_Speech")
hit_rates <- pivot_longer(hit_rates, cols=c("ITD50_Noise","ITD500_Noise","ILD70n_Noise","ILD10_Noise","ITD50_Speech","ITD500_Speech","ILD70n_Speech","ILD10_Speech"),
                          names_to = c("Condition","Masker"), names_sep = "_", values_to = "HitRate")

# Organize Factors
to.factor <- c('S','Masker','Condition')
hit_rates[, to.factor] <- lapply(hit_rates[, to.factor], as.factor)

# Summary Statistics
hit_rates %>% group_by(Condition, Masker) %>% get_summary_stats(HitRate, type = "mean_sd")

# Boxplot
bxp <- ggboxplot(hit_rates, x = "Condition", y = "HitRate", color = "Masker", palette = "jco")
bxp

# Check for normality, remove outliers
hit_rates %>% group_by(Condition, Masker) %>% shapiro_test(HitRate)

#hit_rates_no_outliers <- hit_rates %>%
#  group_by(Condition, Masker) %>%
#  identify_outliers("HitRate") %>%
#  filter(!is.outlier)



# Create a QQ plot
ggqqplot(hit_rates, "HitRate", ggtheme = theme_bw()) + facet_grid(Condition ~ Masker, labeller = "label_both")

# Run ANOVA
res.aov <- anova_test(data = hit_rates, dv = HitRate, wid = S, within = c(Condition, Masker))
get_anova_table(res.aov)

# Pairwise comparisons between maskers
pwc_masker <- hit_rates  %>% pairwise_t_test(HitRate ~ Masker, paired = TRUE, p.adjust.method = "bonferroni")
print(pwc_masker)

# Pairwise comparisons between maskers
pwc_condition <- hit_rates  %>% pairwise_t_test(HitRate ~ Condition, paired = TRUE, comparisons = list(c("ITD50","ITD500"),c("ILD70n","ILD10")), p.adjust.method = "bonferroni")
print(pwc_condition)

