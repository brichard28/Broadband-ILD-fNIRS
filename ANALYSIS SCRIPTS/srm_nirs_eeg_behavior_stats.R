# Author: Benjamin Richardson
# Uses information from srm_nirs_eeg_analyze_behavior.m


library(tidyverse)
library(ggpubr)
library(ggplot2)
library(rstatix)
library(afex)
library(dplyr)

####################################################
##    Hit Rates    ##
####################################################

hit_rates <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\RESULTS DATA\\SRM-NIRS-EEG-1_Hit_Rates.csv")

# Remove unneeded columns, put in long format
hit_rates$OriginalVariableNames <- array(1:31)
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

# Pairwise comparisons between maskers, given conditions
pwc_masker <- hit_rates  %>% pairwise_t_test(HitRate ~ Masker, paired = TRUE, p.adjust.method = "bonferroni")
print(pwc_masker)

# Pairwise comparisons between conditions, given masker
pwc_condition <- hit_rates  %>% pairwise_t_test(HitRate ~ Condition, paired = TRUE,  p.adjust.method = "bonferroni") # comparisons = list(c("ITD50","ITD500"),c("ILD70n","ILD10")),
print(pwc_condition)












####################################################
##    FA Rates    ##
####################################################

FA_rates <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\RESULTS DATA\\SRM-NIRS-EEG-1_FA_Rates.csv")

# Remove unneeded columns, put in long format
FA_rates$OriginalVariableNames <- array(1:31)
colnames(FA_rates) <- c("S","ITD50","ITD500","ILD70n","ILD10")
FA_rates <- pivot_longer(FA_rates, cols=c("ITD50","ITD500","ILD70n","ILD10"),
                          names_to = c("Condition"), values_to = "FARate")

# Organize Factors
to.factor <- c('S','Condition')
FA_rates[, to.factor] <- lapply(FA_rates[, to.factor], as.factor)

# Summary Statistics
FA_rates %>% group_by(Condition) %>% get_summary_stats(FARate, type = "mean_sd")

# Boxplot
bxp <- ggboxplot(FA_rates, x = "Condition", y = "FARate", palette = "jco")
bxp

# Check for normality, remove outliers
FA_rates %>% group_by(Condition) %>% shapiro_test(FARate)

# Create a QQ plot
ggqqplot(FA_rates, "FARate", ggtheme = theme_bw()) # + facet_grid(Condition) #, labeller = "label_both")

# Run ANOVA
res.aov <- anova_test(data = FA_rates, dv = FARate, wid = S, within = c(Condition))
get_anova_table(res.aov)

# Pairwise comparisons between conditions
pwc_masker <- FA_rates  %>% pairwise_t_test(FARate ~ Condition, paired = TRUE, p.adjust.method = "bonferroni")
print(pwc_masker)






####################################################
##    D primes    ##
####################################################

d_primes <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\RESULTS DATA\\SRM-NIRS-EEG-1_d_primes.csv")

# Remove unneeded columns, put in long format
d_primes$OriginalVariableNames <- array(1:31)
colnames(d_primes) <- c("S","ITD50","ITD500","ILD70n","ILD10")
d_primes <- pivot_longer(d_primes, cols=c("ITD50","ITD500","ILD70n","ILD10"),
                         names_to = c("Condition"), values_to = "DPrime")

# Organize Factors
to.factor <- c('S','Condition')
d_primes[, to.factor] <- lapply(d_primes[, to.factor], as.factor)

# Summary Statistics
d_primes %>% group_by(Condition) %>% get_summary_stats(DPrime, type = "mean_sd")

# Boxplot
bxp <- ggboxplot(d_primes, x = "Condition", y = "DPrime", palette = "jco")
bxp

# Check for normality, remove outliers
d_primes %>% group_by(Condition) %>% shapiro_test(DPrime)

# Create a QQ plot
ggqqplot(d_primes, "DPrime", ggtheme = theme_bw()) # + facet_grid(Condition) #, labeller = "label_both")

# Run ANOVA
res.aov <- anova_test(data = d_primes, dv = DPrime, wid = S, within = c(Condition))
get_anova_table(res.aov)

# Pairwise comparisons between conditions
pwc_masker <- d_primes  %>% pairwise_t_test(DPrime ~ Condition, paired = TRUE, p.adjust.method = "bonferroni")
print(pwc_masker)




















####################################################
##    Object Rates    ##
####################################################

object_rates <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\SRM-NIRS-EEG\\RESULTS DATA\\SRM-NIRS-EEG-1_object_Rates.csv")

# Remove unneeded columns, put in long format
object_rates$OriginalVariableNames <- array(1:31)
colnames(object_rates) <- c("S","ITD50_Noise","ITD500_Noise","ILD70n_Noise","ILD10_Noise","ITD50_Speech","ITD500_Speech","ILD70n_Speech","ILD10_Speech")
object_rates <- pivot_longer(object_rates, cols=c("ITD50_Noise","ITD500_Noise","ILD70n_Noise","ILD10_Noise","ITD50_Speech","ITD500_Speech","ILD70n_Speech","ILD10_Speech"),
                          names_to = c("Condition","Masker"), names_sep = "_", values_to = "ObjectRate")

# Organize Factors
to.factor <- c('S','Masker','Condition')
object_rates[, to.factor] <- lapply(object_rates[, to.factor], as.factor)

# Summary Statistics
object_rates %>% group_by(Condition, Masker) %>% get_summary_stats(ObjectRate, type = "mean_sd")

# Boxplot
bxp <- ggboxplot(object_rates, x = "Condition", y = "ObjectRate", color = "Masker", palette = "jco")
bxp

# Check for normality, remove outliers
object_rates %>% group_by(Condition, Masker) %>% shapiro_test(ObjectRate)


# Create a QQ plot
ggqqplot(object_rates, "ObjectRate", ggtheme = theme_bw()) + facet_grid(Condition ~ Masker, labeller = "label_both")

# Run ANOVA
res.aov <- anova_test(data = object_rates, dv = ObjectRate, wid = S, within = c(Condition, Masker))
get_anova_table(res.aov)

# Pairwise comparisons between maskers, given conditions
pwc_masker <- object_rates  %>% pairwise_t_test(ObjectRate ~ Masker, paired = TRUE, p.adjust.method = "bonferroni")
print(pwc_masker)

# Pairwise comparisons between conditions, given masker
pwc_condition <- object_rates  %>% pairwise_t_test(ObjectRate ~ Condition, paired = TRUE,  p.adjust.method = "bonferroni") # comparisons = list(c("ITD50","ITD500"),c("ILD70n","ILD10")),
print(pwc_condition)
