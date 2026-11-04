################################################################################
# Title:        Perrot et al. (2025) - multimodal analysis
# Subtitle:     effect_analysis.R
# Author:       Alain Luong (CentraleSupélec)
################################################################################

library(lmerTest)


dat_ms <- droplevels(subset(
  dat,
  group == "MS" & !is.na(Treatment)
))

dat_ms$Treatment <- relevel(
  factor(dat_ms$Treatment),
  ref = "None"
)

dat_ms$Sex <- as.character(dat_ms$Sex)
dat_ms$Sex[dat_ms$Sex %in% "H"] <- "M"
dat_ms$Sex <- factor(dat_ms$Sex, levels = c("F", "M"))
dat_ms$Subject_ID <- factor(dat_ms$Subject_ID)
dat_ms$Age_c <- dat_ms$Age - 40




print('### Linear mixed model to assess the effect of treatment on comp1 score ###')
modele_traitement <- lmer(
  comp1 ~ Treatment + (1 | Subject_ID),
  data = dat_ms
)
show(summary(modele_traitement))


print('### Linear mixed model to assess the effect of sex on comp1 score ###')
modele_sex <- lmer(
  comp1 ~ Sex + (1 | Subject_ID),
  data = dat_ms
)
show(summary(modele_sex))


print('### Linear mixed model to assess the effect of age on comp1 score ###')
modele_age <- lmer(
  comp1 ~ Age_c + (1 | Subject_ID),
  data = dat_ms
)
show(summary(modele_age))




