#### All Results Figure & Table from Models
# attempt to stop manually copying things over and build the visuals together


####  Packages  ####
library(lme4)
library(lmerTest)
library(emmeans)
library(merTools)
library(tidyquant)
library(tidyverse)
library(gmRi)
library(ggeffects)
library(scales)
library(performance)

# Theme
theme_set(theme_gmri(rect = element_rect(fill = "white", color = NA)))


# ####  Source Everything to Ensure they are up-to-date  ####
# source(here::here("Code/R/modeling/lmer_models/large_community_medlen_lmer.R"))
# source(here::here("Code/R/modeling/lmer_models/large_community_lenspectra_lmer.R"))
# source(here::here("Code/R/modeling/lmer_models/wigley_community_medlen_lmer.R"))
# source(here::here("Code/R/modeling/lmer_models/wigley_community_medwt_lmer.R"))
# source(here::here("Code/R/modeling/lmer_models/wigley_community_lenspectra_lmer.R"))
# source(here::here("Code/R/modeling/lmer_models/wigley_community_bmspectra_lmer.R"))



#### Load Models/Data for Summary Plots  ####
sc_wt <- read_rds(here::here("Data/models_and_results/sc_medwt.RDS"))
sc_len <- read_rds(here::here("Data/models_and_results/sc_medlen.RDS"))
sc_lenb <- read_rds(here::here("Data/models_and_results/sc_lenspectra.RDS"))
sc_wtb <- read_rds(here::here("Data/models_and_results/sc_bmspectra.RDS"))
lc_len <- read_rds(here::here("Data/models_and_results/lc_medlen.RDS"))
lc_lenb <- read_rds(here::here("Data/models_and_results/lc_lenspectra.RDS"))




####  Median Length  ####
medlen_emmeans <- map_dfr(
  list(
    "Data-Limited Community" = lc_len$pass_1,
    "Well-Studied Community" = sc_len$pass_1
    ),
  function(x){
    region_posthoc <- emmeans(
    x, 
    list(pairwise ~ survey_area), 
    adjust = "tukey",
    type = "response")
  region_posthoc[[1]] %>% as.data.frame()
  }, .id = "community"
)


# Plot
ggplot(
  medlen_emmeans,
  aes(survey_area, emmean, ymin = lower.CL, ymax = upper.CL,  color = community)) +
  geom_pointrange(
    size = 1,
    position = position_dodge(width = 0.5), 
    alpha = 1) +
  scale_color_gmri() +
  labs(
    y = "Length (cm)",
    x = "Region",
    color = "Finfish Community",
    title = "EMMMean - Median Length")








#### Length Spectra ####

lenb_emmeans <- map_dfr(
  list(
    "Data-Limited Community" = lc_lenb$pass_1,
    "Well-Studied Community" = sc_lenb$pass_1
  ),
  function(x){
    region_posthoc <- emmeans(
      x, 
      list(pairwise ~ survey_area), 
      adjust = "tukey",
      type = "response")
    region_posthoc[[1]] %>% as.data.frame()
  }, .id = "community"
)



# Plot
ggplot(
  lenb_emmeans,
  aes(survey_area, emmean, ymin = lower.CL, ymax = upper.CL,  color = community)) +
  geom_pointrange(
    size = 1,
    position = position_dodge(width = 0.5), 
    alpha = 1) +
  scale_color_gmri() +
  labs(
    y = "Spectra Slope (length)",
    x = "Region",
    color = "Finfish Community",
    title = "EMMMean - Length Spectra")





####  Summary Table  ####

# Summary table - emmeans
pass1_mods <- list(
  "Well-Studied_Median Weight"    = sc_wt$pass_1,
  "Well-Studied_Median Length"    = sc_len$pass_1,
  "Well-Studied_Length Spectra"   = sc_lenb$pass_1,
  "Well-Studied_Bodymass Spectra" = sc_wtb$pass_1,
  "Data-Limited_Median Length"    = lc_len$pass_1,
  "Data-Limited_Length Spectra"   = lc_lenb$pass_1)


pass2_mods <- list(
  "Well-Studied_Median Weight"    = sc_wt$pass_2,
  "Well-Studied_Median Length"    = sc_len$pass_2,
  "Well-Studied_Length Spectra"   = sc_lenb$pass_2,
  "Well-Studied_Bodymass Spectra" = sc_wtb$pass_2,
  "Data-Limited_Median Length"    = lc_len$pass_2,
  "Data-Limited_Length Spectra"   = lc_lenb$pass_2)



# Summary Tables


# Get the regional means
