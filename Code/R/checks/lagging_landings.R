#####  Checking Landings and Lags  ####

####  Packages  ####
library(emmeans)
library(tidyquant)
library(tidyverse)
library(gmRi)
library(ggeffects)
library(scales)
library(performance)
library(patchwork)
library(gtsummary)
library(gt)

# Theme
theme_set(theme_gmri(rect = element_rect(fill = "white", color = NA), legend.position = "bottom"))

#### Load existing models and Data  ####

# Median Length - Large Community
lc_len <- read_rds(here::here("Data/models_and_results/lc_medlen.RDS"))

# Length Spectra - Wigley Community
sc_lenb <- read_rds(here::here("Data/models_and_results/sc_lenspectra.RDS"))

# Median Weight - Small Community
sc_wt <- read_rds(here::here("Data/models_and_results/sc_medwt.RDS"))





#### Performing Lags  ####

data_l <- list(
  "Median Length" = lc_len$model_data
  #"Median Weight" = sc_wt$model_data
  #"Length Spectra" = sc_lenb$model_data
)



### Apply lags up to 4 years
data_l <- data_l %>% 
  map(function(x){
    x %>% 
      group_by(survey_area, season) %>% 
      arrange(yr_num) %>% 
      mutate(
        l1 = lag(landings, 1),
        l2 = lag(landings, 2),
        l3 = lag(landings, 3),
        l4 = lag(landings, 4),
        l5 = lag(landings, 5),
        l6 = lag(landings, 6),
        l7 = lag(landings, 7),
        l8 = lag(landings, 8),
        l9 = lag(landings, 9)
      ) %>% 
      ungroup()})




### Run Linear Models, Report significance and direction of landings
lag_models <- map(
  data_l,
  function(x){
    
    l1_mod <- lm(med_len_cm ~ survey_area + scale(l1) + scale(bot_temp), data = x)
    l2_mod <- lm(med_len_cm ~ survey_area + scale(l2) + scale(bot_temp), data = x)
    l3_mod <- lm(med_len_cm ~ survey_area + scale(l3) + scale(bot_temp), data = x)
    l4_mod <- lm(med_len_cm ~ survey_area + scale(l4) + scale(bot_temp), data = x)
    l5_mod <- lm(med_len_cm ~ survey_area + scale(l5) + scale(bot_temp), data = x)
    l6_mod <- lm(med_len_cm ~ survey_area + scale(l6) + scale(bot_temp), data = x)
    l7_mod <- lm(med_len_cm ~ survey_area + scale(l7) + scale(bot_temp), data = x)
    l8_mod <- lm(med_len_cm ~ survey_area + scale(l8) + scale(bot_temp), data = x)
    l9_mod <- lm(med_len_cm ~ survey_area + scale(l9) + scale(bot_temp), data = x)
    lags = list(
      "1" = l1_mod,
      "2" = l2_mod,
      "3" = l3_mod,
      "4" = l4_mod,
      "5" = l5_mod,
      "6" = l6_mod,
      "7" = l7_mod,
      "8" = l8_mod,
      "9" = l9_mod
    )
    
    
    
    
  }
)


# Do Median Length
med_len_lags <- map_dfr(
  lag_models$`Median Length`,
  function(x){
    coef_val <- coef(x)[5]
    data.frame(
      "coef_val" = coef_val,
      "model_rmse" = performance::rmse(x),
      "aic" = AIC(x)
    )
  }, .id = "n_years_lagged"
)


# Results plot
med_len_lags %>% 
  pivot_longer(cols = 2:4, names_to = "metric", values_to = "vals") %>% 
  ggplot(aes(n_years_lagged, vals)) +
  geom_line(aes(group = metric)) +
  geom_point() +
  facet_wrap(~metric, ncol = 1, scales = "free")
