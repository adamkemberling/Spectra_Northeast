####  Mean/Median Body Size


####  Packages  ####
library(tidyverse)
library(gmRi)
library(sizeSpectra)
library(tidyquant)


# Load processing functions
source(here::here("Code/R/Processing_Functions.R"))





####  Load Tidy Trawl Data  ####
finfish_trawl <- read_csv(here::here("Data/processed/finfish_trawl_data.csv"))
trawl_wigley <- read_csv(here::here("Data/processed/wigley_species_trawl_data.csv"))




####  Estimate Mean/Median Body Size  ####

# Lengths only for broader finfish community
finfish_lengths <- finfish_trawl %>% 
  group_size_metrics(
    size_data = ., 
    .abund_col = "numlen_adj", 
    .group_cols = c("survey_area", "est_year", "season"),
    has_weights = FALSE)


# Length and Weights for Wigley Species 
wigley_sizes <- trawl_wigley %>% 
  group_size_metrics(
    size_data = ., 
    .abund_col = "numlen_adj", 
    .group_cols = c("survey_area", "est_year", "season"),
    has_weights = TRUE, 
    .weight_col = "ind_weight_kg")




####  Plot Results  ####

# Finfish Community median length
finfish_trawl %>% distinct(comname) %>% nrow()
finfish_lengths %>% 
  mutate(yr_num = as.numeric(as.character(est_year))) %>% 
  ggplot(aes(yr_num, med_len_cm, color = season)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_ma(aes(linetype = "5-Year Moving Average"), n = 5, ma_fun = SMA) +
  geom_smooth(method = "lm", linewidth = 1, se = F, 
              aes(linetype = "Regression Fit")) +
  facet_wrap(~survey_area) +
  scale_color_gmri() +
  labs(
    title = "Median Length - Finfish Community",
     subtitle = "435 Finfish Species",
     y = "Length (cm)",
     x = "Year",
     color = "Season")





# Wigley Species Median Length
trawl_wigley %>% distinct(comname) %>% nrow()
wigley_sizes %>% 
  mutate(yr_num = as.numeric(as.character(est_year))) %>% 
  ggplot(aes(yr_num, med_len_cm, color = season)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_ma(aes(linetype = "5-Year Moving Average"), n = 5, ma_fun = SMA) +
  geom_smooth(method = "lm", linewidth = 1, se = F, 
              aes(linetype = "Regression Fit")) +
  facet_wrap(~survey_area) +
  scale_color_gmri() +
  labs(
    title = "Median Length - Wigley Community",
    subtitle = "74 Finfish Species",
    y = "Length (cm)",
    x = "Year",
    color = "Season")


# Wigley Species Median Length
wigley_sizes %>% 
  mutate(yr_num = as.numeric(as.character(est_year))) %>% 
  ggplot(aes(yr_num, med_wt_kg, color = season)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_ma(aes(linetype = "5-Year Moving Average"), n = 5, ma_fun = SMA) +
  geom_smooth(method = "lm", linewidth = 1, se = F, 
              aes(linetype = "Regression Fit")) +
  facet_wrap(~survey_area) +
  scale_color_gmri() +
  labs(
    title = "Median Weight - Wigley Community",
    subtitle = "74 Finfish Species",
    y = "Weight (kg)",
    x = "Year",
    color = "Season")



####  Save Length and Weight Metrics  ####
write_csv(finfish_lengths, here::here("Data/processed/finfish_species_length_summary.csv"))
write_csv(wigley_sizes, here::here("Data/processed/wigley_species_size_summary.csv"))

