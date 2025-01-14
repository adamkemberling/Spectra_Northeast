# Process Spectra using Northeast Trawl Survey
# NUMLEN un-adjusted


####  Packages  ####
library(tidyverse)
library(gmRi)
library(sizeSpectra)
library(tidyquant)
library(conflicted)
library(patchwork)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)


# Load processing functions
source(here::here("Code/R/Processing_Functions.R"))



####  Load Tidy Trawl Data  ####

# All finfish species, length only spectra
finfish_trawl <- read_csv(here::here("Data/processed/finfish_trawl_data.csv"))

# Species with allometric l-w relationships
trawl_wigley <- read_csv(here::here("Data/processed/wigley_species_trawl_data.csv"))


# Table to join for swapping shorthand for long-hand names
area_df <- data.frame(
  area = c("Scotian Shelf", "Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight", "Northeast Shelf"),
  survey_area = c("SS", "GoM", "GB", "SNE", "MAB", "Northeast Shelf"),
  area_titles = c("Scotian Shelf", "Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight", "Northeast Shelf"))


#### Length Spectra  ####



##### All Species  ####
# All finfish length_spectra
length_binspectra <- group_binspecies_spectra(
  ss_input = finfish_trawl,
  grouping_vars = c("est_year", "season", "survey_area"),
  abundance_vals = "numlen",
  length_vals = "length_cm",
  use_weight = FALSE,
  isd_xmin = 1,
  global_min = TRUE,
  isd_xmax = NULL,
  global_max = FALSE,
  bin_width = 1, 
  vdiff = 4)



# Plot it
length_binspectra %>% 
  mutate(yr_num = as.numeric(as.character(est_year))) %>% 
  ggplot(aes(yr_num, b, color = season)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_ma(aes(linetype = "5-Year Moving Average"), n = 5, ma_fun = SMA) +
  geom_smooth(method = "lm", linewidth = 1, se = F, 
              aes(linetype = "Regression Fit")) +
  facet_wrap(~survey_area) +
  scale_color_gmri() +
  labs(title = "Length Spectra Slope - MLE Bins Method",
       subtitle = "Enforced xmin = 1, xmax = max(length_cm + 1)",
       y = "b",
       x = "Year",
       color = "Season")



# All finfish length_spectra
length_binspectra_shelf <- group_binspecies_spectra(
  ss_input = mutate(finfish_trawl, survey_area = "Northeast Shelf"),
  grouping_vars = c("est_year", "season", "survey_area"),
  abundance_vals = "numlen",
  length_vals = "length_cm",
  use_weight = FALSE,
  isd_xmin = 1,
  global_min = TRUE,
  isd_xmax = NULL,
  global_max = FALSE,
  bin_width = 1, 
  vdiff = 4)



# Plot it
length_binspectra_shelf %>% 
  mutate(yr_num = as.numeric(as.character(est_year))) %>% 
  ggplot(aes(yr_num, b, color = season)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_ma(aes(linetype = "5-Year Moving Average"),n = 5, ma_fun = SMA) +
  geom_smooth(method = "lm", linewidth = 1, se = F, 
              aes(linetype = "Regression Fit")) +
  facet_wrap(~survey_area) +
  scale_color_gmri() +
  labs(title = "Length Spectra Slope - MLE Bins Method",
       subtitle = "Enforced xmin = 1, xmax = max(length_cm + 1)",
       y = "b",
       x = "Year",
       color = "Season")





#####  Wigley Species   ####

# Wigley species length spectra
length_binspectra_wigley <- group_binspecies_spectra(
  ss_input = trawl_wigley,
  grouping_vars = c("est_year", "season", "survey_area"),
  abundance_vals = "numlen",
  length_vals = "length_cm",
  use_weight = FALSE,
  global_min = TRUE,
  isd_xmin = 1,
  isd_xmax = NULL,
  global_max = FALSE,
  bin_width = 4)




# Plot them
length_binspectra_wigley %>% 
  mutate(yr_num = as.numeric(as.character(est_year))) %>% 
  ggplot(aes(yr_num, b, color = season)) +
  geom_point(size = 1, alpha = 0.6) +
  #geom_ma(aes(linetype = "5-Year Moving Average"),n = 5, ma_fun = SMA) +
  geom_smooth(method = "lm", linewidth = 1, se = F, 
              aes(linetype = "Regression Fit")) +
  facet_wrap(~survey_area) +
  scale_color_gmri() +
  labs(title = "Length Spectra Slope - MLE Bins Method Wigley",
       subtitle = "Enforced xmin = 1, xmax = max(length_cm + 1)",
       y = "b",
       x = "Year",
       color = "Season")




####  Mass Spectra   ####


# For the Wigley species we can estimate individual weight from length
# then estimate size spectra using individual weights
bodymass_binspectra_wigley <- group_binspecies_spectra(
  ss_input = trawl_wigley,
  grouping_vars = c("est_year", "season", "survey_area"),
  abundance_vals = "numlen",
  length_vals = "length_cm",
  use_weight = TRUE,
  isd_xmin = 1,
  global_min = TRUE,
  isd_xmax = NULL,
  global_max = FALSE,
  bin_width = 1,
  vdiff = 2)



# And do a plot check
bodymass_binspectra_wigley %>% 
  filter(season %in% c("Spring", "Fall")) %>% 
  mutate(yr_num = as.numeric(as.character(est_year))) %>% 
  ggplot(aes(yr_num, b, color = season)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_smooth(method = "lm", linewidth = 1, se = F, 
              aes(linetype = "Regression Fit")) +
  facet_wrap(~survey_area) +
  scale_color_gmri() +
  labs(title = "Bodymass Spectra - MLE Bins Method Wigley",
       subtitle = "Enforced xmin = 1g, xmax = max(ind_weight_g)",
       y = "b",
       x = "Year",
       color = "Season")


# For the Wigley species we can estimate individual weight from length
# then estimate size spectra using individual weights
bodymass_binspectra_wigley_shelf <- group_binspecies_spectra(
  ss_input = mutate(trawl_wigley, survey_area = "Northeast Shelf"),
  grouping_vars = c("est_year", "season", "survey_area"),
  abundance_vals = "numlen",
  length_vals = "length_cm",
  use_weight = TRUE,
  isd_xmin = 1,
  global_min = TRUE,
  isd_xmax = NULL,
  global_max = FALSE,
  bin_width = 1,
  vdiff = 2)



# And do a plot check
bodymass_binspectra_wigley_shelf %>% 
  filter(season %in% c("Spring", "Fall")) %>% 
  mutate(yr_num = as.numeric(as.character(est_year))) %>% 
  ggplot(aes(yr_num, b, color = season)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_ma(aes(linetype = "5-Year Moving Average"),n = 5, ma_fun = SMA) +
  geom_smooth(method = "lm", linewidth = 1, se = F, 
              aes(linetype = "Regression Fit")) +
  facet_wrap(~survey_area) +
  scale_color_gmri() +
  labs(title = "Bodymass Spectra - MLE Bins Method Wigley",
       subtitle = "Enforced xmin = 1g, xmax = max(ind_weight_g)",
       y = "b",
       x = "Year",
       color = "Season")







####__________________________####
####  Fixed 16g WMIN Size Cutoff  ####


# Join the data back to the original so we can run the shelf at the same time
wigley_in_16g <- trawl_wigley %>%
  mutate(area = "Northeast Shelf", area_titles = "Northeast Shelf", survey_area = "Northeast Shelf") %>%
  bind_rows(left_join(trawl_wigley, area_df)) %>%
  filter(wmin_g >= 16)




# and again for 16g
mle_results_16g <- group_binspecies_spectra(
  ss_input = wigley_in_16g,
  grouping_vars = c("est_year", "season", "survey_area"),
  abundance_vals = "numlen",
  length_vals = "length_cm",
  use_weight = TRUE,
  isd_xmin = 16,
  global_min = TRUE,
  isd_xmax = NULL,
  global_max = FALSE,
  bin_width = 1,
  vdiff = 2) %>%
  mutate(est_year = as.numeric(est_year)) %>%
  left_join(area_df)


# Do we still have survey area etc.
mle_results_16g %>% distinct(survey_area, area, area_titles)


####  Minimum Cutoff Comparison:  #####



min_1_plot <- bind_rows(bodymass_binspectra_wigley, bodymass_binspectra_wigley_shelf) %>% glimpse()
  filter(season %in% c("Spring", "Fall")) %>% 
  mutate(yr_num = as.numeric(as.character(est_year)),
         survey_area = factor(survey_area, area_levels_all)) %>% 
  ggplot(aes(yr_num, b, color = season)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_ma(aes(linetype = "5-Year Moving Average"),n = 5, ma_fun = SMA) +
  geom_smooth(method = "lm", linewidth = 1, se = F, 
              aes(linetype = "Regression Fit")) +
  facet_wrap(~survey_area, ncol = 2) +
  scale_color_gmri() +
  labs(title = "Bodymass Spectra - MLE Bins Method Wigley",
       subtitle = "Enforced xmin = 1g, xmax = max(ind_weight_g)",
       y = "b",
       x = "Year",
       color = "Season")






# And do a plot check
min_16_plot <- mle_results_16g %>% 
  filter(season %in% c("Spring", "Fall")) %>% 
  mutate(
    yr_num = as.numeric(as.character(est_year)),
    survey_area = factor(
      survey_area, 
      levels = c("Northeast Shelf", "GoM", "GB", "SNE", "MAB"))) %>% 
  
  ggplot(aes(yr_num, b, color = season)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_ma(aes(linetype = "5-Year Moving Average"),n = 5, ma_fun = SMA) +
  geom_smooth(
    method = "lm", 
    linewidth = 1, 
    se = F, 
    aes(linetype = "Regression Fit")) +
  facet_wrap(~survey_area, ncol = 2, scales = "free") +
  scale_color_gmri() +
  labs(title = "Bodymass Spectra - MLE Bins Method Wigley",
       subtitle = "Fixed xmin = 16g, xmax = max(ind_weight_g)",
       y = "b",
       x = "Year",
       color = "Season")





library(patchwork)
min_1_plot | min_16_plot



####__________________________####


####  Saving Spectra Exponent Results  ####

# Length Based
write_csv(length_binspectra, here::here("Data/processed/finfish_length_spectra_numlen.csv"))
write_csv(length_binspectra_wigley, here::here("Data/processed/wigley_species_length_spectra_numlen.csv"))
write_csv(length_binspectra_shelf, here::here("Data/processed/shelfwide_finfish_length_spectra_numlen.csv"))


# Weight Based
write_csv(bodymass_binspectra_wigley, here::here("Data/processed/wigley_species_bodymass_spectra_numlen.csv"))
write_csv(bodymass_binspectra_wigley_shelf, here::here("Data/processed/shelfwide_wigley_species_bodymass_spectra_numlen.csv"))
# universal 16g wmin:
write_csv(mle_results_16g, here::here("Data/processed/wigley_species_min16_bodymass_spectra.csv"))


