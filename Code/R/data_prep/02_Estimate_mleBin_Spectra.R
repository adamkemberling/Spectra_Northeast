# Process Spectra using Northeast Trawl Survey


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
finfish_trawl <- read_csv(here::here("Data/processed/finfish_trawl_data.csv"))
trawl_wigley <- read_csv(here::here("Data/processed/wigley_species_trawl_data.csv"))



#### Length Spectra All Species - By Region  ####


# Test 1 subgroup
finfish_trawl %>% 
  filter(est_year == 1970,
         survey_area == "GB",
         season == "Fall") %>% 
  # run the function
  group_binspecies_spectra(
    ss_input = .,
    grouping_vars = c("est_year", "season", "survey_area"),
    abundance_vals = "numlen_adj",
    length_vals = "length_cm",
    use_weight = FALSE,
    global_min = TRUE,
    isd_xmin = 1,
    isd_xmax = NULL,
    global_max = FALSE,
    bin_width = 1, 
    vdiff = 4)




# All finfish length_spectra
length_binspectra <- group_binspecies_spectra(
  ss_input = finfish_trawl,
  grouping_vars = c("est_year", "season", "survey_area"),
  abundance_vals = "numlen_adj",
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



#### Full Shelf ####

# All finfish length_spectra
length_binspectra_shelf <- group_binspecies_spectra(
  ss_input = mutate(finfish_trawl, survey_area = "Northeast Shelf"),
  grouping_vars = c("est_year", "season", "survey_area"),
  abundance_vals = "numlen_adj",
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





####  Length Spectra for Wigley Species - By Region  ####

# Wigley species length spectra
length_binspectra_wigley <- group_binspecies_spectra(
  ss_input = trawl_wigley,
  grouping_vars = c("est_year", "season", "survey_area"),
  abundance_vals = "numlen_adj",
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

# # Consistent with what we got last time? yes, cool
# lbw_save <- read_csv(here::here("Data/processed/wigley_species_length_spectra.csv"))
# lbw_save %>% 
#   mutate(yr_num = as.numeric(as.character(est_year))) %>% 
#   ggplot(aes(yr_num, b, color = season)) +
#   geom_point(size = 1, alpha = 0.6) +
#   #geom_ma(aes(linetype = "5-Year Moving Average"),n = 5, ma_fun = SMA) +
#   geom_smooth(method = "lm", linewidth = 1, se = F, 
#               aes(linetype = "Regression Fit")) +
#   facet_wrap(~survey_area) +
#   scale_color_gmri() +
#   labs(title = "Length Spectra Slope - MLE Bins Method Wigley - saved",
#        subtitle = "Enforced xmin = 1, xmax = max(length_cm + 1)",
#        y = "b",
#        x = "Year",
#        color = "Season")



####  Biomass Spectra for Wigley Species - By Region  ####


# For the Wigley species we can estimate individual weight from length
# then estimate size spectra using individual weights
bodymass_binspectra_wigley <- group_binspecies_spectra(
  ss_input = trawl_wigley,
  grouping_vars = c("est_year", "season", "survey_area"),
  abundance_vals = "numlen_adj",
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


####  Full Shelf Mass Spectra  ####

# For the Wigley species we can estimate individual weight from length
# then estimate size spectra using individual weights
bodymass_binspectra_wigley_shelf <- group_binspecies_spectra(
  ss_input = mutate(trawl_wigley, survey_area = "Northeast Shelf"),
  grouping_vars = c("est_year", "season", "survey_area"),
  abundance_vals = "numlen_adj",
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













# Plot comparison

bind_rows(
  mutate(fine_params, flag = "Okay"),
  mutate(weird_params, flag = "Off")) %>% 
ggplot(aes(est_year, val, color = flag)) +
  geom_point(alpha = 0.5) +
  facet_grid(var ~ survey_area, scales = "free_y")









####  Save Length/Mass Spectra  ####

# Save regional
write_csv(length_binspectra, here::here("Data/processed/finfish_length_spectra.csv"))
write_csv(length_binspectra_wigley, here::here("Data/processed/wigley_species_length_spectra.csv"))
write_csv(bodymass_binspectra_wigley, here::here("Data/processed/wigley_species_bodymass_spectra.csv"))

# Save Shelf
write_csv(length_binspectra_shelf, here::here("Data/processed/shelfwide_finfish_length_spectra.csv"))
write_csv(bodymass_binspectra_wigley_shelf, here::here("Data/processed/shelfwide_wigley_species_bodymass_spectra.csv"))




#### Setting Minimum Sizes at Abundance Peaks  ####

# In the quarto doc here:Chasing_Peasks.qmd
# Followed recommendations to set minimum size for distribution
# based on location of normalized abundance peak

# This is the key for the peaks and the minimum size to use:


# This is code to run all of those:


# Make the key
min_size_key <- read_csv(here::here("Data/processed/wigley_species_l2peaks_key.csv"))



# # Join the data back to the original
# wigley_in_new <- wigley_in %>%
#   mutate(area = "Northeast Shelf") %>%
#   bind_rows(wigley_in) %>%
#   left_join(min_size_key, join_by("est_year", "area", "season")) %>%
#   filter(wmin_g >= min_cutoff_g)
# 
# 
# 
# # Re-Run the MLE method to get new data
# 
# # and again for 16g to 50kg
# peak_chase_results <- group_binspecies_spectra(
#     ss_input       = wigley_in_new,
#     grouping_vars  = c("est_year", "season", "survey_area"),
#     abundance_vals = "numlen_adj",
#     length_vals    = "length_cm",
#     use_weight     = TRUE,
#     isd_xmin       = NULL,
#     global_min     = FALSE,
#     isd_xmax       = NULL,
#     global_max     = FALSE,
#     bin_width      = 1,
#     vdiff          = 2) %>% 
#   mutate(est_year = as.numeric(est_year)) %>% 
#   left_join(area_df)
# 
# 
# shelf_peak_results  <- group_binspecies_spectra(
#     ss_input       = filter(wigley_in_new, area == "Northeast Shelf"),
#     grouping_vars  = c("est_year", "season"),
#     abundance_vals = "numlen_adj",
#     length_vals    = "length_cm",
#     use_weight     = TRUE,
#     isd_xmin       = NULL,
#     global_min     = FALSE,
#     isd_xmax       = NULL,
#     global_max     = FALSE,
#     bin_width      = 1,
#     vdiff          = 2) %>% 
#   mutate(
#     est_year = as.numeric(est_year),
#     area = "Northeast Shelf") %>% 
#   left_join(area_df)
# 
# 
# # Join those together for one file
# moving_peak_spectra <- bind_rows(peak_chase_results, shelf_peak_results)


# # Save them
# write_csv(moving_peak_spectra, here::here("Data/processed/wigley_species_l2peaks_bmspectra.csv"))








#####  New Minimum Size Cutoffs?  ####



# Full shelf 10g cutoff
# For the Wigley species we can estimate individual weight from length
# then estimate size spectra using individual weights
bodymass_binspectra_wigley_shelf_100 <- group_binspecies_spectra(
  ss_input = mutate(trawl_wigley, survey_area = "Northeast Shelf") %>% 
    filter(wmin_g >= 100, wmax_g <=10000),
  grouping_vars = c("est_year", "season", "survey_area"),
  abundance_vals = "numlen_adj",
  length_vals = "length_cm",
  use_weight = TRUE,
  isd_xmin = 1,
  global_min = TRUE,
  isd_xmax = NULL,
  global_max = FALSE,
  bin_width = 1,
  vdiff = 2)



# And do a plot check
bodymass_binspectra_wigley_shelf_100 %>% 
  filter(season %in% c("Spring", "Fall")) %>% 
  mutate(yr_num = as.numeric(as.character(est_year))) %>% 
  ggplot(aes(yr_num, b, color = fct_rev(season))) +
  geom_point(size = 1, alpha = 0.6) +
  geom_ma(aes(linetype = "5-Year Moving Average"), n = 5, ma_fun = SMA) +
  geom_smooth(
    method = "lm", linewidth = 1, se = F, 
    aes(group = 1, linetype = "Regression Fit")) +
  facet_wrap(~survey_area) +
  scale_color_gmri() +
  labs(title = "Bodymass Spectra - MLE Bins Method Wigley",
       subtitle = "Enforced xmin = 100g, xmax = max(ind_weight_g), mass >10kg removed",
       y = "b",
       x = "Year",
       color = "Season")




\
