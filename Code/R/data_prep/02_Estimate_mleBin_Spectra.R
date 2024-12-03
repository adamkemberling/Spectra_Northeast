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


# # Test 1 subgroup
# finfish_trawl %>% 
#   filter(est_year == 1970,
#          survey_area == "GB",
#          season == "Fall") %>% 
#   # run the function
#   group_binspecies_spectra(
#     ss_input = .,
#     grouping_vars = c("est_year", "season", "survey_area"),
#     abundance_vals = "numlen_adj",
#     length_vals = "length_cm",
#     use_weight = FALSE,
#     global_min = TRUE,
#     isd_xmin = 1,
#     isd_xmax = NULL,
#     global_max = FALSE,
#     bin_width = 1, 
#     vdiff = 4)



##### All Species  ####
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





#####  Wigley Species   ####

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




####  Mass Spectra   ####


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
  abundance_vals = "numlen_adj",
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


# And do a plot check
mle_results_16g %>% 
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
  facet_wrap(~survey_area, ncol = 1, scales = "free") +
  scale_color_gmri() +
  labs(title = "Bodymass Spectra - MLE Bins Method Wigley",
       subtitle = "Floating xmin = 1g, xmax = max(ind_weight_g)",
       y = "b",
       x = "Year",
       color = "Season")






# # Save those
# write_csv(
#   mle_results_16g,
#   here::here("Data/processed/wigley_species_min16_bodymass_spectra.csv"))



####__________________________####
#### Dynamic WMIN Locations  ####

# In the quarto doc here:Chasing_Peasks.qmd
# Followed recommendations to set minimum size for distribution
# based on location of normalized abundance peak

# This is the key for the peaks and the minimum size to use:
# this was made in chasing_peaks.qmd
min_size_key <- read_csv(here::here("Data/processed/wigley_species_l2peaks_key.csv"))

# This is code to run all of those:

# Join the data back to the original so we can run the shelf at the same time
# Join the data back to the original
wigley_in_new <- trawl_wigley %>%
  mutate(area = "Northeast Shelf", area_titles = "Northeast Shelf", survey_area = "Northeast Shelf") %>%
  bind_rows(left_join(trawl_wigley, area_df)) %>%
  left_join(min_size_key, join_by("est_year", "area", "season")) %>%
  filter(wmin_g >= min_cutoff_g)






# Re-Run the MLE method to get new data with those cutoffs
# drops survey_areas not in the pre-ordained list
peak_chase_results <- group_binspecies_spectra(
    ss_input       = wigley_in_new,
    grouping_vars  = c("est_year", "season", "survey_area"),
    abundance_vals = "numlen_adj",
    length_vals    = "length_cm",
    use_weight     = TRUE,
    isd_xmin       = NULL,
    global_min     = FALSE,
    isd_xmax       = NULL,
    global_max     = FALSE,
    bin_width      = 1,
    vdiff          = 2) %>%
  mutate(est_year = as.numeric(est_year)) %>%
  left_join(area_df)

# do we still have survey area etc.
peak_chase_results %>% distinct(survey_area, area, area_titles)



# And do a plot check
peak_chase_results %>% 
  filter(season %in% c("Spring", "Fall")) %>% 
  mutate(
    yr_num = as.numeric(as.character(est_year)),
    survey_area = factor(
      survey_area, 
      levels = c("Northeast Shelf", "GoM", "GB", "SNE", "MAB"))) %>% 
  
  ggplot(aes(yr_num, b, color = season)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_ma(aes(linetype = "5-Year Moving Average"),n = 5, ma_fun = SMA) +
  geom_smooth(method = "lm", linewidth = 1, se = F, 
              aes(linetype = "Regression Fit")) +
  facet_wrap(~survey_area, ncol = 1, scales = "free") +
  scale_color_gmri() +
  labs(title = "Bodymass Spectra - MLE Bins Method Wigley",
       subtitle = "Floating xmin = 1g, xmax = max(ind_weight_g)",
       y = "b",
       x = "Year",
       color = "Season")




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






####__________________________####


####  Saving Spectra Exponent Results  ####

# Save Regional subsets
write_csv(length_binspectra, here::here("Data/processed/finfish_length_spectra.csv"))
write_csv(length_binspectra_wigley, here::here("Data/processed/wigley_species_length_spectra.csv"))
write_csv(bodymass_binspectra_wigley, here::here("Data/processed/wigley_species_bodymass_spectra.csv"))

# Save Shelfwide subsets
write_csv(length_binspectra_shelf, here::here("Data/processed/shelfwide_finfish_length_spectra.csv"))
write_csv(bodymass_binspectra_wigley_shelf, here::here("Data/processed/shelfwide_wigley_species_bodymass_spectra.csv"))

# Saving Shifted Wmin Spectra:

# dynamic wmin spectra:
write_csv(peak_chase_results, here::here("Data/processed/wigley_species_l2peaks_bmspectra.csv"))

# Save the data used for shifting peak fits
write_csv(wigley_in_new, here::here("Data/processed/wigley_species_trawl_wmin_filtered.csv"))

# universal 16g wmin:
write_csv(mle_results_16g, here::here("Data/processed/wigley_species_min16_bodymass_spectra.csv"))
