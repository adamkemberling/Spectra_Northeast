# Process Spectra using Northeast Trawl Survey


####  Packages  ####
library(tidyverse)
library(gmRi)
library(sizeSpectra)
library(tidyquant)
library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)


# Load processing functions
source(here::here("Code/R/Processing_Functions.R"))



####  Load Tidy Trawl Data  ####
finfish_trawl <- read_csv(here::here("Data/processed/finfish_trawl_data.csv"))
trawl_wigley <- read_csv(here::here("Data/processed/wigley_species_trawl_data.csv"))



#### Length Spectra All Species  ####


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


####  Length Spectra for Wigley Species  ####

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



####  Biomass Spectra for Wigley Species  ####
# Not functional here 5/23/2024

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


# Why are there so many -1's? - that's an issue, vecdiff?
weird_fits <- bodymass_binspectra_wigley %>% 
  filter(between(b, -1.005, -0.9995)) %>% 
  distinct(est_year, survey_area, season) %>% 
  mutate(est_year = as.numeric(est_year))


# Look at the size distributions for those years
weird_inputs <- trawl_wigley %>% inner_join(weird_fits)
fine_inputs <- trawl_wigley %>% anti_join(weird_fits)



weird_params <- weird_inputs %>% 
  unite("group_var", est_year, survey_area, season, sep = "-") %>% 
  split(.$group_var) %>% 
  map_dfr(function(x){
    data.frame(
      min_weight = min(x$wmin_g),
      max_weight = max(x$wmin_g),
      n = sum(x$numlen_adj)
    )
  }, .id = "group_var") %>% 
  separate(
    col = "group_var", 
    into = c("est_year", "survey_area", "season"), sep = "-") %>% 
  mutate(est_year = as.numeric(est_year)) %>% 
  pivot_longer(cols = c(min_weight, max_weight, n), names_to = "var", values_to = "val")

fine_params <- fine_inputs %>% 
  unite("group_var", est_year, survey_area, season, sep = "-") %>% 
  split(.$group_var) %>% 
  map_dfr(function(x){
    data.frame(
      min_weight = min(x$wmin_g),
      max_weight = max(x$wmin_g),
      n = sum(x$numlen_adj)
    )
  }, .id = "group_var") %>% 
  separate(col = "group_var", into = c("est_year", "survey_area", "season"), sep = "-") %>% 
  mutate(est_year = as.numeric(est_year)) %>% 
  pivot_longer(cols = c(min_weight, max_weight, n), names_to = "var", values_to = "val")




# Plot comparison
library(patchwork)
bind_rows(
  mutate(fine_params, flag = "Okay"),
  mutate(weird_params, flag = "Off")) %>% 
ggplot(aes(est_year, val, color = flag)) +
  geom_point(alpha = 0.5) +
  facet_grid(var ~ survey_area, scales = "free_y")









####  Save Length Spectra  ####
write_csv(length_binspectra, here::here("Data/processed/finfish_length_spectra.csv"))
write_csv(length_binspectra_wigley, here::here("Data/processed/wigley_species_length_spectra.csv"))
write_csv(bodymass_binspectra_wigley, here::here("Data/processed/wigley_species_bodymass_spectra.csv"))


