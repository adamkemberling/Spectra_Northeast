# Process Spectra using Northeast Trawl Survey


####  Packages  ####
library(tidyverse)
library(gmRi)
library(sizeSpectra)
library(tidyquant)


# Load processing functions
source(here::here("Code/R/Processing_Functions.R"))



####  Load + Tidy Data  ####

# Load data from local directory
#Path where "survdat" & spp_keys.csv data exists
data_path <- here::here("Data/")


# Just read it in, no* species filtering yet
# General tidying only, removal of strata outside our study area 
# (removes inshore and rarely/inconsistently sampled strata etc.)
trawl_basic <- tidy_nefsc_trawl(data_path = data_path)

# try removing no abundance or no biomass records
# we drop 17545 records dropping cases where positive abundanance/biomass
# has no value of biomass/abundance




# Minor tidying:
# Filter to spring and fall only
#  Drop duplicate length*species*, why do these exist? 34 rows duplicated - dropped
trawl_basic <- trawl_basic %>% 
  filter(season %in% c("Spring", "Fall")) %>% 
  distinct(id, comname, catchsex, length_cm, .keep_all = T) 


# Create a second dataset containing
# species with l-w coeffficients from wigley 06
trawl_wigley <- add_wigley_lw(trawl_basic, data_path = data_path)


# Get weight at length, and weight at length + 1
trawl_wigley <- trawl_wigley %>% 
  mutate(
    lngt_max     = length_cm + 1,    # maximum length for a fish of wmin
    ln_wmax      = (ln_a + b * log(lngt_max)),  # max weight for fish of that wmin
    ind_weight_g = ind_weight_kg * 1000,    # grams
    wmin_g       = ind_weight_g,               # minimum weight of individual in grams
    wmax_g       = exp(ln_wmax) * 1000         # max weight of individual in grams
  )    



####  Filter to Only Finfish Species for Analysis  ####

# --- quick approach, remove things measured in mm

# So if the species is a decapod or shellfish then they estimate length to the nearest mm. 
# Otherwise length is to the nearest cm. So this is important to distinguish for lobster vs. finfish predators.
crusts <- trawl_basic %>%
  mutate(length_char = as.character(length_cm)) %>%
  filter(grepl("[.]", length_char)) %>%
  distinct(comname) %>% 
  pull(comname)
crusts

# Drop crustaceans
# Wigley is all fish, not necessary
finfish_trawl <- filter(trawl_basic, comname %not in% crusts)


# OR be thorough -- join taxa info, subset finfish families





#### Estimate Length Spectra  ####




# All finfish length_spectra
length_binspectra <- group_binspecies_spectra(
  ss_input = trawl_wigley,
  grouping_vars = c("est_year", "season", "survey_area"),
  abundance_vals = "numlen_adj",
  length_vals = "length_cm",
  use_weight = FALSE,
  isd_xmin = 1,
  isd_xmax = NULL,
  global_min = TRUE,
  global_max = FALSE,
  bin_width = 1)



# Plot it
length_binspectra %>% 
  filter(season %in% c("Spring", "Fall")) %>% 
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




# Wigley species length spectra
length_binspectra_wigley <- group_binspecies_spectra(
  ss_input = trawl_wigley,
  grouping_vars = c("est_year", "season", "survey_area"),
  abundance_vals = "numlen_adj",
  length_vals = "length_cm",
  use_weight = FALSE,
  isd_xmin = 1,
  isd_xmax = NULL,
  global_min = TRUE,
  global_max = FALSE,
  bin_width = 1)



# Plot them
length_binspectra_wigley %>% 
  filter(season %in% c("Spring", "Fall")) %>% 
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





####  Save Length Spectra  ####
write_csv(length_binspectra, here::here("Data/processed/finfish_length_spectra.csv"))
write_csv(length_binspectra_wigley, here::here("Data/processed/wigley_species_length_spectra.csv"))




####  Biomass Spectra  ####
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
  #geom_ma(aes(linetype = "5-Year Moving Average"),n = 5, ma_fun = SMA) +
  geom_smooth(method = "lm", linewidth = 1, se = F, 
              aes(linetype = "Regression Fit")) +
  facet_wrap(~survey_area) +
  scale_color_gmri() +
  labs(title = "Bodymass Spectra - MLE Bins Method Wigley",
       subtitle = "Enforced xmin = 1g, xmax = max(ind_weight_g)",
       y = "b",
       x = "Year",
       color = "Season")





# Why are there so many -1's? - that's an issue
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
  separate(col = "group_var", into = c("est_year", "survey_area", "season"), sep = "-") %>% 
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


# Plot the data for ISD as density?
bind_rows(
  mutate(fine_inputs, flag = "Okay"),
  mutate(weird_inputs, flag = "Off")) %>% 
  group_by(flag, est_year, survey_area, season, comname, length_cm, wmin_g, wmax_g) %>% 
  summarise(
    Number = sum(numlen_adj),
    .groups = "drop") %>% 
  ggplot(aes(wmin_g, Number, color = flag)) +
  geom_point(alpha = 0.5) +
  scale_y_log10() +
  facet_grid(est_year ~ survey_area, scales = "free_y")
