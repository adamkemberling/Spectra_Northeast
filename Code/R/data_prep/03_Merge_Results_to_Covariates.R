#### Joining Length/Biomass Spectra with Temp + Landings


####  Packages  ####
library(gmRi)
library(tidyverse)




# Larger finfish community
finfish_length_spectra <- read_csv(here::here("Data/processed/finfish_length_spectra.csv"))
finfish_sizes <- read_csv(here::here("Data/processed/finfish_species_length_summary.csv"))

# Wigley Species Community
wigley_sizes <- read_csv(here::here("Data/processed/wigley_species_size_summary.csv"))
wigley_length_spectra <- read_csv(here::here("Data/processed/wigley_species_length_spectra.csv"))
wigley_bodymass_spectra <- read_csv(here::here("Data/processed/wigley_species_bodymass_spectra.csv"))

# Shifting wmin
wigley_shifting_spectra <- read_csv(here::here("Data/processed/wigley_species_l2peaks_bmspectra.csv"))



####  Covariate Data  ####

# vectors for factor levels
area_levels <- c("Northeast Shelf", "GoM", "GB", "SNE", "MAB")
area_levels_long <- c("Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight")

# table to join for swapping shorthand for long-hand names
area_df <- data.frame(
  area = c("Scotian Shelf", "Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight", "Northeast Shelf"),
  survey_area = c("SS", "GoM", "GB", "SNE", "MAB", "Northeast Shelf"),
  area_titles = c("Scotian Shelf", "Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight", "Northeast Shelf"))




####  1. GARFO Landings  ####


# Read in GARFO Landings, averaged over ~survey_areas

# # From 2022 GMRI inventory
# landings_annual <- read_csv(here::here("Data/processed/GARFO_regional_finfish_landings.csv")) %>% 

# From Andy Beet
landings_annual <- read_csv(here::here("Data/processed/BEET_GARFO_regional_finfish_landings.csv")) %>% 
  rename(area = survey_area,
         est_year = year) %>% 
  left_join(area_df) %>% 
  select(survey_area, area, area_titles, est_year, 
         total_weight_lb = total_live_lb)

# Get shelf totals
landings_annual <- landings_annual %>% 
  group_by(est_year) %>% 
  summarise(total_weight_lb = sum(total_weight_lb, na.rm = T)) %>%
  mutate(
    survey_area = "Northeast Shelf",
    area = "Northeast Shelf",
    area_titles = "Northeast Shelf") %>% 
  bind_rows(landings_annual)







####  2. Seasonal Bottom Temperatures  ####

# Read in du pontavice bottom temperatures, 
# these are averaged within survey_areas in Code/py/Annual_BT_Processing.ipynb

# Timeseries are produced regionally in: 
# 00_Bottom_Temp_Processing.R
# 01_Average_Seasonal_BT.R

# ## Annual Bottom temp
# bot_temps <- read_csv(here::here("Data/processed", "trawl_region_bottom_temps.csv")) %>%
#   rename(area = survey_area,
#          est_year = year) %>%
#   left_join(area_df) %>%
#  select(survey_area, area_titles, est_year, bot_temp) 

# Seasonal Bottom Temp
bot_temps <- read_csv(here::here("Data/processed", "trawl_region_seasonal_bottom_temps.csv")) %>% 
  rename(area = survey_area,
         est_year = year) %>% 
  left_join(area_df) %>% 
  select(survey_area, area_titles, est_year, season, bot_temp, overall_avg_btemp, bot_temp_anom)




####  3. Zooplankton Community Indices  ####

# Load zooplankton community size indices from {ecodata}

# Small large index not present anymore, has "regimes"
zoo_regimes <- ecodata::zoo_regime
forage_index <- ecodata::forage_index
ecodata::plot_forage_anomaly(report = "NewEngland")









#### Create Modeling Datasets  ####
# These weren't intended to be modeled at shelf scale, so shelf gets dropped here

# 1. Large Community Median Length + length spectra
ffish_medlen_model_df     <- left_join(finfish_sizes, bot_temps) %>% left_join(landings_annual)
ffish_lenspectra_model_df <- left_join(finfish_length_spectra, bot_temps) %>% left_join(landings_annual)


# 2. Smaller Community Length/Weight/Spectra
wigley_medlen_model_df     <- left_join(wigley_sizes, bot_temps) %>% left_join(landings_annual)
wigley_lenspectra_model_df <- left_join(wigley_length_spectra, bot_temps) %>% left_join(landings_annual)
wigley_bmspectra_model_df  <- left_join(wigley_bodymass_spectra, bot_temps) %>% left_join(landings_annual)


# 3. Shifting WMIN bodymass spectra
wigley_shifting_bmspectra_model_df <- left_join(wigley_shifting_spectra, bot_temps) %>% left_join(landings_annual)


ggplot(wigley_shifting_bmspectra_model_df, aes(bot_temp, b, color = season)) + geom_point() + facet_wrap(~area)
ggplot(wigley_shifting_bmspectra_model_df, aes(bot_temp, total_weight_lb, color = season)) + geom_point() + facet_wrap(~survey_area)


#### Save Modeling Data  ####

# data limited community
write_csv(ffish_medlen_model_df, here::here("Data/model_ready/large_community_medlength_mod.csv"))
write_csv(ffish_lenspectra_model_df, here::here("Data/model_ready/large_community_lenspectra_mod.csv"))

# well studied community
write_csv(wigley_medlen_model_df, here::here("Data/model_ready/wigley_community_medsize_mod.csv"))
write_csv(wigley_lenspectra_model_df, here::here("Data/model_ready/wigley_community_lenspectra_mod.csv"))
write_csv(wigley_bmspectra_model_df, here::here("Data/model_ready/wigley_community_bmspectra_mod.csv"))

# shifting wmin
write_csv(wigley_shifting_bmspectra_model_df, here::here("Data/model_ready/wigley_community_shifting_bmspectra_mod.csv"))



