#### Tidy Trawl Data
# All basic cleaning before median sizes or spectra slopes are estimated



####  Packages  ####
####  Packages  ####
library(tidyverse)
library(gmRi)


# Load processing functions
source(here::here("Code/R/Processing_Functions.R"))


####  Load and Prep Data  ####

#Path where "survdat" & spp_keys.csv data exists
data_path <- here::here("Data/raw")


# Just read it in, no* species filtering yet
# General tidying only, removal of strata outside our study area 
# (removes inshore and rarely/inconsistently sampled strata etc.)
trawl_basic <- tidy_nefsc_trawl(data_path = data_path)


# Minor tidying:
# Filter to spring and fall only
#  Drop duplicate length * species *, why do these exist? 34 rows duplicated - dropped
trawl_basic <- trawl_basic %>% 
  filter(season %in% c("Spring", "Fall"),
         est_year >= 1970, 
         est_year <= 2019) %>% 
  distinct(id, comname, catchsex, length_cm, .keep_all = T) 




# So if the species is a decapod or shellfish then they estimate length to the nearest mm. 
# Otherwise length is to the nearest cm. So this is important to distinguish for lobster vs. finfish predators.
crusts <- trawl_basic %>%
  mutate(length_char = as.character(length_cm)) %>%
  filter(grepl("[.]", length_char)) %>%
  distinct(comname) %>% 
  pull(comname)


# Drop crustaceans
# Wigley is all fish, not necessary
finfish_trawl <- filter(trawl_basic, comname %not in% crusts)


####  Load Wigley Species  ####

# Create a second dataset containing
# species with l-w coeffficients from wigley 06
trawl_wigley <- add_wigley_lw(trawl_basic, data_path = data_path)


# Get weight at length, and weight at length + 1
trawl_wigley <- trawl_wigley %>% 
  mutate(
    lngt_max     = length_cm + 1,              # maximum length for a fish of wmin
    ln_wmax_kg   = (ln_a + b * log(lngt_max)), # ln(max weight for fish of that wmin
    ind_weight_g = ind_weight_kg * 1000,
    wmin_g       = ind_weight_g,                # minimum weight of individual in grams
    wmax_g       = exp(ln_wmax_kg) * 1000       # max weight of individual in grams
  )    






####  Save Tidy Trawl Datasets  ####
write_csv(trawl_basic, here::here("Data/processed/trawl_clean_all_data.csv"))
write_csv(finfish_trawl, here::here("Data/processed/finfish_trawl_data.csv"))
write_csv(trawl_wigley, here::here("Data/processed/wigley_species_trawl_data.csv"))

