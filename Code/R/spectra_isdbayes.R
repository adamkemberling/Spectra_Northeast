####  Spectra Slope Verifications using ISDBayes ####



####  Packages  ####
library(tidyverse)
library(gmRi)
library(tidyquant)
library(isdbayes)
library(brms)
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
