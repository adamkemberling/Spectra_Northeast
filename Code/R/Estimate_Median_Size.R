####  Mean/Median Body Size


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


# Minor tidying:
# Filter to spring and fall only
#  Drop duplicate length*species*, why do these exist? 34 rows duplicated - dropped
trawl_basic <- trawl_basic%>% 
  filter(season %in% c("Spring", "Fall")) %>% 
  distinct(id, comname, catchsex, length_cm, .keep_all = T) 


# Create a second dataset containing
# species with l-w coeffficients from wigley 06
trawl_wigley <- add_wigley_lw(trawl_basic, data_path = data_path)


# Get weight at length, and weight at length + 1
trawl_wigley <- trawl_wigley %>% 
  mutate(
    lngt_max     = length_cm + 1,              # maximum length for a fish of wmin
    ln_wmax      = (ln_a + b * log(lngt_max)), # max weight for fish of that wmin
    ind_weight_g = ind_weight_kg * 1000     
  )    





####  Estimate Mean/Median Body Size  ####



# Get weighted mean lengths and weights using 
# tow-level and total stratified abundances
# Note: uses numlen because individuals were ID'd and measured
group_size_metrics <- function(
    size_data, 
    .group_cols = "Year", 
    weighting_col = "numlen_adj"){
  
  
  # 1. Build group_level from desired group columns
  group_size_data <- size_data %>% 
    unite(
      col = "group_var", 
      {{.group_cols}}, 
      sep = "-", 
      remove = FALSE, 
      na.rm = FALSE)
  
  
  # Run Min/Max/Avg. Size for the group
  group_results <- group_size_data %>% 
    split(.$group_var) %>% 
    imap_dfr(function(group_data, group_name){
      
      # Length (measured for all)
      mean_len    <- weighted.mean(group_data[, "length_cm"], group_data[, weighting_col], na.rm = T)
      med_len     <- matrixStats::weightedMedian(group_data[, "length_cm"], group_data[, weighting_col], na.rm = T)
      min_len     <- min(group_data[, "length_cm"], na.rm = T)
      max_len     <- max(group_data[, "length_cm"], na.rm = T)
      
      # Weights (estimated from length w/ wigley data)
      mean_weight <- weighted.mean(group_data[,"ind_weight_kg"], group_data[, weighting_col], na.rm = T)
      med_weight  <- matrixStats::weightedMedian(group_data[, "ind_weight_kg"], group_data[, weighting_col], na.rm = T)
      min_weight  <- min(group_data[, "ind_weight_kg"], na.rm = T)
      max_weight  <- max(group_data[, "ind_weight_kg"], na.rm = T)
      
      # Abundance totals
      total_abund <- sum(group_data[, "numlen_adj"], na.rm = T)
      strat_abund <- sum(group_data[, "strat_total_abund_s"], na.rm = T)
      
      # number of species
      num_species <- group_data %>% 
        filter(is.na(comname) == FALSE) %>% 
        distinct(comname) %>% 
        length()
      
      # Put in table
      table_out <- data.frame(
        "group_var"    = group_name,
        "n_species"    = num_species,  
        "survey_abund" = total_abund,
        "strat_abund"  = strat_abund,
        "mean_len_cm"  = mean_len,
        "med_len_cm"   = med_len,
        "min_len_cm"   = min_len,
        "max_len_cm"   = max_len,
        "mean_wt_kg"   = mean_weight,
        "med_wt_kg"    = med_weight,
        "min_wt_kg"    = min_weight,
        "max_wt_kg"    = max_weight) %>% 
        
        # replace Inf with NA
        mutate(across(.cols = survey_abund:max_wt_kg, .fns = ~ifelse(is.infinite(abs(.x)), NA, .x)))
      
      # return the table
      return(table_out)
    })
  

  
  # Return the results
  return(group_results)
  
  
}