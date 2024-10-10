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


####  Processing Function  ####


# debugging the function that does it for groups
group_size_metrics <- function(
    size_data, 
    .group_cols = "Year", 
    .abund_col = "numlen_adj",
    .length_col = "length_cm",
    has_weights = FALSE,
    .weight_col = "ind_weight_kg"){
  
  
  # 1. Build group_level from desired group columns
  group_size_data <- size_data %>% 
    drop_na({{.abund_col}}, {{.length_col}}) %>% 
    unite(
      col = "group_var", 
      {{.group_cols}}, 
      sep = "-", 
      remove = TRUE, 
      na.rm = TRUE)
  
  
  # Run Min/Max/Avg. Size for the group
  group_results <- group_size_data %>% 
    split(.$group_var) %>% 
    imap_dfr(function(group_data, group_name){
      
      # Length (measured for all)
      mean_len    <- weighted.mean(
        group_data[, .length_col], 
        group_data[, .abund_col], 
        na.rm = T)
      med_len   <- matrixStats::weightedMedian(
        group_data[, .length_col], 
        group_data[, .abund_col], 
        na.rm = T)
      min_len     <- min(group_data[, .length_col], na.rm = T)
      max_len     <- max(group_data[, .length_col], na.rm = T)
      
      
      
      # Weights (estimated from length w/ wigley data)
      if(has_weights == TRUE){
        group_data <- drop_na(group_data, {{.weight_col}})
        
        mean_weight <- weighted.mean(
          group_data[,.weight_col], 
          group_data[, .abund_col], 
          na.rm = T)
        med_weight  <- matrixStats::weightedMedian(
          group_data[, .weight_col], 
          group_data[, .abund_col], 
          na.rm = T)
        min_weight  <- min(group_data[, .weight_col], na.rm = T)
        max_weight  <- max(group_data[, .weight_col], na.rm = T)
      }
      
      # Abundance totals
      total_abund <- sum(group_data[, .abund_col], na.rm = T)
      
      # Total number of species
      num_species <- group_data %>% 
        filter(is.na(comname) == FALSE) %>% 
        distinct(comname) %>% 
        nrow()
      
      # Put in table
      table_out <- data.frame(
        "group_var"        = group_name,
        "n_species"        = num_species,  
        "numlen_adj_total" = total_abund,
        "mean_len_cm"      = mean_len,
        "med_len_cm"       = med_len,
        "min_len_cm"       = min_len,
        "max_len_cm"       = max_len)
      
      # Add weight information
      if(has_weights == TRUE){
        table_out <- bind_cols(
          table_out,
          data.frame(
            "mean_wt_kg"       = mean_weight,
            "med_wt_kg"        = med_weight,
            "min_wt_kg"        = min_weight,
            "max_wt_kg"        = max_weight)) %>% 
          
          # Replace Inf with NA
          mutate(across(
            .cols = numlen_adj_total:max_wt_kg, 
            .fns = ~ifelse(is.infinite(abs(.x)), NA, .x)))
        
      }
      
      
      
      # return the table
      return(table_out)
    }) %>% 
    # Decompose the group ID back into original columns
    separate(
      group_var, 
      sep = "-", 
      into = .group_cols, 
      remove = T)
  
  
  
  # Return the results
  return(group_results)
  
  
}



####  Estimate Mean/Median Body Size  ####


#### Finfish Community  ####



# Lengths only for broader finfish community
finfish_lengths <- finfish_trawl %>% 
  group_size_metrics(
    size_data = ., 
    .abund_col = "numlen_adj", 
    .group_cols = c("survey_area", "est_year", "season"),
    has_weights = FALSE)


finfish_lengths_shelf <- trawl_wigley %>% 
  mutate(survey_area = "Northeast Shelf") %>% 
  group_size_metrics(
    size_data = ., 
    .abund_col = "numlen_adj", 
    .group_cols = c("survey_area", "est_year", "season"),
    has_weights = FALSE)



#### Wigley Species  ####

# Length and Weights for Wigley Species 
wigley_sizes <- trawl_wigley %>% 
  group_size_metrics(
    size_data = ., 
    .abund_col = "numlen_adj", 
    .group_cols = c("survey_area", "est_year", "season"),
    has_weights = TRUE, 
    .weight_col = "ind_weight_kg")

wigley_sizes_shelf <- trawl_wigley %>% 
  mutate(survey_area = "Northeast Shelf") %>% 
  group_size_metrics(
    size_data = ., 
    .abund_col = "numlen_adj", 
    .group_cols = c("survey_area", "est_year", "season"),
    has_weights = TRUE, 
    .weight_col = "ind_weight_kg")




####  Plot Results  ####

# Finfish Community median length
finfish_trawl %>% distinct(comname) %>% nrow()
bind_rows(finfish_lengths, finfish_lengths_shelf) %>% 
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
bind_rows(wigley_sizes, wigley_sizes_shelf) %>% 
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
write_csv(finfish_lengths_shelf, here::here("Data/processed/shelfwide_finfish_species_length_summary.csv"))
write_csv(wigley_sizes_shelf, here::here("Data/processed/shelfwide_wigley_species_size_summary.csv"))

