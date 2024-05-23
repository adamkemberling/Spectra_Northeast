####  Landings Data EDA  ####
# 6/7/2022



####  Packages  ####
{
  library(readxl)
  library(here)
  library(tidyverse)
  library(gmRi)
  library(janitor)
  library(gt)
  library(sf)
  library(targets)
  library(scales)
}


####  Data  ####

# Different Sheets

# # Organized by port
# by_port <- read_xlsx(
#   path = here("Data/KMills_landings by area 1964-2021 - FINFISH ONLY_MAY 2022.xlsx"), 
#   sheet = 1, 
#   skip = 0) %>% 
#   clean_names()
# 
# # Organized by stat zone
# by_szone <- read_xlsx(
#   path = here("Data/KMills_landings by area 1964-2021 - FINFISH ONLY_MAY 2022.xlsx"), 
#   sheet = 2, 
#   skip = 0) %>% 
#   clean_names()
# 
# # List of species
# spec_list <- read_xlsx(
#   path = here("Data/KMills_landings by area 1964-2021 - FINFISH ONLY_MAY 2022.xlsx"), 
#   sheet = 3, 
#   skip = 0) %>% 
#   clean_names()
# 
# 
# # Landings of finfish* sheet 5
# res_path <- cs_path("res")
# landings <- read_xlsx(
#   path = here("Data/KMills_landings by area 1964-2021_JUN 2022.xlsx"), sheet = 5) %>% 
#   rename_all(tolower)




####  Assign Regions to Follow Strata-Aggregates  ####

landings <- read_xlsx(
  path = here("Data/KMills_landings by area 1964-2021_JUN 2022.xlsx"), sheet = 5) %>% 
  rename_all(tolower)




# Make a list of zones to roughly match the survey areas:
fish_zones <- list(
  "Gulf of Maine"        = c(511:515, 464, 465),
  "Georges Bank"         = c(521, 522, 525, 561, 562),
  "Southern New England" = c(611, 612, 613, 616, 526, 537, 538, 539),
  "Mid-Atlantic Bight"   = c(614:615, 621, 622, 625, 626, 631, 632))



# # Why no GB Landings in 2011?
# landings %>% filter(
#   `stat area` %in% fish_zones$"Georges Bank", year == 2011)




# Add the labels into the landings data and remove what we don't need there:
landings <- landings %>% 
  mutate(
    survey_area = case_when(
      `stat area` %in% fish_zones$"Gulf of Maine" ~ "Gulf of Maine",
      `stat area` %in% fish_zones$"Georges Bank" ~ "Georges Bank",
      `stat area` %in% fish_zones$"Southern New England" ~ "Southern New England",
      `stat area` %in% fish_zones$"Mid-Atlantic Bight" ~ "Mid-Atlantic Bight")) %>% 
  filter(survey_area %in% c("Georges Bank", "Gulf of Maine", "Southern New England", "Mid-Atlantic Bight"))




# Get Summaries
landings_summ <- landings %>% 
  rename(
    "weight_lb" = `landed lbs`,
    "live_lb"   = `live lbs`) %>% 
  group_by(year, survey_area) %>% 
  summarise( 
    across(
      .cols = c(value, weight_lb, live_lb), 
      .fns = list(mean = ~mean(.x , na.rm = T), 
                  total = ~sum(.x , na.rm = T)), 
      .names = "{.fn}_{.col}"), 
    .groups = "drop")  


# Save those
write_csv(landings_summ, here::here("Data/processed/GARFO_regional_finfish_landings.csv"))




