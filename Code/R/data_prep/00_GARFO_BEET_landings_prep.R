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
  library(scales)
}


####  Data  ####

# Landings data obtained from Andy Beet with the NEFSC


####  Assign Regions to Follow Strata-Aggregates  ####

landings <- read_rds(here::here("Data/raw/LandingsByYearAreaSpecies.rds")) %>% 
  rename_all(tolower)




# Make a list of zones to roughly match the survey areas:
fish_zones <- list(
  "Gulf of Maine"        = c(511:515, 464, 465),
  "Georges Bank"         = c(521, 522, 525, 561, 562),
  "Southern New England" = c(611, 612, 613, 616, 526, 537, 538, 539),
  "Mid-Atlantic Bight"   = c(614:615, 621, 622, 625, 626, 631, 632))



# # Why no GB Landings in 2011? - none in the data
# landings %>% filter(
#   `stat area` %in% fish_zones$"Georges Bank", year == 2011)

# Add the labels into the landings data and remove what we don't need there:
landings <- landings %>% 
  mutate(
    survey_area = case_when(
      `area` %in% fish_zones$"Gulf of Maine" ~ "Gulf of Maine",
      `area` %in% fish_zones$"Georges Bank" ~ "Georges Bank",
      `area` %in% fish_zones$"Southern New England" ~ "Southern New England",
      `area` %in% fish_zones$"Mid-Atlantic Bight" ~ "Mid-Atlantic Bight")) %>% 
  filter(survey_area %in% c("Georges Bank", "Gulf of Maine", "Southern New England", "Mid-Atlantic Bight"))


# Map them to be sure

# # Shapefiles for the fisheries stat zones
# res_path <- cs_path("res")
# stat_zones <- read_sf(str_c(res_path, "Shapefiles/Statistical_Areas/Statistical_Areas_2010_withNames.shp"))

# stat_zones <- stat_zones %>% 
#   mutate(
#     survey_area = case_when(
#       `Id` %in% fish_zones$"Gulf of Maine" ~ "Gulf of Maine",
#       `Id` %in% fish_zones$"Georges Bank" ~ "Georges Bank",
#       `Id` %in% fish_zones$"Southern New England" ~ "Southern New England",
#       `Id` %in% fish_zones$"Mid-Atlantic Bight" ~ "Mid-Atlantic Bight")) %>% 
#   filter(survey_area %in% c("Georges Bank", "Gulf of Maine", "Southern New England", "Mid-Atlantic Bight"))
# 
# ggplot(stat_zones) +
#   geom_sf(aes(fill = survey_area))
# write_sf(stat_zones, here::here("Data/raw/nmfs_statareas_collection.geojson"))



####  Remove Freshwater Species if they're there  ####
landings <- landings %>% mutate(sppname = tolower(common_name)) 
landings %>% distinct(sppname) %>% pull() %>% sort()

# Filter the freshwater species out
landings <- landings %>% 
  filter(
    !str_detect(sppname, "catfish"),
    !str_detect(sppname, "carp"),
    !str_detect(sppname, "crappie"),
    !str_detect(sppname, "snakehead"),
    sppname != "perch, white",
    sppname != "perch, yellow",
    sppname != "starfish",
    sppname != "sunfish"
  )




#### Create Annual Summaries  ####

# Get Summaries
landings_summ <- landings %>% 
  mutate("live_lb"   = mtlive * 2204.62) %>% 
  group_by(year, survey_area) %>% 
  summarise( 
    across(
      .cols = c(mtlive, live_lb), 
      .fns = list(mean = ~mean(.x , na.rm = T), 
                  total = ~sum(.x , na.rm = T)), 
      .names = "{.fn}_{.col}"), 
    .groups = "drop")  


# Save those
write_csv(landings_summ, here::here("Data/processed/BEET_GARFO_regional_finfish_landings.csv"))





####  Compare  ####
kmills_garfo <- read_csv(here::here("Data/processed/GARFO_regional_finfish_landings.csv"))


ggplot() +
  geom_line(
    data = landings_summ, aes(year, total_live_lb, color = "Andy Landings")) +
  geom_line(
    data = kmills_garfo, aes(year, total_live_lb, color = "Old Landings")) +
  facet_wrap(~survey_area)



