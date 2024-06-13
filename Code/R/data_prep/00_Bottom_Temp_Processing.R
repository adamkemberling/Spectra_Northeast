####  Bottom Temperature Regional Timeseries
# Prep the dupontavice BT for use in models


library(raster)
library(sf)
library(tidyverse)
library(gmRi)



# Crop bottom temperatures to the survey regions and make yearly



# Load polygons for aggregate regions based on trawl survey strata
regions <- read_sf(here::here("Data/raw", "nmfs_trawl_regions_collection.geojson"))



# load BT Netcdf - Obtained directly from Hubert du Pontavice via email
# This is Jan 1 Temperatures only...
ponta <- raster::stack(here::here("Data/raw","annual_bottom_temp_combined_product_1959_2020.nc"))







# Masking Function to clip to the study area
mask_nc <- function(ras_obj, mask_shape, rotate = TRUE){
  
  # First Mask using the polygon, rotate if needed
  if(rotate == TRUE){
    m1 <- raster::mask(raster::rotate(ras_obj), mask_shape)
  } else {
    m1 <- raster::mask(ras_obj, mask_shape)
  }
  
  # Then Crop
  m1 <- crop(m1, mask_shape)
  return(m1)
}


# Run it for each of them


# Prepare the regions as single polygons:
gb <- filter(regions, area == "Georges Bank") %>% st_union() %>% st_as_sf()
gom <- filter(regions, area == "Gulf of Maine") %>% st_union() %>% st_as_sf()
mab <- filter(regions, area == "Mid-Atlantic Bight") %>% st_union() %>% st_as_sf()
sne <- filter(regions, area == "Southern New England") %>% st_union() %>% st_as_sf()


# Use those polygons to mask
gb_bt <- raster::mask(ponta, gb) %>% crop(., gb) %>% cellStats(mean) %>% as.data.frame() %>% rownames_to_column("year")
gom_bt <- raster::mask(ponta, gom) %>% crop(., gom) %>% cellStats(mean) %>% as.data.frame() %>% rownames_to_column("year")
mab_bt <- raster::mask(ponta, mab) %>% crop(., mab) %>% cellStats(mean) %>% as.data.frame() %>% rownames_to_column("year")
sne_bt <- raster::mask(ponta, sne) %>% crop(., sne) %>% cellStats(mean) %>% as.data.frame() %>% rownames_to_column("year")


# Put them in a list for using purrr::map etc. or bind_rows
all_bt <- list(
  "Georges Bank" = gb_bt,
  "Gulf of Maine" = gom_bt,
  "Southern New England" = sne_bt,
  "Mid-Atlantic Bight" = mab_bt) %>% 
  bind_rows(.id = "survey_area") %>% 
  setNames(c("survey_area", "year", "bot_temp")) %>% 
  mutate(year = str_sub(year, 2, 5),
         year = as.numeric(as.character(year)))


# Plot the annual trends
ggplot(all_bt, aes(year, bot_temp)) +
  geom_line(aes(color = survey_area))


# Save them out
write_csv(all_bt, here::here("Data/processed", "trawl_region_bottom_temps.csv"))



# Bonus: save them back to box:
bpath <- cs_path("res", "Du_Pontavice_Combined_BT/RegionalTimeseries")
write_csv(all_bt, str_c(bpath, "trawl_region_bottom_temps.csv"))
