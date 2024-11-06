####  Seasonal Bottom temperature Differences
#### Seasonal Bottom Temperature Anomalies


####  Packages  ####
library(tidyverse)
library(gmRi)



# Theme
theme_set(theme_gmri(rect = element_rect(fill = "white", color = NA)))

# vectors for factor levels
area_levels_long <- c("Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight")



####  Load Data

# trawl dta - to confirm dates
trawl <- read_csv(here::here("Data/processed/finfish_trawl_data.csv"))


# monthly bottom temperature
bpath <- cs_path("res", "Du_Pontavice_Combined_BT/RegionalTimeseries")
bt_monthly <- read_csv(str_c(bpath, "trawl_region_monthly_bottom_temps.csv"))



# When are the seasons
trawl %>% 
  split(.$season) %>% 
  map(function(x){
    drange <- pull(x, est_towdate) %>% range(na.rm = T)
    })

# Spring is March 12 - May - 12
# Fall is Sept 5 - November 12




#### Comparisons  ####

# Are the regions consistently different temperatures, or only in one season?
# month weighting so we don't have to go back to daily when averaging
month_wts <- data.frame(
    "month" =  c("03", "04", "05", "09", "10", "11"),
    "n_days" = days_in_month(
      as.Date(c("2000-03-01","2000-04-01","2000-05-01",
                "2000-09-01","2000-10-01","2000-11-01")))
    )


# Calculate average for the season
seasonal_bt <- bt_monthly %>% 
  left_join(month_wts) %>% 
  mutate(
    season = ifelse(month %in% str_pad(c(3:5), side = "left", pad = "0", width = 2), "Spring", NA),
    season = ifelse(month %in% str_pad(c(9:11), side = "left", pad = "0", width = 2), "Fall", season),
    survey_area = factor(survey_area, levels = c("Northeast Shelf", area_levels_long))) %>% 
  drop_na() %>% 
  group_by(survey_area, year, season) %>% 
  summarise(
    bot_temp = weighted.mean(bot_temp, w = n_days),
    .groups = "drop"
  )

# Plot Seasonal BT
seasonal_bt %>% 
  ggplot(aes(survey_area, bot_temp, color = season)) +
  geom_boxplot(position = position_dodge(width = 0.5))


# Which Months are coldest/hottest
bt_monthly %>% 
  mutate(survey_area = factor(survey_area, levels = c("Northeast Shelf", area_levels_long))) %>% 
  drop_na() %>% 
  ggplot(aes(survey_area, bot_temp, color = month)) +
  geom_boxplot(position = position_dodge(width = 0.5)) +
  facet_wrap(~survey_area, scales = "free")



####  Seasonal BT Anomalies  ####

# Don't pick a climatology, just do all years
bt_seas_clim <- seasonal_bt %>% 
  group_by(survey_area, season) %>% 
  summarise(overall_avg_btemp = mean(bot_temp),
            .groups = "drop") %>% 
  right_join(seasonal_bt) %>% 
  mutate(bot_temp_anom = bot_temp - overall_avg_btemp)




# #### Saving Seasonal Bottom Temperature  ####
write_csv(bt_seas_clim, here::here("Data/processed", "trawl_region_seasonal_bottom_temps.csv"))







