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





####  Covariate Data  ####

# vectors for factor levels
area_levels <- c("GoM", "GB", "SNE", "MAB")
area_levels_long <- c("Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight")

# table to join for swapping shorthand for long-hand names
area_df <- data.frame(
  area = c("Scotian Shelf", "Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight", "All"),
  survey_area = c("SS", "GoM", "GB", "SNE", "MAB", "Northeast Shelf"),
  area_titles = c("Scotian Shelf", "Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight", "Northeast Shelf"))


# Read in GARFO Landings, averaged over ~survey_areas
landings_annual <- read_csv(here::here("Data/processed/GARFO_regional_finfish_landings.csv")) %>% 
  rename(area = survey_area,
         est_year = year) %>% 
  left_join(area_df) %>% 
  select(survey_area, est_year, total_weight_lb)


# Read in du pontavice bottom temperatures, averaged within survey_areas
bot_temps <- read_csv(here::here("Data", "trawl_region_bottom_temps.csv")) %>% 
  rename(area = survey_area,
         est_year = year) %>% 
  left_join(area_df) %>% 
 select(survey_area, area_titles, est_year, bot_temp) %>% 
  mutate(
    survey_area = factor(
      survey_area, 
      levels = area_levels),
    region = str_replace_all(area_titles, "-| ", "_"),
    region = factor(
      region, 
      levels = area_levels_long))





#### Create Modeling Datasets  ####

# 1. Large Community Median Length + length spectra
ffish_medlen_model_df <- left_join(finfish_sizes, bot_temps) %>% left_join(landings_annual)
ffish_lenspectra_model_df <- left_join(finfish_length_spectra, bot_temps) %>% left_join(landings_annual)

# 2. Smaller Community Length/Weight/Spectra
wigley_medlen_model_df <- left_join(wigley_sizes, bot_temps) %>% left_join(landings_annual)
wigley_lenspectra_model_df <- left_join(wigley_length_spectra, bot_temps) %>% left_join(landings_annual)
wigley_bmspectra_model_df <- left_join(wigley_bodymass_spectra, bot_temps) %>% left_join(landings_annual)





####  LMER With Covariates  ####








#### Modeling Length Spectra

library(lme4)
library(lmerTest)
library(emmeans)
library(merTools)



#lm first

# Drop NA
model_df <- drop_na(finfish_length_spectra, landings_raw, bot_temp, b) %>% 
  filter(season %in% c("Spring", "Fall"))

# remove 

# Proportion/residuals of landings not explained by bot temp
landings_temp_lm <- lm(landings_raw ~ bot_temp, data = model_df)
model_df$land_resid <- resid(landings_temp_lm)




# lmer
elmer <- lmerTest::lmer(
  formula = b ~ scale(bot_temp) * region * season + scale(land_resid) * region + (1 | year),
  data = model_df)

summary(elmer)


# ggpredictions
library(ggeffects)
simple_preds <- as.data.frame(
  ggpredict(elmer, ~ bot_temp + season + region) )

# Plot over observed data
simple_preds %>% 
  mutate(
    season = factor(group, levels = c("Spring", "Fall")),
    region = factor(facet, levels = str_replace_all(area_levels_long, " |-", "_"))) %>% 
  ggplot() +
  geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = season), alpha = 0.1) +
  geom_line(
    aes(x, predicted, color = season, linetype = "Predicted"), 
    linewidth = 1) +
  facet_wrap(~region) +
  geom_point(
    data = model_df, 
    aes(bot_temp, b, color = season),
    alpha = 0.6) +
  geom_line(
    data = model_df,
    aes(bot_temp, b, color = season, linetype = "Observed"),
    alpha = 0.3, linewidth = 0.5) +
  scale_linetype_manual(values = c(3, 1)) +
  labs(y = "b")



# not pairwise
slope_phoc <- emtrends(
  object = elmer, 
  specs =  ~ region * season,
  var = "bot_temp",
  adjust = "sidak")



slope_phoc %>% 
  as_tibble() %>% 
  mutate(flag_alpha = ifelse(lower.CL > 0 | upper.CL <0, 1, 0.3)) %>% 
  ggplot(aes(region, bot_temp.trend, ymin = lower.CL, ymax = upper.CL)) +
  geom_hline(yintercept = 0, linetype = 1, color = "gray30", linewidth = 1) +
  geom_pointrange(aes(alpha = I(flag_alpha)), color = gmri_cols("blue"), size = 1) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  facet_wrap(~fct_rev(season), ncol = 2) +
  labs(
    title = "Region * Seasonal Slopes",
    x = NULL,
    y = "yr_num Trend")






# Landings residuals
# not pairwise
slope_phoc <- emtrends(
  object = elmer, 
  specs =  ~ region,
  var = "land_resid",
  adjust = "sidak")

slope_phoc %>% 
  as_tibble() %>% 
  mutate(flag_alpha = ifelse(lower.CL > 0 | upper.CL <0, 1, 0.3)) %>% 
  ggplot(aes(region, land_resid.trend, ymin = lower.CL, ymax = upper.CL)) +
  geom_hline(yintercept = 0, linetype = 1, color = "gray30", linewidth = 1) +
  geom_pointrange(aes(alpha = I(flag_alpha)), color = gmri_cols("blue"), size = 1) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(
    title = "Region Slopes",
    x = NULL,
    y = "yr_num Trend")