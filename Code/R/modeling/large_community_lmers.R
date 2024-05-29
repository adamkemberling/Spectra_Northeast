
#### Modeling Length Spectra  ####

library(lme4)
library(lmerTest)
library(emmeans)
library(merTools)



#### Load Data  ####
ffish_medlen_df <- read_csv(here::here("Data/model_ready/large_community_medlength_mod.csv"))
ffish_lenspectra_df <- read_csv(here::here("Data/model_ready/large_community_lenspectra_mod.csv"))

#lm first

# Drop NA's
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