#### Modeling Large Community Length Spectra  ####

library(lme4)
library(lmerTest)
library(emmeans)
library(merTools)
library(tidyquant)
library(tidyverse)
library(gmRi)
library(ggeffects)
library(scales)
library(performance)
library(gtsummary)

# theme set
theme_set(theme_gmri(rect = element_rect(fill = "white", color = NA)))

#### Load Data  ####
ffish_lenspectra_df <- read_csv(here::here("Data/model_ready/large_community_lenspectra_mod.csv"))

# vectors for factor levels
area_levels <- c("GoM", "GB", "SNE", "MAB")
area_levels_long <- c("Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight")

# Degree symbol
deg_c <- "\u00b0C"


# Drop NA's
lenb_model_df <- drop_na(ffish_lenspectra_df, total_weight_lb, bot_temp, b) %>% 
  mutate(yr_num = as.numeric(est_year),
         yr_fac = factor(est_year),
         survey_area = factor(survey_area, levels = area_levels),
         season = factor(season, levels = c("Spring", "Fall")),
         landings = total_weight_lb)




####  Overview Plots  ####


# Spectra Slope
ggplot(lenb_model_df, aes(yr_num, b, color = season)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_ma(aes(linetype = "5-Year Moving Average"),n = 5, ma_fun = SMA) +
  geom_smooth(method = "lm", linewidth = 1, se = F, 
              aes(linetype = "Regression Fit")) +
  facet_wrap(~survey_area) +
  scale_color_gmri() +
  labs(title = "Length Spectra Slope - MLE Bins Method",
       subtitle = "Finfish Community\nEnforced xmin = 1, xmax = max(length_cm + 1)",
       y = "b",
       x = "Year",
       color = "Season")


# Bottom Temperatures
# Plot the annual trends
ggplot(lenb_model_df, aes(yr_num, bot_temp, color = season)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_ma(aes(linetype = "5-Year Moving Average"),n = 5, ma_fun = SMA) +
  geom_smooth(method = "lm", linewidth = 1, se = F, 
              aes(linetype = "Regression Fit")) +
  facet_wrap(~survey_area) +
  scale_color_gmri() +
  scale_y_continuous(labels = label_number(suffix = deg_c)) +
  labs(title = "Regional Bottom Temperatures",
       y = "Temperature",
       x = "Year",
       color = "Area")




####____________________####
#### Length Spectra Models  ####
####__  Pass 1: Change/Time  __####

#### 1. Length Spectra Model  ####

ffish_b_mod <- lm(
  formula = b ~ survey_area * yr_num * season,
  data = lenb_model_df)


# Check important predictors
tbl_regression(ffish_b_mod)  %>% 
  add_global_p(keep = T) %>% 
  bold_p(t = 0.05)  %>%
  bold_labels() %>%
  italicize_levels() 






##### b. Model Predictions  ####

# Plot the predictions over data
ffish_b_mod_preds <- as.data.frame(
  ggpredict(ffish_b_mod, ~ yr_num + survey_area) )


# Plot over observed data
ffish_b_mod_preds %>% 
  mutate(survey_area = factor(group, levels = area_levels)) %>% 
  ggplot() +
  geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = survey_area), alpha = 0.1) +
  geom_line(
    aes(x, predicted, color = survey_area), 
    linewidth = 1) +
  geom_point(
    data = lenb_model_df,
    aes(yr_num, b, color = survey_area),
    alpha = 0.4,
    size = 1) +
  facet_wrap(~survey_area, scales = "free") +
  scale_color_gmri() +
  labs(y = "b",
       title = "large community, length spectra",
       x = "Year")




##### c. Intercept Post-hoc  ####
# Use emmeans for post-hoc testing for factors

# Regions - significant
(region_phoc_b <- emmeans(ffish_b_mod, list(pairwise ~ survey_area), adjust = "tukey"))


# Custom Plot
(lc_p1_regionemmeans_b <- region_phoc_b$`emmeans of survey_area` %>% 
    as_tibble() %>%
    ggplot(aes(survey_area, emmean, ymin = lower.CL, ymax = upper.CL)) +
    geom_pointrange(position = position_dodge(width = 0.25), alpha = 0.8) +
    scale_color_gmri() +
    labs(
      y = "EMMean - Length Spectra",
      x = "Region",
      title = "Large Finfish Community",
      subtitle = "Length Spectra Regional Post-Hoc Comparison"))


# Save
ggsave(
  plot = lc_p1_regionemmeans_b, 
  filename = here::here("Figs/large_community/lc_lenb_region_emmeans.png"))




##### d. Trend Posthoc  ####

# Slope Comparisons
# Just survey area and year
sloperegion_phoc_lenb <- emtrends(
  object = ffish_b_mod, 
  specs =  ~ survey_area,
  var = "yr_num",
  adjust = "bonf")

# Values
sloperegion_phoc_lenb


# Plotting Slope 
(lc_p1_yearemtrends_b <- sloperegion_phoc_lenb %>% 
    as_tibble() %>% 
    mutate(flag_alpha = ifelse(lower.CL > 0 | upper.CL <0, 1, 0.4)) %>% 
    ggplot(aes(survey_area, yr_num.trend, ymin = lower.CL, ymax = upper.CL)) +
    geom_hline(yintercept = 0, linetype = 2, color = "black") +
    geom_pointrange(aes(alpha = I(flag_alpha))) +
    #geom_pointrange(color = "gray30") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    labs(
      title = "Large Finfish Community",
      subtitle = "Length Spectra Annual Trend Coefficients",
      x = NULL,
      y = "Trend Coefficient"))


# Save
ggsave(
  plot = lc_p1_yearemtrends_b, 
  filename = here::here("Figs/large_community/lc_lenb_year_emtrends.png"))




###__ Pass 2: Covariates  __####



#### 2. Length Spectra Model  ####
ffish_b_mod2 <- lm(
  formula = b ~ survey_area + log10(landings) + scale(bot_temp),
  data = lenb_model_df)

# vif check - seems still bad
plot(performance::check_collinearity(ffish_b_mod2)) +
  coord_flip() +
  geom_hline(yintercept = 3, linetype = 3, aes(color = "Zuur Threshold"), linewidth = 1)


# Check important predictors
tbl_regression(ffish_b_mod2)  %>% 
  add_global_p(keep = T) %>% 
  bold_p(t = 0.05)  %>%
  bold_labels() %>%
  italicize_levels() 





##### b. Model Predictions  ####


# Plot marginal effects plots over observed data for:
# Bottom Temperature
bt_preds_b <- as.data.frame(
  ggpredict(ffish_b_mod2, ~ bot_temp + survey_area) )


# Plotting over bottom temp differences
(lc_p2_btemp_region_margeffect_b <- bt_preds_b %>% 
    mutate(
      survey_area = factor(group, levels = area_levels)) %>% 
    ggplot() +
    geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = survey_area), alpha = 0.1) +
    geom_line(
      aes(x, predicted, group = survey_area), 
      linewidth = 1) +
    geom_point(
      data = lenb_model_df, 
      aes(bot_temp, b, color = season),
      alpha = 0.6,
      size = 1) +
    scale_color_gmri() +
    facet_wrap(~survey_area) +
    scale_x_continuous(labels = scales::number_format(suffix = deg_c))+
    labs(
      y = "b", 
      x = "Bottom Temperature",
      title = "Large Finfish Community",
      subtitle = "Length Spectra and Bottom Temperature Marginal Mean Predictions"
    ))


# Save
ggsave(
  plot = lc_p2_btemp_region_margeffect_b, 
  filename = here::here("Figs/large_community/lc_lenb_btemp_regseason_margeffects.png"))





# Plot marginal effects plots over observed data for:
# Bottom Temperature
land_preds_b <- as.data.frame(
  ggpredict(
    model = ffish_b_mod2,
    terms = list("landings" = 10^c(1:10), 
                 "survey_area" = area_levels))
)



# Plotting over landings differences
(lc_p2_land_region_margeffect_b <- land_preds_b %>%
    mutate(
      survey_area = factor(group, levels = area_levels)) %>%
    ggplot() +
    geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = survey_area), alpha = 0.1) +
    geom_line(
      aes(x, predicted, group = survey_area),
      linewidth = 1) +
    geom_point(
      data = lenb_model_df,
      aes(landings, b, color = season),
      alpha = 0.6,
      size = 1) +
    scale_color_gmri() +
    facet_wrap(~survey_area) +
    scale_x_log10(labels = scales::label_log(10))+
    labs(
      y = "Length Spectra Slope (b)",
      x = "Total Landings (lb.)",
      title = "Large Finfish Community",
      subtitle = "Median Length and Lendings Marginal Mean Predictions"
    ))


# Save
ggsave(
  plot = lc_p2_btemp_region_margeffect_b,
  filename = here::here("Figs/large_community/lc_lenb_l10landings_regseason_margeffects.png"))









##### c.  Intercept Post-Hoc  ####

# Region and Season Interaction - Significant
(regseas_phoc_b2 <- emmeans(
  object = ffish_b_mod2,
  specs = list(pairwise ~ survey_area),
  adjust = "tukey"))

# Plot of the region & seasonal differences
regseas_phoc_b2$`emmeans of survey_area` %>%
  as_tibble() %>%
  ggplot(aes(survey_area, emmean, ymin = lower.CL, ymax = upper.CL)) +
  geom_pointrange(position = position_dodge(width = 0.25), alpha = 0.8) +
  scale_color_gmri() +
  labs(
    y = "EMMean - length spectra (b)",
    x = "Region",
    title = "Factor Fixed Effects - Post-Hoc")




##### d.  Trend Post-Hoc  ####

# Significant Trend in Bot temp
emtrends(
  object = ffish_b_mod2, 
  specs =  ~ survey_area,
  var = "bot_temp",
  adjust = "bonf")



# Bottom Temperature
(lc_p2_btemp_margeffect_b <- as.data.frame(
  ggpredict(ffish_b_mod2, ~ bot_temp) ) %>% 
    ggplot(aes(x, predicted, ymin = conf.low, ymax = conf.high)) +
    geom_ribbon(alpha = 0.1) +
    geom_line() +
    scale_x_continuous(labels = label_number(suffix = deg_c)) +
    labs(
      y = "Length Spectra (b)", 
      x = "Bottom Temperature",
      title = "Large Finfish Community",
      subtitle = "Median Length and Bottom Temperature Marginal Mean Effect"
    ))



# Save
ggsave(
  plot = lc_p2_btemp_margeffect_b, 
  filename = here::here("Figs/large_community/lc_lenb_btemp_region_margeffects.png"))



# Significant Trend in landings
emtrends(
  object = ffish_b_mod2, 
  specs =  ~ survey_area,
  var = "landings",
  adjust = "bonf")



# Plot marginal effects plots over observed data for:
# Landings Effect
(lc_p2_land_margeffect_b <- as.data.frame(
  ggpredict(
    model = ffish_b_mod2, 
    terms = list("landings" = 10^c(1:10)))) %>% 
    ggplot(aes(x, predicted, ymin = conf.low, ymax = conf.high)) +
    geom_ribbon(alpha = 0.1) +
    geom_line() +
    scale_x_log10(labels = scales::label_log(10),
                  limits = c(10^6, 10^9))+
    labs(
      y = "Length Spectra (b)", 
      x = "Total Landings (lb.)",
      title = "Large Finfish Community",
      subtitle = "Median Length and Landings Marginal Mean Effect"
    ))


# Save
ggsave(
  plot = lc_p2_land_margeffect_b, 
  filename = here::here("Figs/large_community/lc_lenb_landings_margeffects.png"))



####  Export Results  ####
model_pack <- list(
  "pass_1" = ffish_b_mod,
  "pass_2" = ffish_b_mod2,
  "model_data" = lenb_model_df)

# Save them
saveRDS(object = model_pack, file = here::here("Data/models_and_results/lc_lenspectra.RDS"))
read_rds(here::here("Data/models_and_results/lc_lenspectra.RDS"))

