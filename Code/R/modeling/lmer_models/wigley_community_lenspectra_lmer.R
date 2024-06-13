
#### Modeling Length Spectra  ####
# Reduced finfish community

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

theme_set(theme_gmri())

# Degree symbol
deg_c <- "\u00b0C"

#### Load Data  ####
wigley_medlen_df <- read_csv(here::here("Data/model_ready/wigley_community_medsize_mod.csv"))
wigley_lenspectra_df <- read_csv(here::here("Data/model_ready/wigley_community_lenspectra_mod.csv"))
wigley_bmspectra_df <- read_csv(here::here("Data/model_ready/wigley_community_bmspectra_mod.csv"))



# vectors for factor levels
area_levels <- c("GoM", "GB", "SNE", "MAB")
area_levels_long <- c("Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight")



# Drop NA's
lenb_model_df <- drop_na(wigley_lenspectra_df, total_weight_lb, bot_temp, b) %>% 
  mutate(yr_num = as.numeric(est_year),
         yr_fac = factor(est_year),
         survey_area = factor(survey_area, levels = area_levels),
         season = factor(season, levels = c("Spring", "Fall")),
         landings = total_weight_lb)




#### Colinearity Avoidance  ####

# Approach 1:
# Proportion/residuals of landings not explained by bot temp

# Approach 2:
# Differencing of tthe two predictors to remove trend


# 1. length spectra model df
landings_temp_lm <- lm(total_weight_lb ~ bot_temp, data = lenb_model_df)
lenb_model_df$land_resid <- resid(landings_temp_lm)
lenb_model_df <- lenb_model_df %>% 
  mutate(
    delta_bt = lag(bot_temp) - bot_temp,
    delta_land = lag(total_weight_lb) - total_weight_lb
  )



####  Quick Plots


# Median Length
ggplot(lenb_model_df, aes(yr_num, b, color = season)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_ma(aes(linetype = "5-Year Moving Average"),n = 5, ma_fun = SMA) +
  geom_smooth(method = "lm", linewidth = 1, se = T, 
              aes(linetype = "Regression Fit")) +
  facet_wrap(~survey_area, scales = "free") +
  scale_color_gmri() +
  labs(title = "Median Length",
       subtitle = "Finfish Community",
       y = "Length (cm)",
       x = "Year",
       color = "Season")


# Length Spectra Slope
ggplot(lenb_model_df, aes(yr_num, b, color = season)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_ma(aes(linetype = "5-Year Moving Average"),n = 5, ma_fun = SMA) +
  geom_smooth(method = "lm", linewidth = 1, se = F, 
              aes(linetype = "Regression Fit")) +
  facet_wrap(~survey_area, scales = "free") +
  scale_color_gmri() +
  labs(title = "Length Spectra Slope - MLE Bins Method",
       subtitle = "Finfish Community\nEnforced xmin = 1, xmax = max(length_cm + 1)",
       y = "b",
       x = "Year",
       color = "Season")

# Bodymass Spectra Slope
ggplot(lenb_model_df, aes(yr_num, b, color = season)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_ma(aes(linetype = "5-Year Moving Average"),n = 5, ma_fun = SMA) +
  geom_smooth(method = "lm", linewidth = 1, se = F, 
              aes(linetype = "Regression Fit")) +
  facet_wrap(~survey_area) +
  scale_color_gmri() +
  labs(title = "Bodymass Spectra Slope - MLE Bins Method",
       subtitle = "Finfish Community\nEnforced xmin = 1, xmax = max(length_cm + 1)",
       y = "b",
       x = "Year",
       color = "Season")

# Covariates change over time
lenb_model_df %>% 
  select(est_year, survey_area, total_weight_lb, bot_temp) %>% 
  mutate(bot_temp = scale(bot_temp)[,1],
         total_weight_lb = scale(total_weight_lb)[,1]) %>% 
  pivot_longer(-c(est_year, survey_area), names_to = "Covariate", values_to = "value") %>% 
  ggplot(aes(est_year, value, color = Covariate)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_ma(aes(linetype = "5-Year Moving Average"),n = 5, ma_fun = SMA,alpha = 0.6) +
  geom_smooth(method = "lm", linewidth = 1, se = F, 
              aes(linetype = "Regression Fit")) +
  facet_wrap(~survey_area) +
  scale_color_gmri() +
  labs(title = "Predictor Trends",
       subtitle = "Exploring colinearity with time",
       y = "b",
       x = "Year",
       color = "Covariate")




####____________________####

####__  Pass 1: Change/Time  __####

####  1. Body Length Spectra Model  ####
wig_lenb_mod <- lmerTest::lmer(
  formula = b ~ survey_area * yr_num * season + (1 | yr_fac),
  data = lenb_model_df)



# Check important predictors
summary(wig_lenb_mod)
# survey area 
# year
# area * year 
# area * season
# area * year * season


# Diagnostics
check_model(wig_lenb_mod)
check_outliers(wig_lenb_mod)


# Calculate intra-class correlation
icc(wig_lenb_mod)
# 6.1% of total variance is attributable to between year variance






##### b. Model Predictions  ####


# Plot the predictions over data

# Full Model
# Plot the predictions over data
wig_lenb_mod_preds <- as.data.frame(
  ggpredict(wig_lenb_mod, ~ yr_num + survey_area + season) )

# Plot over observed data
wig_lenb_mod_preds %>% 
  mutate(
    survey_area = factor(group, levels = area_levels),
    season = factor(facet, levels = c("Spring", "Fall"))) %>% 
  ggplot() +
  geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = season), alpha = 0.1) +
  geom_line(
    aes(x, predicted, color = season), 
    linewidth = 1) +
  geom_point(
    data = lenb_model_df,
    aes(yr_num, b, color = season),
    alpha = 0.4,
    size = 1) +
  facet_wrap(~survey_area, scales = "free") +
  scale_color_gmri() +
  labs(y = "Body Length Spectra Slope (b)",
       title = "Small community, Body Length Spectra",
       x = "Year")




##### c. Intercept Post-hoc  ####
# Use emmeans for post-hoc testing for factors



# # Regions - significant
# (region_phoc_lenb <- emmeans(
#   wig_lenb_mod, 
#   list(pairwise ~ survey_area), 
#   adjust = "tukey",
#   type = "response"))
# 
# 
# 
# # Custom Plot
# (sc_p1_regionemmeans <- region_phoc_lenb$`emmeans of survey_area` %>% 
#     as_tibble() %>%
#     ggplot(aes(survey_area, emmean, ymin = lower.CL, ymax = upper.CL)) +
#     geom_pointrange(position = position_dodge(width = 0.25), alpha = 0.8) +
#     scale_color_gmri() +
#     labs(
#       y = "EMMean - Body Length Spectra Slope (b)",
#       x = "Region",
#       title = "Small Finfish Community",
#       subtitle = "Body Length Spectra Regional Post-Hoc Comparison"))
# 
# 
# # Save
# ggsave(
#   plot = sc_p1_regionemmeans, 
#   filename = here::here("Figs/small_community/sc_lenb_region_emmeans.png"))


# Regions - significant
(regseas_phoc_lenb <- emmeans(
  wig_lenb_mod, 
  list(pairwise ~ survey_area + season), 
  adjust = "tukey",
  type = "response"))



# Custom Plot
(sc_p1_regseasemmeans <- regseas_phoc_lenb$`emmeans of survey_area, season` %>% 
    as_tibble() %>%
    ggplot(aes(survey_area, emmean, ymin = lower.CL, ymax = upper.CL, color = season)) +
    geom_pointrange(position = position_dodge(width = 0.25), alpha = 0.8) +
    scale_color_gmri() +
    labs(
      y = "EMMean - Body Length Spectra Slope (b)",
      x = "Group",
      title = "Small Finfish Community",
      subtitle = "Body Length Spectra Regional Post-Hoc Comparison"))


# Save
ggsave(
  plot = sc_p1_regseasemmeans, 
  filename = here::here("Figs/small_community/sc_lenb_regionseason_emmeans.png"))




##### d. Trend Posthoc  ####


# # Just survey area and year
# (sloperegion_phoc_lenb <- emtrends(
#   object = wig_lenb_mod, 
#   specs =  ~ survey_area,
#   var = "yr_num",
#   adjust = "sidak"))
# 
# 
# # Plotting Slope 
# (sc_p1_yearemtrends <- sloperegion_phoc_lenb %>% 
#     as_tibble() %>% 
#     mutate(flag_alpha = ifelse(lower.CL > 0 | upper.CL <0, 1, 0.4)) %>% 
#     ggplot(aes(survey_area, yr_num.trend, ymin = lower.CL, ymax = upper.CL)) +
#     geom_hline(yintercept = 0, linetype = 2, color = "black") +
#     geom_pointrange(aes(alpha = I(flag_alpha))) +
#     # geom_pointrange(color = "gray30") +
#     scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
#     labs(
#       title = "Small Finfish Community",
#       subtitle = "Body Length Spectra Annual Trend Coefficients",
#       x = NULL,
#       y = "Trend Coefficient"))
# 
# 
# # Save
# ggsave(
#   plot = sc_p1_yearemtrends, 
#   filename = here::here("Figs/small_community/sc_lenb_year_emtrends.png"))


# Survey area * season * year
(sloperegseason_phoc_lenb <- emtrends(
  object = wig_lenb_mod, 
  specs =  ~ survey_area + season,
  var = "yr_num",
  adjust = "sidak"))


# Plotting Slope 
(sc_p1_yearemtrends <- sloperegseason_phoc_lenb %>% 
    as_tibble() %>% 
    mutate(flag_alpha = ifelse(lower.CL > 0 | upper.CL <0, 1, 0.4)) %>% 
    ggplot(
      aes(survey_area, yr_num.trend, ymin = lower.CL, ymax = upper.CL, color = season)) +
    geom_hline(yintercept = 0, linetype = 2, color = "black") +
    geom_pointrange(
      aes(alpha = I(flag_alpha)),
      position = position_dodge(width = 0.25)) +
    scale_color_gmri() +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    labs(
      title = "Small Finfish Community",
      subtitle = "Body Length Spectra Annual Trend Coefficients",
      x = NULL,
      y = "Trend Coefficient"))


# Save
ggsave(
  plot = sc_p1_yearemtrends, 
  filename = here::here("Figs/small_community/sc_lenb_year_emtrends.png"))









###__ Pass 2: Covariates  __####





####  2. Body Length Spectra Model  ####


## Pass 2: No year term, introduce covariates

# No Interactions version - scaled btemp, log10 landings
# Singular fits for year intercepts b/c temp and landings are annual
wig_lenb_mod2 <- lmerTest::lmer(
  formula = b ~ survey_area * season + log10(landings) + scale(bot_temp) + (1 | yr_fac),
  data = lenb_model_df)

# vif check - seems tolerable
plot(performance::check_collinearity(wig_lenb_mod2)) +
  coord_flip() +
  geom_hline(yintercept = 3, linetype = 3, aes(color = "Zuur Threshold"), linewidth = 1)


# Summary 
summary(wig_lenb_mod2)
# log10(landings)
# scale(bot_temp)
# survey area
# season
# season:area



check_model(wig_lenb_mod2)
check_outliers(wig_lenb_mod2)



# Calculate intra-class correlation
# ratio of the random intercept variance (between year variance)
# to the total variance
icc(wig_lenb_mod2)
# 3%




##### b. Model Predictions  ####


# Plot marginal effects plots over observed data for:
# Bottom Temperature
bt_preds <- as.data.frame(
  ggpredict(wig_lenb_mod2, ~ bot_temp + survey_area + season) )


# Plotting over bottom temp differences
(sc_p2_btemp_region_margeffect <- bt_preds %>%
    mutate(
      season = factor(facet, levels = c("Spring", "Fall")),
      survey_area = factor(group, levels = area_levels)) %>%
    ggplot() +
    geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = season), alpha = 0.1) +
    geom_line(
      aes(x, predicted, color = season),
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
      y = "Body Length Spectra Slope (b)",
      x = "Bottom Temperature",
      title = "Small Finfish Community",
      subtitle = "Body Length Spectra and Bottom Temperature Marginal Mean Predictions"
    ))


# # Save
# ggsave(
#   plot = sc_p2_btemp_region_margeffect, 
#   filename = here::here("Figs/small_community/sc_lenb_btemp_region_margeffects.png"))




# Plot marginal effects plots over observed data for:
# Bottom Temperature
land_preds <- as.data.frame(
  ggpredict(
    model = wig_lenb_mod2,
    #terms = ~ landings + survey_area + season)
    terms = list("landings" = 10^c(3:10), 
                 "survey_area" = area_levels,
                 season = c("Spring", "Fall")))
)



# Plotting over landings differences
(sc_p2_land_region_margeffect <- land_preds %>%
    mutate(
      season = factor(facet, levels = c("Spring", "Fall")),
      survey_area = factor(group, levels = area_levels)) %>%
    ggplot() +
    geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = season), alpha = 0.1) +
    geom_line(
      aes(x, predicted, color = season),
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
      y = "Body Length Spectra Slope (b)",
      x = "Total Landings (lb.)",
      title = "Small Finfish Community",
      subtitle = "Body Length Spectra and Lendings Marginal Mean Predictions"
    ))


# # Save
# ggsave(
#   plot = sc_p2_btemp_region_margeffect,
#   filename = here::here("Figs/small_community/sc_lenb_l10landings_region_margeffects.png"))





##### a.  Intercept Post-Hoc  ####

# Summary 
summary(wig_lenb_mod2)
# area
# season
# area:season


# Region and Season Interaction - Significant

# emmean coefs
regseas_phoc_lenb2 <- emmeans(
  object = wig_lenb_mod2,
  specs = list(pairwise ~ survey_area * season),
  adjust = "tukey")

regseas_phoc_lenb2


# Response scale
regseas_phoc_lenb2 <- emmeans(
  object = wig_lenb_mod2,
  specs = list(pairwise ~ survey_area * season),
  adjust = "tukey", 
  type = "response")


# Plot of the region & seasonal differences
(sc_p2_regionseason_emmeans <- regseas_phoc_lenb2$`emmeans of survey_area, season` %>%
    as_tibble() %>%
    ggplot(aes(survey_area, emmean, ymin = lower.CL, ymax = upper.CL)) +
    geom_pointrange(aes(color = season), position = position_dodge(width = 0.25), alpha = 0.8) +
    scale_color_gmri() +
    labs(
      y = "Body Length Spectra Slope (b)",
      x = "Region",
      title = "Small Finfish Community",
      color = "Season",
      subtitle = "Body Length Spectra Region x Season Fixed Effects - Post-Hoc"))


# Save
ggsave(
  plot = sc_p2_regionseason_emmeans, 
  filename = here::here("Figs/small_community/sc_lenb_regionseason_emmeans.png"))




##### b.  Trend Marginal Effects  ####

# Significant Trend in Bot temp
emtrends(
  object = wig_lenb_mod2, 
  specs =  ~ survey_area,
  var = "bot_temp",
  adjust = "sidak")


# Just temp
# Plot marginal effects plots over observed data for:
# Bottom Temperature
(sc_p2_btemp_margeffect <- as.data.frame(
  ggpredict(wig_lenb_mod2, ~ bot_temp) ) %>%
    ggplot(aes(x, predicted, ymin = conf.low, ymax = conf.high)) +
    geom_ribbon(alpha = 0.1) +
    geom_line() +
    labs(
      y = "Body Length Spectra Slope (b)",
      x = "Bottom Temperature",
      title = "Small Finfish Community",
      subtitle = "Body Length Spectra and Bottom Temperature Marginal Mean Effect"
    ))


# Save
ggsave(
  plot = sc_p2_btemp_margeffect,
  filename = here::here("Figs/small_community/sc_lenb_btemp_margeffects.png"))




# Significant Trend in landings
emtrends(
  object = wig_lenb_mod2, 
  specs =  ~ survey_area,
  var = "landings",
  adjust = "sidak")


# Plot marginal effects plots over observed data for:
# Landings Effect
(sc_p2_land_margeffect <- as.data.frame(
  ggpredict(
    model = wig_lenb_mod2, 
    terms = list("landings" = 10^c(1:10)))) %>% 
    ggplot(aes(x, predicted, ymin = conf.low, ymax = conf.high)) +
    geom_ribbon(alpha = 0.1) +
    geom_line() +
    scale_x_log10(labels = scales::label_log(10),
                  limits = c(10^1, 10^10))+
    labs(
      y = "Body Length Spectra Slope (b)", 
      x = "Total Landings (lb.)",
      title = "Small Finfish Community",
      subtitle = "Body Length Spectra and Landings Marginal Mean Effect"
    ))


# Save
ggsave(
  plot = sc_p2_land_margeffect, 
  filename = here::here("Figs/small_community/sc_lenb_landings_margeffects.png"))

