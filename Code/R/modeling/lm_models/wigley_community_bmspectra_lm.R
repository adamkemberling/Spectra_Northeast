
#### Modeling Bodymass Spectra  ####
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

# plot theme
theme_set(theme_gmri(rect = element_rect(fill = "white", color = NA)))


# Degree symbol
deg_c <- "\u00b0C"

#### Load Data  ####
wigley_bmspectra_df <- read_csv(here::here("Data/model_ready/wigley_community_bmspectra_mod.csv"))



# vectors for factor levels
area_levels <- c("GoM", "GB", "SNE", "MAB")
area_levels_long <- c("Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight")



# Drop NA's
wtb_model_df <- drop_na(wigley_bmspectra_df, total_weight_lb, bot_temp, b) %>% 
  mutate(yr_num = as.numeric(est_year),
         yr_fac = factor(est_year),
         survey_area = factor(survey_area, levels = area_levels),
         season = factor(season, levels = c("Spring", "Fall")),
         landings = total_weight_lb)




#### Colinearity Avoidance  ####

# Approach 1:
# Proportion/residuals of landings not explained by bot temp

# 1. bodymass spectra model df
landings_temp_lm <- lm(total_weight_lb ~ bot_temp, data = wtb_model_df)
wtb_model_df$land_resid <- resid(landings_temp_lm)
wtb_model_df <- wtb_model_df %>% 
  mutate(
    delta_bt = lag(bot_temp) - bot_temp,
    delta_land = lag(total_weight_lb) - total_weight_lb
  )



####  Quick Plots


# Median Length
ggplot(wtb_model_df, aes(yr_num, b, color = season)) +
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
ggplot(wtb_model_df, aes(yr_num, b, color = season)) +
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
ggplot(wtb_model_df, aes(yr_num, b, color = season)) +
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




####________________####
####  Pass 1: Change Over Time  ####

####____________________####

####__  Pass 1: Change/Time  __####

####  1. Body Mass Spectra Model  ####
wig_wtb_mod <- lm(
  formula = b ~ survey_area * yr_num * season,
  data = wtb_model_df)



# Check important predictors
# Check important predictors
tbl_regression(wig_wtb_mod)  %>% 
  add_global_p(keep = T) %>% 
  bold_p(t = 0.05)  %>%
  bold_labels() %>%
  italicize_levels() 

# Diagnostics
check_model(wig_wtb_mod)
check_outliers(wig_wtb_mod)



##### b. Model Predictions  ####


# Plot the predictions over data

# Full Model
# Plot the predictions over data
wig_wtb_mod_preds <- as.data.frame(
  ggpredict(wig_wtb_mod, ~ yr_num + survey_area + season) )

# Plot over observed data
wig_wtb_mod_preds %>% 
  mutate(
    survey_area = factor(group, levels = area_levels),
    season = factor(facet, levels = c("Spring", "Fall"))) %>% 
  ggplot() +
  geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = season), alpha = 0.1) +
  geom_line(
    aes(x, predicted, color = season), 
    linewidth = 1) +
  geom_point(
    data = wtb_model_df,
    aes(yr_num, b, color = season),
    alpha = 0.4,
    size = 1) +
  facet_wrap(~survey_area, scales = "free") +
  scale_color_gmri() +
  labs(y = "Body Mass Spectra Slope (b)",
       title = "Small community, Body Mass Spectra",
       x = "Year")




##### c. Intercept Post-hoc  ####
# Use emmeans for post-hoc testing for factors



# Regions - significant
(regseas_phoc_lenb <- emmeans(
  wig_wtb_mod, 
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
      y = "EMMean - Body Mass Spectra Slope (b)",
      x = "Group",
      title = "Small Finfish Community",
      subtitle = "Body Mass Spectra Regional Post-Hoc Comparison"))


# Save
ggsave(
  plot = sc_p1_regseasemmeans, 
  filename = here::here("Figs/small_community/sc_wtb_regseason_emmeans.png"))


##### d. Trend Posthoc  ####

# Slope Comparisons from 0

# # Just survey area and year
# (sloperegion_phoc_wtb <- emtrends(
#   object = wig_wtb_mod, 
#   specs =  ~ survey_area,
#   var = "yr_num",
#   adjust = "bonf"))
# 
# 
# # Plotting Slope 
# (sc_p1_yearemtrends <- sloperegion_phoc_wtb %>% 
#     as_tibble() %>% 
#     mutate(flag_alpha = ifelse(lower.CL > 0 | upper.CL <0, 1, 0.4)) %>% 
#     ggplot(aes(survey_area, yr_num.trend, ymin = lower.CL, ymax = upper.CL)) +
#     geom_hline(yintercept = 0, linetype = 2, color = "black") +
#     geom_pointrange(aes(alpha = I(flag_alpha))) +
#     # geom_pointrange(color = "gray30") +
#     scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
#     labs(
#       title = "Small Finfish Community",
#       subtitle = "Body Mass Spectra Annual Trend Coefficients",
#       x = NULL,
#       y = "Trend Coefficient"))
# 
# 
# # Save
# ggsave(
#   plot = sc_p1_yearemtrends, 
#   filename = here::here("Figs/small_community/sc_wtb_year_emtrends.png"))


# Survey area * season * year
(sloperegseason_phoc_wtb <- emtrends(
  object = wig_wtb_mod, 
  specs =  ~ survey_area + season,
  var = "yr_num",
  adjust = "bonf"))


# Plotting Slope 
(sc_p1_yearemtrends <- sloperegseason_phoc_wtb %>% 
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
      subtitle = "Body Mass Spectra Annual Trend Coefficients",
      x = NULL,
      y = "Trend Coefficient"))


# Save
ggsave(
  plot = sc_p1_yearemtrends, 
  filename = here::here("Figs/small_community/sc_wtb_year_emtrends.png"))








###__ Pass 2: Covariates  __####





####  2. Body Mass Spectra Model  ####


## Pass 2: No year term, introduce covariates

# No Interactions version - scaled btemp, log10 landings
# Singular fits for year intercepts b/c temp and landings are annual
wig_wtb_mod2 <- lm(
  formula = b ~ survey_area + log10(landings) + scale(bot_temp) ,
  data = wtb_model_df)

# vif check - seems tolerable
plot(performance::check_collinearity(wig_wtb_mod2)) +
  coord_flip() +
  geom_hline(yintercept = 3, linetype = 3, aes(color = "Zuur Threshold"), linewidth = 1)


# Check important predictors
tbl_regression(wig_wtb_mod2)  %>% 
  add_global_p(keep = T) %>% 
  bold_p(t = 0.05)  %>%
  bold_labels() %>%
  italicize_levels() 


performance::check_model(wig_wtb_mod2)




##### a. Model Predictions  ####

# Not significant in log(weight) model

# Plot marginal effects plots over observed data for:


# Bottom Temperature
bt_preds <- as.data.frame(
  ggpredict(wig_wtb_mod2, ~ bot_temp + survey_area) )


# Plotting over bottom temp differences
(sc_p2_btemp_region_margeffect <- bt_preds %>%
    mutate(
      survey_area = factor(group, levels = area_levels)) %>%
    ggplot() +
    geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = survey_area), alpha = 0.1) +
    geom_line(
      aes(x, predicted, group = survey_area),
      linewidth = 1) +
    geom_point(
      data = wtb_model_df,
      aes(bot_temp, b, color = season),
      alpha = 0.6,
      size = 1) +
    scale_color_gmri() +
    facet_wrap(~survey_area) +
    scale_x_continuous(labels = scales::number_format(suffix = deg_c))+
    labs(
      y = "Body Mass Spectra Slope (b)",
      x = "Bottom Temperature",
      title = "Small Finfish Community",
      subtitle = "Body Mass Spectra and Bottom Temperature Marginal Mean Predictions"
    ))


# Save
ggsave(
  plot = sc_p2_btemp_region_margeffect,
  filename = here::here("Figs/small_community/sc_wtb_btemp_region_margeffects.png"))




# Plot marginal effects plots over observed data for:
# Bottom Temperature
land_preds <- as.data.frame(
  ggpredict(
    model = wig_wtb_mod2,
   terms = list("landings" = 10^c(3:10), 
                 "survey_area" = area_levels))
)



# Plotting over landings differences
(sc_p2_land_region_margeffect <- land_preds %>%
    mutate(
      survey_area = factor(group, levels = area_levels)) %>%
    ggplot() +
    geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = survey_area), alpha = 0.1) +
    geom_line(
      aes(x, predicted, group = survey_area),
      linewidth = 1) +
    geom_point(
      data = wtb_model_df,
      aes(landings, b, color = season),
      alpha = 0.6,
      size = 1) +
    scale_color_gmri() +
    facet_wrap(~survey_area) +
    scale_x_log10(labels = scales::label_log(10))+
    labs(
      y = "Body Mass Spectra Slope (b)",
      x = "Total Landings (lb.)",
      title = "Small Finfish Community",
      subtitle = "Body Mass Spectra and Lendings Marginal Mean Predictions"
    ))


# Save
ggsave(
  plot = sc_p2_land_region_margeffect,
  filename = here::here("Figs/small_community/sc_wtb_l10landings_region_margeffects.png"))





##### a.  Intercept Post-Hoc  ####

# Summary 
summary(wig_wtb_mod2)
# area
# season
# area:season


# Region and Season Interaction - Significant

# Response scale
(regseas_phoc_wtb2 <- emmeans(
  object = wig_wtb_mod2,
  specs = list(pairwise ~ survey_area),
  adjust = "tukey", 
  type = "response"))


# Plot of the region & seasonal differences
(sc_p2_regionseason_emmeans <- regseas_phoc_wtb2$`emmeans of survey_area` %>%
    as_tibble() %>%
    ggplot(aes(survey_area, emmean, ymin = lower.CL, ymax = upper.CL)) +
    geom_pointrange(position = position_dodge(width = 0.25), alpha = 0.8) +
    scale_color_gmri() +
    labs(
      y = "Body Mass Spectra Slope (b)",
      x = "Region",
      title = "Small Finfish Community",
      color = "Season",
      subtitle = "Body Mass Spectra Region x Season Fixed Effects - Post-Hoc"))



##### b.  Trend Marginal Effects  ####


emtrends(
  object = wig_wtb_mod2,
  specs =  ~ survey_area,
  var = "bot_temp",
  adjust = "bonf")


# Just temp
# Plot marginal effects plots over observed data for:
# Bottom Temperature
(sc_p2_btemp_margeffect <- as.data.frame(
  ggpredict(wig_wtb_mod2, ~ bot_temp) ) %>%
    ggplot(aes(x, predicted, ymin = conf.low, ymax = conf.high)) +
    geom_ribbon(alpha = 0.1) +
    geom_line() +
    labs(
      y = "Body Mass Spectra Slope (b)",
      x = "Bottom Temperature",
      title = "Small Finfish Community",
      subtitle = "Body Mass Spectra and Bottom Temperature Marginal Mean Effect"
    ))


# Save
ggsave(
  plot = sc_p2_btemp_margeffect,
  filename = here::here("Figs/small_community/sc_wtb_btemp_margeffects.png"))




# Significant Trend in landings
emtrends(
  object = wig_wtb_mod2, 
  specs =  ~ survey_area,
  var = "landings",
  adjust = "bonf")


# Plot marginal effects plots over observed data for:
# Landings Effect
(sc_p2_land_margeffect <- as.data.frame(
  ggpredict(
    model = wig_wtb_mod2, 
    terms = list("landings" = 10^c(1:10)))) %>% 
    ggplot(aes(x, predicted, ymin = conf.low, ymax = conf.high)) +
    geom_ribbon(alpha = 0.1) +
    geom_line() +
    scale_x_log10(labels = scales::label_log(10),
                  limits = c(10^1, 10^10))+
    labs(
      y = "Body Mass Spectra Slope (b)", 
      x = "Total Landings (lb.)",
      title = "Small Finfish Community",
      subtitle = "Body Mass Spectra and Landings Marginal Mean Effect"
    ))


# Save
ggsave(
  plot = sc_p2_land_margeffect, 
  filename = here::here("Figs/small_community/sc_wtb_landings_margeffects.png"))





####  Export Results  ####
model_pack <- list(
  "pass_1" = wig_wtb_mod,
  "pass_2" = wig_wtb_mod2,
  "model_data" = wtb_model_df)

# Save them
saveRDS(object = model_pack, file = here::here("Data/models_and_results/sc_bmspectra.RDS"))




