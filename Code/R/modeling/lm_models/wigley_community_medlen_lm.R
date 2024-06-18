
#### Modeling Length Spectra  ####
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


#### Load Data  ####
wigley_medlen_df <- read_csv(here::here("Data/model_ready/wigley_community_medsize_mod.csv"))



# vectors for factor levels
area_levels <- c("GoM", "GB", "SNE", "MAB")
area_levels_long <- c("Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight")

# Degree symbol
deg_c <- "\u00b0C"



# Drop NA's
len_model_df <- drop_na(wigley_medlen_df, total_weight_lb, bot_temp, med_len_cm) %>% 
  mutate(yr_num = as.numeric(est_year),
         yr_fac = factor(est_year),
         survey_area = factor(survey_area, levels = area_levels),
         season = factor(season, levels = c("Spring", "Fall")))





#### Colinearity Avoidance  ####



# 1. length model df
landings_temp_lm <- lm(total_weight_lb ~ bot_temp, data = len_model_df)
len_model_df$land_resid <- resid(landings_temp_lm)
len_model_df <- len_model_df %>% 
  mutate(yr_num = as.numeric(est_year),
         yr_fac = factor(est_year),
         survey_area = factor(survey_area, levels = area_levels),
         season = factor(season, levels = c("Spring", "Fall")),
         landings = total_weight_lb)



####  Quick Plots


# Median Length
ggplot(len_model_df, aes(yr_num, med_len_cm, color = season)) +
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



# Covariates over time time
len_model_df %>% 
  ggplot(aes(est_year, bot_temp, color = season)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_ma(aes(linetype = "5-Year Moving Average", group = season),
          n = 5, ma_fun = SMA,alpha = 0.6) +
  geom_smooth(method = "lm", linewidth = 1, se = F, 
              aes(linetype = "Regression Fit", group = season)) +
  facet_wrap(~survey_area) +
  scale_color_gmri() +
  labs(title = "Predictor Trends",
       subtitle = "Exploring colinearity with time",
       x = "Year",
       color = "Season")


# Landings
len_model_df %>% 
  ggplot(aes(est_year, landings, color = survey_area)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_ma(aes(linetype = "5-Year Moving Average"),
          n = 5, ma_fun = SMA,alpha = 0.6) +
  geom_smooth(method = "lm", linewidth = 1, se = F, 
              aes(linetype = "Regression Fit")) +
  facet_wrap(~survey_area) +
  scale_color_gmri() +
  labs(
    title = "Predictor Trends",
    subtitle = "Exploring colinearity with time",
    x = "Year")




####____________________####

####__  Pass 1: Change/Time  __####

####  1. Median Length Model  ####
wig_len_mod <- lm(
  formula = med_len_cm ~ survey_area * yr_num * season,
  data = len_model_df)

# Summary
summary(wig_len_mod)


# Check important predictors
# This is how gtsummary is doing effect testing:
# car::Anova(wig_wt_mod)
tbl_regression(wig_len_mod)  %>% 
  add_global_p(keep = T) %>% 
  bold_p(t = 0.05)  %>%
  bold_labels() %>%
  italicize_levels() 


# Diagnostics
# check_model(wig_len_mod)
check_outliers(wig_len_mod)
check_collinearity(wig_len_mod)
check_autocorrelation(wig_len_mod)







##### b. Model Predictions  ####


# Plot the predictions over data


# Full Model
# Plot the predictions over data
wig_len_mod_preds <- as.data.frame(
  ggpredict(wig_len_mod, ~ yr_num + survey_area + season) )

# Plot over observed data
wig_len_mod_preds %>% 
  mutate(
    survey_area = factor(group, levels = area_levels),
    season = factor(facet, levels = c("Spring", "Fall"))) %>% 
  ggplot() +
  geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = season), alpha = 0.1) +
  geom_line(
    aes(x, predicted, color = season), 
    linewidth = 1) +
  geom_point(
    data = len_model_df,
    aes(yr_num, med_len_cm, color = season),
    alpha = 0.4,
    size = 1) +
  facet_wrap(~survey_area, scales = "free") +
  scale_color_gmri() +
  labs(y = "Median Length (cm)",
       title = "Small community, median weight",
       x = "Year")




##### c. Intercept Post-hoc  ####
# Use emmeans for post-hoc testing for factors

# Regions - significant
(region_phoc_len <- emmeans(
  wig_len_mod, 
  list(pairwise ~ survey_area), 
  adjust = "tukey",
  type = "response"))



# Custom Plot
(sc_p1_regionemmeans <- region_phoc_len$`emmeans of survey_area` %>% 
    as_tibble() %>%
    ggplot(aes(survey_area, emmean, ymin = lower.CL, ymax = upper.CL)) +
    geom_pointrange(position = position_dodge(width = 0.25), alpha = 0.8) +
    scale_color_gmri() +
    labs(
      y = "EMMean - Median Length (cm)",
      x = "Region",
      title = "Small Finfish Community",
      subtitle = "Median Length Regional Post-Hoc Comparison"))


# Save
ggsave(
  plot = sc_p1_regionemmeans, 
  filename = here::here("Figs/small_community/sc_medlen_region_emmeans.png"))




##### d. Trend Posthoc  ####

# Slope Comparisons from 0

# Just survey area and year
(sloperegion_phoc_len <- emtrends(
  object = wig_len_mod, 
  specs =  ~ survey_area,
  var = "yr_num",
  adjust = "bonf"))


# Plotting Slope 
(sc_p1_yearemtrends <- sloperegion_phoc_len %>% 
    as_tibble() %>% 
    mutate(flag_alpha = ifelse(lower.CL > 0 | upper.CL <0, 1, 0.4)) %>% 
    ggplot(aes(survey_area, yr_num.trend, ymin = lower.CL, ymax = upper.CL)) +
    geom_hline(yintercept = 0, linetype = 2, color = "black") +
    geom_pointrange(aes(alpha = I(flag_alpha))) +
    # geom_pointrange(color = "gray30") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    labs(
      title = "Small Finfish Community",
      subtitle = "Median Length Annual Trend Coefficients",
      x = NULL,
      y = "Trend Coefficient"))


# Save
ggsave(
  plot = sc_p1_yearemtrends, 
  filename = here::here("Figs/small_community/sc_medlen_year_emtrends.png"))









###__ Pass 2: Covariates  __####





####  2. Median Length Model  ####


## Pass 2: No year term, introduce covariates

# No Interactions version - scaled btemp, log10 landings
# Singular fits for year intercepts b/c temp and landings are annual
wig_len_mod2 <- lm(
  formula = med_len_cm ~ survey_area + log10(landings) + scale(bot_temp),
  data = len_model_df)

# vif check - seems tolerable
plot(performance::check_collinearity(wig_len_mod2)) +
  coord_flip() +
  geom_hline(yintercept = 3, linetype = 3, aes(color = "Zuur Threshold"), linewidth = 1)


# Diagnostics
performance::check_model(wig_len_mod2)



# Summary 
tbl_regression(wig_len_mod2)  %>% 
  add_global_p(keep = T) %>% 
  bold_p(t = 0.05)  %>%
  bold_labels() %>%
  italicize_levels() 




##### b. Model Predictions  ####


# Plot marginal effects plots over observed data for:
# Bottom Temperature
bt_preds <- as.data.frame(
  ggpredict(wig_len_mod2, ~ bot_temp + survey_area) )


# Plotting over bottom temp differences
(sc_p2_btemp_region_margeffect <- bt_preds %>%
    mutate(survey_area = factor(group, levels = area_levels)) %>%
    ggplot() +
    geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = survey_area), alpha = 0.1) +
    geom_line(
      aes(x, predicted, group = survey_area),
      linewidth = 1) +
    geom_point(
      data = len_model_df,
      aes(bot_temp, med_len_cm, color = season),
      alpha = 0.6,
      size = 1) +
    scale_color_gmri() +
    facet_wrap(~survey_area) +
    scale_x_continuous(labels = scales::number_format(suffix = deg_c))+
    labs(
      y = "Median Length (cm)",
      x = "Bottom Temperature",
      title = "Small Finfish Community",
      subtitle = "Median Length and Bottom Temperature Marginal Mean Predictions"
    ))


# # Save
# ggsave(
#   plot = sc_p2_btemp_region_margeffect, 
#   filename = here::here("Figs/small_community/sc_medlen_btemp_region_margeffects.png"))




# Plot marginal effects plots over observed data for:
# Bottom Temperature
land_preds <- as.data.frame(
  ggpredict(
    model = wig_len_mod2,
    terms = list("landings" = 10^c(2:10), 
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
      data = len_model_df,
      aes(landings, med_len_cm, color = season),
      alpha = 0.6,
      size = 1) +
    scale_color_gmri() +
    facet_wrap(~survey_area) +
    scale_x_log10(labels = scales::label_log(10))+
    labs(
      y = "Median Length (cm)",
      x = "Total Landings (lb.)",
      title = "Small Finfish Community",
      subtitle = "Median Length and Lendings Marginal Mean Predictions"
    ))


# # Save
# ggsave(
#   plot = sc_p2_btemp_region_margeffect,
#   filename = here::here("Figs/small_community/sc_medlen_l10landings_region_margeffects.png"))





##### a.  Intercept Post-Hoc  ####

# Summary 
summary(wig_len_mod2)
# area
# season
# area:season



# Response scale
regseas_phoc_len2 <- emmeans(
  object = wig_len_mod2,
  specs = list(pairwise ~ survey_area),
  adjust = "tukey", 
  type = "response")


# Plot of the region & seasonal differences
(sc_p2_regionseason_emmeans <- regseas_phoc_len2$`emmeans of survey_area` %>%
    as_tibble() %>%
    ggplot(aes(survey_area, emmean, ymin = lower.CL, ymax = upper.CL)) +
    geom_pointrange(position = position_dodge(width = 0.25), alpha = 0.8) +
    scale_color_gmri() +
    labs(
      y = "Median Length (cm)",
      x = "Region",
      title = "Small Finfish Community",
      color = "Season",
      subtitle = "Median Length Region x Season Fixed Effects - Post-Hoc"))







##### b.  Trend Marginal Effects  ####

# No Significant Trend in Bot temp
emtrends(
  object = wig_len_mod2, 
  specs =  ~ survey_area,
  var = "bot_temp",
  adjust = "bonf")


# Just temp
# Plot marginal effects plots over observed data for:
# Bottom Temperature
(sc_p2_btemp_margeffect <- as.data.frame(
  ggpredict(wig_len_mod2, ~ bot_temp) ) %>%
    ggplot(aes(x, predicted, ymin = conf.low, ymax = conf.high)) +
    geom_ribbon(alpha = 0.1) +
    geom_line() +
    labs(
      y = "Median Length (cm)",
      x = "Bottom Temperature",
      title = "Small Finfish Community",
      subtitle = "Median Length and Bottom Temperature Marginal Mean Effect"
    ))


# Save
ggsave(
  plot = sc_p2_btemp_margeffect,
  filename = here::here("Figs/small_community/sc_medlen_btemp_margeffects.png"))




# Significant Trend in landings
emtrends(
  object = wig_len_mod2, 
  specs =  ~ survey_area,
  var = "landings",
  adjust = "bonf")


# Plot marginal effects plots over observed data for:
# Landings Effect
(sc_p2_land_margeffect <- as.data.frame(
  ggpredict(
    model = wig_len_mod2, 
    terms = list("landings" = 10^c(1:10)))) %>% 
    ggplot(aes(x, predicted, ymin = conf.low, ymax = conf.high)) +
    geom_ribbon(alpha = 0.1) +
    geom_line() +
    scale_x_log10(labels = scales::label_log(10),
                  limits = c(10^6, 10^9))+
    labs(
      y = "Median Length (cm)", 
      x = "Total Landings (lb.)",
      title = "Small Finfish Community",
      subtitle = "Median Length and Landings Marginal Mean Effect"
    ))


# Save
ggsave(
  plot = sc_p2_land_margeffect, 
  filename = here::here("Figs/small_community/sc_medlen_landings_margeffects.png"))



####  Export Results  ####
model_pack <- list(
  "pass_1" = wig_len_mod,
  "pass_2" = wig_len_mod2,
  "model_data" = len_model_df)

# Save them
saveRDS(object = model_pack, file = here::here("Data/models_and_results/sc_medlen.RDS"))




