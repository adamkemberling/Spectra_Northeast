#### Modeling Bodymass Distributions - Wigley ####


####  Packages  ####
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
wigley_medwt_df <- read_csv(here::here("Data/model_ready/wigley_community_medsize_mod.csv"))
wigley_bmspectra_df <- read_csv(here::here("Data/model_ready/wigley_community_bmspectra_mod.csv"))

# vectors for factor levels
area_levels <- c("GoM", "GB", "SNE", "MAB")
area_levels_long <- c("Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight")

# Degree symbol
deg_c <- "\u00b0C"


# Drop NA's
wt_model_df <- drop_na(wigley_medwt_df, total_weight_lb, bot_temp, med_wt_kg) %>% 
  mutate(yr_num = as.numeric(est_year),
         yr_fac = factor(est_year),
         survey_area = factor(survey_area, levels = area_levels),
         season = factor(season, levels = c("Spring", "Fall")),
         landings = total_weight_lb)



#### Colinearity Avoidance  ####

# Approach 1:
# Differencing of tthe two predictors to remove trend


# 1. weight model df
landings_temp_lm <- lm(landings ~ bot_temp, data = wt_model_df)
wt_model_df$land_resid <- resid(landings_temp_lm)
wt_model_df <- wt_model_df %>% 
  mutate(
    delta_bt = lag(bot_temp) - bot_temp,
    delta_land = lag(landings) - landings
  )




####  Quick Plots


# Median weight
ggplot(wt_model_df, aes(yr_num, med_wt_kg, color = season)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_ma(aes(linetype = "5-Year Moving Average"),n = 5, ma_fun = SMA) +
  geom_smooth(method = "lm", linewidth = 1, se = F,
              aes(linetype = "Regression Fit")) +
  facet_wrap(~survey_area, scales = "free") +
  scale_color_gmri() +
  labs(title = "Median Weight",
       subtitle = "Finfish Community",
       y = "Weight (kg)",
       x = "Year",
       color = "Season")




# Covariates over time time
wt_model_df %>% 
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
wt_model_df %>% 
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

####  1. Median Weight Model  ####


# Pass 1 Model
# wig_wt_mod <- glm(
#   formula = med_wt_kg ~ survey_area * yr_num * season,
#   family = gaussian(link = "log"),
#   data = wt_model_df)
wig_wt_mod <- lm(
  formula = log10(med_wt_kg) ~ survey_area * yr_num * season,
  data = wt_model_df)


# Check overall predictor significance

# This is how gtsummary is doing effect testing:
# car::Anova(wig_wt_mod)
tbl_regression(wig_wt_mod)  %>% 
  add_global_p(keep = T) %>% 
  bold_p(t = 0.05)  %>%
  bold_labels() %>%
  italicize_levels() 



# Model Diagnostics
check_model(wig_wt_mod)
check_outliers(wig_wt_mod)








##### b. Model Predictions  ####


# Plot the predictions over data


# Full Model
# Plot the predictions over data
wig_wt_mod_preds <- as.data.frame(
  ggpredict(wig_wt_mod, ~ yr_num + survey_area + season) )

# Plot over observed data
wig_wt_mod_preds %>% 
  mutate(
    survey_area = factor(group, levels = area_levels),
    season = factor(facet, levels = c("Spring", "Fall"))) %>% 
  ggplot() +
  geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = season), alpha = 0.1) +
  geom_line(
    aes(x, predicted, color = season), 
    linewidth = 1) +
  geom_point(
    data = wt_model_df,
    aes(yr_num, med_wt_kg, color = season),
    alpha = 0.4,
    size = 1) +
  facet_wrap(~survey_area, scales = "free") +
  scale_color_gmri() +
  labs(y = "Median Weight (kg)",
       title = "Small community, median weight",
       x = "Year")


##### c. Intercept Post-hoc  ####
# Use emmeans for post-hoc testing for factors

# Regions - significant
region_phoc_wt <- emmeans(
  wig_wt_mod, 
  list(pairwise ~ survey_area), 
  adjust = "tukey",
  type = "response")


# Custom Plot
(sc_p1_regionemmeans <- region_phoc_wt$`emmeans of survey_area` %>% 
    as_tibble() %>%
    ggplot(aes(survey_area, response, ymin = lower.CL, ymax = upper.CL)) +
    geom_pointrange(position = position_dodge(width = 0.25), alpha = 0.8) +
    scale_color_gmri() +
    labs(
      y = "EMMean - Median Weight (kg)",
      x = "Region",
      title = "Well-Studied Community",
      subtitle = "Median Weight Regional Post-Hoc Comparison"))



# Save
ggsave(
  plot = sc_p1_regionemmeans, 
  filename = here::here("Figs/small_community/sc_medwt_region_emmeans.png"))




##### d. Trend Posthoc  ####

# Slope Comparisons from 0

# Just survey area and year
(sloperegion_phoc_wt <- emtrends(
  object = wig_wt_mod, 
  specs =  ~ survey_area,
  var = "yr_num",
  adjust = "bonf"))


# Plotting Slope 
(sc_p1_yearemtrends <- sloperegion_phoc_wt %>% 
    as_tibble() %>% 
    mutate(flag_alpha = ifelse(lower.CL > 0 | upper.CL <0, 1, 0.4)) %>% 
    ggplot(aes(survey_area, yr_num.trend, ymin = lower.CL, ymax = upper.CL)) +
    geom_hline(yintercept = 0, linetype = 2, color = "black") +
    geom_pointrange(aes(alpha = I(flag_alpha))) +
    # geom_pointrange(color = "gray30") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    labs(
      title = "Well-Studied Community",
      subtitle = "Median Weight Annual Trend Coefficients",
      x = NULL,
      y = "Trend Coefficient"))


# Save
ggsave(
  plot = sc_p1_yearemtrends, 
  filename = here::here("Figs/small_community/sc_medwt_year_emtrends.png"))



# Pairwise testing of slopes:
(sloperegion_phoc_wt <- emtrends(
  object = wig_wt_mod, 
  specs =  pairwise ~ survey_area,
  var = "yr_num",
  adjust = "bonf"))







###__ Pass 2: Covariates  __####









####  2. Median Weight Model  ####



## Pass 2: No year term, introduce covariates

# No Interactions version - scaled btemp, log10 landings
# Singular fits for year intercepts b/c temp and landings are annual
wig_wt_mod2 <- lm(
  formula = log(med_wt_kg) ~ survey_area + log10(landings) + scale(bot_temp),
  data = wt_model_df)

# vif check - seems tolerable
plot(performance::check_collinearity(wig_wt_mod2)) +
  coord_flip() +
  geom_hline(yintercept = 3, linetype = 3, aes(color = "Zuur Threshold"), linewidth = 1)


# Summary 
gtsummary::tbl_regression(wig_wt_mod2)  %>% 
  add_global_p(keep = T) %>% 
  bold_p(t = 0.05)  %>%
  bold_labels() %>%
  italicize_levels() 



# Check performance
performance::check_model(wig_wt_mod2)



# Residuals Check
performance::check_autocorrelation(wig_wt_mod2)





##### a. Model Predictions  ####

# Not significant in log(weight) model

# Plot marginal effects plots over observed data for:

# Bottom Temperature
bt_preds <- as.data.frame(
  ggpredict(wig_wt_mod2, ~ bot_temp + survey_area) )


# Plotting over bottom temp differences
(sc_p2_btemp_region_margeffect <- bt_preds %>%
    mutate(survey_area = factor(group, levels = area_levels)) %>%
    ggplot() +
    geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high), alpha = 0.1) +
    geom_line(
      aes(x, predicted),
      linewidth = 1) +
    geom_point(
      data = wt_model_df,
      aes(bot_temp, med_wt_kg, color = season),
      alpha = 0.6,
      size = 1) +
    scale_color_gmri() +
    facet_wrap(~survey_area) +
    scale_x_continuous(labels = scales::number_format(suffix = deg_c))+
    labs(
      y = "Median Weight (kg)",
      x = "Bottom Temperature",
      title = "Well-Studied Community",
      subtitle = "Median Weight and Bottom Temperature Marginal Mean Predictions"
    ))


# Save
ggsave(
  plot = sc_p2_btemp_region_margeffect,
  filename = here::here("Figs/small_community/sc_medwt_btemp_region_margeffects.png"))




# Plot marginal effects plots over observed data for:
# Bottom Temperature
land_preds <- as.data.frame(
  ggpredict(
    model = wig_wt_mod2,
    terms = list("landings" = 10^c(2:10), 
                 "survey_area" = area_levels))
)



# Plotting over landings differences
(sc_p2_land_region_margeffect <- land_preds %>%
    mutate(
      survey_area = factor(group, levels = area_levels)) %>%
    ggplot() +
    geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high), alpha = 0.1) +
    geom_line(
      aes(x, predicted),
      linewidth = 1) +
    geom_point(
      data = wt_model_df,
      aes(landings, med_wt_kg, color = season),
      alpha = 0.6,
      size = 1) +
    scale_color_gmri() +
    facet_wrap(~survey_area) +
    scale_x_log10(labels = scales::label_log(10))+
    labs(
      y = "Median Weight (kg)",
      x = "Total Landings (lb.)",
      title = "Well-Studied Community",
      subtitle = "Median Weight and Lendings Marginal Mean Predictions"
    ))


# Save
ggsave(
  plot = sc_p2_btemp_region_margeffect,
  filename = here::here("Figs/small_community/sc_medwt_l10landings_region_margeffects.png"))





##### b.  Intercept Post-Hoc  ####

# Summary 
tbl_regression(wig_wt_mod2, exponentiate = F)
# area



# Region and Season Interaction - Significant

# emmean coefs
(regs_phoc_wt2 <- emmeans(
  object = wig_wt_mod2,
  specs = list(pairwise ~ survey_area),
  adjust = "tukey"))


# Response scale
(reg_phoc_wt2 <- emmeans(
  object = wig_wt_mod2,
  specs = list(pairwise ~ survey_area),
  adjust = "tukey", 
  type = "response"))


# Plot of the region & seasonal differences
(sc_p2_region_emmeans <- reg_phoc_wt2$`emmeans of survey_area` %>%
    as_tibble() %>%
    ggplot(aes(survey_area, response, ymin = lower.CL, ymax = upper.CL)) +
    geom_pointrange(position = position_dodge(width = 0.25), alpha = 0.8) +
    scale_color_gmri() +
    labs(
      y = "Median Weight (kg)",
      x = "Region",
      title = "Well-Studied Community",
      color = "Season",
      subtitle = "Median Weight Controlling for Temperature & Landings"))





##### b.  Trend Marginal Effects  ####
tbl_regression(wig_wt_mod2)  %>% 
  # add_global_p() %>% 
  bold_p(t = 0.10)  %>%
  bold_labels() %>%
  italicize_levels() 


# Just temp
# Plot marginal effects plots over observed data for:
# Bottom Temperature
(sc_p2_btemp_margeffect <- as.data.frame(
  ggpredict(wig_wt_mod2, ~ bot_temp) ) %>%
    ggplot(aes(x, predicted, ymin = conf.low, ymax = conf.high)) +
    geom_ribbon(alpha = 0.1) +
    geom_line() +
    labs(
      y = "Median Weight (kg)",
      x = "Bottom Temperature",
      title = "Well-Studied Community",
      subtitle = "Median Weight and Bottom Temperature Marginal Mean Effect"
    ))


# Save
ggsave(
  plot = sc_p2_btemp_margeffect,
  filename = here::here("Figs/small_community/sc_medwt_btemp_margeffects.png"))




# Significant Trend in landings


# Plot marginal effects plots over observed data for:
# Landings Effect
(sc_p2_land_margeffect <- as.data.frame(
  ggpredict(
    model = wig_wt_mod2, 
    terms = list("landings" = 10^c(1:10)))) %>% 
    ggplot(aes(x, predicted, ymin = conf.low, ymax = conf.high)) +
    geom_ribbon(alpha = 0.1) +
    geom_line() +
    scale_x_log10(labels = scales::label_log(10),
                  limits = c(10^1, 10^10))+
    labs(
      y = "Median Weight (kg)", 
      x = "Total Landings (lb.)",
      title = "Well-Studied Community",
      subtitle = "Median Weight and Landings Marginal Mean Effect"
    ))


# Save
ggsave(
  plot = sc_p2_land_margeffect,
  filename = here::here("Figs/small_community/sc_medwt_landings_margeffects.png"))






####  Export Results  ####
model_pack <- list(
  "pass_1" = wig_wt_mod,
  "pass_2" = wig_wt_mod2,
  "model_data" = wt_model_df)

# Save them
saveRDS(object = model_pack, file = here::here("Data/models_and_results/sc_medwt.RDS"))


