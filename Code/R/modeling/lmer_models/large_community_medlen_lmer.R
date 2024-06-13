
#### Modeling Median Length - Large Community Length Spectra  ####

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

#### Load Data  ####
ffish_medlen_df <- read_csv(here::here("Data/model_ready/large_community_medlength_mod.csv"))
ffish_lenspectra_df <- read_csv(here::here("Data/model_ready/large_community_lenspectra_mod.csv"))

# vectors for factor levels
area_levels <- c("GoM", "GB", "SNE", "MAB")
area_levels_long <- c("Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight")

# Degree symbol
deg_c <- "\u00b0C"


# Drop NA's
len_model_df <- drop_na(ffish_medlen_df, total_weight_lb, bot_temp, med_len_cm) %>% 
  mutate(yr_num = as.numeric(est_year),
         yr_fac = factor(est_year),
         survey_area = factor(survey_area, levels = area_levels),
         season = factor(season, levels = c("Spring", "Fall")),
         landings = total_weight_lb)
lenb_model_df <- drop_na(ffish_lenspectra_df, total_weight_lb, bot_temp, b) %>% 
  mutate(yr_num = as.numeric(est_year),
         yr_fac = factor(est_year),
         survey_area = factor(survey_area, levels = area_levels),
         season = factor(season, levels = c("Spring", "Fall")),
         landings = total_weight_lb)




####  Overview Plots  ####


# Median Length
ggplot(len_model_df, aes(yr_num, med_len_cm, color = season)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_ma(aes(linetype = "5-Year Moving Average"),n = 5, ma_fun = SMA) +
  geom_smooth(method = "lm", linewidth = 1, se = T, 
              aes(linetype = "Regression Fit")) +
  facet_wrap(~survey_area) +
  scale_color_gmri() +
  labs(title = "Median Length",
       subtitle = "Finfish Community",
       y = "Length (cm)",
       x = "Year",
       color = "Season")


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







#### Predictor Colinearity Avoidance  ####

# Approach 1:
# Proportion/residuals of landings not explained by bot temp

# Approach 2:
# Differencing of the two predictors to remove trend

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}


# 1. length model df
landings_temp_lm <- lm(landings ~ bot_temp, data = len_model_df)
len_model_df$land_resid <- resid(landings_temp_lm)
len_model_df <- len_model_df %>% 
  group_by(survey_area, season) %>%
  arrange(yr_num) %>% 
  mutate(
    delta_bt = lag(bot_temp) - bot_temp,
    delta_land = lag(landings) - landings) %>% 
  ungroup()

# How different are they?
#summary(landings_temp_lm) # almost no correlation
# # plot(len_model_df$land_resid, len_model_df$landings) # basically the same impact

# 2. length spectra model df
landings_temp_lm <- lm(landings ~ bot_temp, data = lenb_model_df)
lenb_model_df$land_resid <- resid(landings_temp_lm)
lenb_model_df <- lenb_model_df %>% 
  group_by(survey_area, season) %>%
  arrange(yr_num) %>% 
  mutate(
    delta_bt = lag(bot_temp) - bot_temp,
    delta_land = lag(landings) - landings) %>% 
  ungroup()



####  Why do these look identical


# Covariates / time
#### No adjustments, only scaling
len_model_df %>% 
  select(est_year, survey_area, landings, bot_temp) %>% 
  mutate(bot_temp = scale_this(bot_temp),
         landings = scale_this(landings)) %>% 
  pivot_longer(
    -c(est_year, survey_area), 
    names_to = "Covariate", 
    values_to = "value") %>% 
  ggplot(aes(est_year, value, color = Covariate)) +
  geom_point(size = 1, alpha = 0.2) +
  #geom_ma(aes(linetype = "5-Year Moving Average"),n = 5, ma_fun = SMA) +
  geom_smooth(method = "lm", linewidth = 1, se = T, 
              aes(linetype = "Regression Fit")) +
  facet_wrap(~survey_area) +
  scale_color_gmri() +
  labs(title = "Trends in scaled predictors",
       x = "Year",
       color = "Covariate")

# # Cor.test - 0.07
# cor.test(len_model_df$landings, len_model_df$bot_temp)


#### Approach 1: resid(landings ~ bot_temp)
# Do we get different plots with landings as residuals?
len_model_df %>% 
  select(est_year, survey_area, land_resid, bot_temp) %>% 
  mutate(
    bot_temp = scale_this(bot_temp),
    land_resid = scale_this(land_resid)) %>% 
  pivot_longer(
    -c(est_year, survey_area), 
    names_to = "Covariate", 
    values_to = "value") %>% 
  ggplot(aes(est_year, value, color = Covariate)) +
  geom_point(size = 1, alpha = 0.2) +
  #geom_ma(aes(linetype = "5-Year Moving Average"),n = 5, ma_fun = SMA) +
  geom_smooth(method = "lm", linewidth = 1, se = T, 
              aes(linetype = "Regression Fit")) +
  facet_wrap(~survey_area) +
  scale_color_gmri() +
  labs(title = "Trends in resid(landings ~ bot_temp) _ scale()",
       x = "Year",
       color = "Covariate")

# wth! - model removed like no original variability
# summary(landings_temp_lm)
# # Cor.test - 
# cor.test(scale(len_model_df$land_resid), scale(len_model_df$bot_temp))
# 


# Differencing
#### Approach 2: diff(landings) & diff(bot_temp)
len_model_df %>% 
  select(est_year, survey_area, delta_land, delta_bt) %>% 
  mutate(
    delta_bt = scale_this(delta_bt),
    delta_land = scale_this(delta_land)) %>% 
  pivot_longer(-c(est_year, survey_area), 
               names_to = "Covariate", 
               values_to = "value") %>% 
  ggplot(aes(est_year, value, color = Covariate)) +
  geom_point(size = 1, alpha = 0.2) +
  #geom_ma(aes(linetype = "5-Year Moving Average"),n = 5, ma_fun = SMA) +
  geom_smooth(method = "lm", linewidth = 1, se = T, 
              aes(linetype = "Regression Fit")) +
  facet_wrap(~survey_area) +
  scale_color_gmri() +
  labs(title = "Exploring effect of differencing",
       x = "Year",
       color = "Covariate")

# Not correlated - cool
#cor.test(len_model_df$delta_bt, len_model_df$delta_land)


####  Log(landings) vs. scale(landings)

# scale
hist(scale(lenb_model_df$landings)) # Raw -hist
skewness(scale(lenb_model_df$landings)) # Raw - skew
range(scale(lenb_model_df$landings)) # Raw - range
hist(log10(lenb_model_df$landings)) # Log10 - hist
skewness(log10(lenb_model_df$landings)) # Log10 - skew
range(log10(lenb_model_df$landings)) # Log10 - range

# Both, loses interpretability, same skew
skewness(scale(log10(lenb_model_df$landings)))
hist(scale(log10(lenb_model_df$landings)))




# Check VIF before any dummy-variables are included -- Seem Fine < 1.2
# https://statisticalhorizons.com/multicollinearity/
collin_lm <- lm(med_len_cm ~ log10(landings) + scale(bot_temp) + yr_num, data = len_model_df)
check_collinearity(collin_lm)
collin_blm <- lm(b ~ log10(landings) + scale(bot_temp) + yr_num, data = lenb_model_df)
check_collinearity(collin_blm)

# going rogue - probably fine actually
# simple lm
rogue_len_mod <- lm(
  formula = med_len_cm ~ survey_area * yr_num * season + scale(bot_temp) + log10(landings),
  data = len_model_df)
check_model(rogue_len_mod)
summary(rogue_len_mod)

# lmer
rogue_len_lmer <- lmerTest::lmer(
  formula = med_len_cm ~ survey_area * yr_num * season + scale(bot_temp) + log10(landings) +
    (1 | yr_fac),
  data = len_model_df)
check_model(rogue_len_lmer)
icc(rogue_len_lmer)
summary(rogue_len_lmer)


####____________________####

####__  Pass 1: Change/Time  __####

####  1. Median Length Model  ####
ffish_len_mod <- lmerTest::lmer(
  formula = med_len_cm ~ survey_area * yr_num * season + (1 | yr_fac),
  # formula = log(med_len_cm) ~ survey_area * yr_num * season + (1 | yr_fac),
  data = len_model_df)


# Check important predictors
summary(ffish_len_mod)
check_model(ffish_len_mod)

# survey area and year, not season


# Calculate intra-class correlation
# https://quantdev.ssri.psu.edu/tutorials/r-bootcamp-introduction-multilevel-model-and-interactions
# Store the random effect variances, which will be the first column of the VarCorr object
RandomEffects <- as.data.frame(VarCorr(ffish_len_mod))
RandomEffects

# Next, compute the ICC. It is the ratio of the random intercept variance (between-year var) over the total variance (between + within var):
ICC_between <- RandomEffects[1,4] / (RandomEffects[1,4] + RandomEffects[2,4]) 
ICC_between
icc(ffish_len_mod)
# 6.5% of total variance is attributable to between year variance






##### b. Model Predictions  ####


# Plot the predictions over data
ffish_len_mod_preds <- as.data.frame(
  ggpredict(ffish_len_mod, ~ yr_num + survey_area) )

# Plot over observed data
ffish_len_mod_preds %>% 
  mutate(survey_area = factor(group, levels = area_levels)) %>% 
  ggplot() +
  geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = survey_area), alpha = 0.1) +
  geom_line(
    aes(x, predicted, color = survey_area), 
    linewidth = 1) +
  geom_point(
    data = len_model_df,
    aes(yr_num, med_len_cm, color = survey_area),
    alpha = 0.4,
    size = 1) +
  facet_wrap(~survey_area, scales = "free") +
  scale_color_gmri() +
  labs(y = "Median Length (cm)",
       title = "large community, median length",
       x = "Year")




##### c. Intercept Post-hoc  ####
# Use emmeans for post-hoc testing for factors

# Regions - significant
(region_phoc_len <- emmeans(
  ffish_len_mod, 
  list(pairwise ~ survey_area), 
  adjust = "tukey"))


# Custom Plot
(lc_p1_regionemmeans <- region_phoc_len$`emmeans of survey_area` %>% 
  as_tibble() %>%
  ggplot(aes(survey_area, emmean, ymin = lower.CL, ymax = upper.CL)) +
  geom_pointrange(position = position_dodge(width = 0.25), alpha = 0.8) +
  scale_color_gmri() +
  labs(
    y = "EMMean - Median Length (cm)",
    x = "Region",
    title = "Large Finfish Community",
    subtitle = "Median Length Regional Post-Hoc Comparison"))


# Save
ggsave(
  plot = lc_p1_regionemmeans, 
  filename = here::here("Figs/large_community/lc_medlen_region_emmeans.png"))




##### d. Trend Posthoc  ####

# Slope Comparisons from 0

# Just survey area and year
(sloperegion_phoc_len <- emtrends(
  object = ffish_len_mod, 
  specs =  ~ survey_area,
  var = "yr_num",
  adjust = "sidak"))


# Plotting Slope 
(lc_p1_yearemtrends <- sloperegion_phoc_len %>% 
  as_tibble() %>% 
  mutate(flag_alpha = ifelse(lower.CL > 0 | upper.CL <0, 1, 0.4)) %>% 
  ggplot(aes(survey_area, yr_num.trend, ymin = lower.CL, ymax = upper.CL)) +
  geom_hline(yintercept = 0, linetype = 2, color = "black") +
  geom_pointrange(aes(alpha = I(flag_alpha))) +
  # geom_pointrange(color = "gray30") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(
    title = "Large Finfish Community",
    subtitle = "Median Length Annual Trend Coefficients",
    x = NULL,
    y = "Trend Coefficient"))


# Save
ggsave(
  plot = lc_p1_yearemtrends, 
  filename = here::here("Figs/large_community/lc_medlen_year_emtrends.png"))






####  Evaluating Pass2 Options  ####

# # # Original version - was using landings residuals
# ffish_len_mod2 <- lmerTest::lmer(
#   formula = med_len_cm ~ scale(bot_temp) * survey_area * season + scale(land_resid) * survey_area + (1 | yr_fac),
#   data = len_model_df)
# 
# # vif check - looks bad, so I'm questioning why use resid()...
# plot(performance::check_collinearity(ffish_len_mod2)) +
#   coord_flip() +
#   geom_hline(yintercept = 3, linetype = 3, aes(color = "Zuur Threshold"), linewidth = 1)


# Revisiting whether we can include all 3, Can we have our cake & eat it? 
# - with simple models yes? complex models no?
# - continue with the 2 pass approach

# # Complex version - scaled actual predictors, everything interacts
# ffish_len_mod2 <- lmerTest::lmer(
#   formula = med_len_cm ~  yr_num * survey_area * season * scale(bot_temp) * scale(landings)  + (1 | yr_fac),
#   data = len_model_df)
# 
# # vif check - seems horrible
# # fails without year as well
# plot(performance::check_collinearity(ffish_len_mod2)) +
#   coord_flip() +
#   geom_hline(yintercept = 3, linetype = 3, aes(color = "Zuur Threshold"), linewidth = 1)
# 
# 
# # Simple version - scaled actual predictors
# ffish_len_mod2 <- lmerTest::lmer(
#   formula = med_len_cm ~  yr_num * survey_area * season + scale(bot_temp) + scale(landings)  + (1 | yr_fac),
#   data = len_model_df)
# 
# # vif check - seems fine for the things we're focused on
# plot(performance::check_collinearity(ffish_len_mod2)) +
#   coord_flip() +
#   geom_hline(yintercept = 3, linetype = 3, aes(color = "Zuur Threshold"), linewidth = 1)


# # Interactions version - scaled actual predictors
# # Terrible VIF
# ffish_len_mod2 <- lmerTest::lmer(
#   formula = med_len_cm ~ survey_area:scale(bot_temp)  + survey_area:scale(landings) +
#     scale(bot_temp) + scale(landings) + (1 | yr_fac),
#   data = len_model_df)
# 
# # vif check - seems not good
# plot(performance::check_collinearity(ffish_len_mod2)) +
#   coord_flip() +
#   geom_hline(yintercept = 3, linetype = 3, aes(color = "Zuur Threshold"), linewidth = 1)


# # Limited Interactions version - scaled actual predictors
# ffish_len_mod2 <- lmerTest::lmer(
#   formula = med_len_cm ~ survey_area:scale(landings) +
#     scale(bot_temp) + scale(landings) + (1 | yr_fac),
#   data = len_model_df)
# 
# # vif check - seems still bad
# plot(performance::check_collinearity(ffish_len_mod2)) +
#   coord_flip() +
#   geom_hline(yintercept = 3, linetype = 3, aes(color = "Zuur Threshold"), linewidth = 1)


# # No Interactions version - scaled actual predictors
# # Singular fits for year intercepts b/c temp and landings are annual
# ffish_len_mod2 <- lmerTest::lmer(
#   formula = med_len_cm ~ survey_area * season + scale(landings) + scale(bot_temp) + (1 | yr_fac),
#   data = len_model_df)
# 
# # vif check - seems tolerable
# plot(performance::check_collinearity(ffish_len_mod2)) +
#   coord_flip() +
#   geom_hline(yintercept = 3, linetype = 3, linewidth = 1)




###__ Pass 2: Covariates  __####





####  2. Median Length Model  ####


## Pass 2: No year term, introduce covariates

# No Interactions version - scaled btemp, log10 landings
# Singular fits for year intercepts b/c temp and landings are annual
ffish_len_mod2 <- lmerTest::lmer(
  formula = med_len_cm ~ survey_area * season + log10(landings) + scale(bot_temp) + (1 | yr_fac),
  data = len_model_df)

# vif check - seems tolerable
plot(performance::check_collinearity(ffish_len_mod2)) +
  coord_flip() +
  geom_hline(yintercept = 3, linetype = 3, aes(color = "Zuur Threshold"), linewidth = 1)


# Summary 
summary(ffish_len_mod2)
# scale(bot_temp)
# log10(landings), but not scale(landings)
# survey area
# season
# season:area



# Calculate intra-class correlation
# ratio of the random intercept variance (between year variance)
# to the total variance
RandomEffects <- as.data.frame(VarCorr(ffish_len_mod2))
ICC_between <- RandomEffects[1,4]/(RandomEffects[1,4]+RandomEffects[2,4]) 
ICC_between
# 5.5% very low




##### b. Model Predictions  ####


# Plot marginal effects plots over observed data for:
# Bottom Temperature
bt_preds <- as.data.frame(
  ggpredict(ffish_len_mod2, ~ bot_temp + survey_area + season) )


# Plotting over bottom temp differences
(lc_p2_btemp_region_margeffect <- bt_preds %>% 
    mutate(
      season = factor(facet, levels = c("Spring", "Fall")),
      survey_area = factor(group, levels = area_levels)) %>% 
    ggplot() +
    geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = season), alpha = 0.1) +
    geom_line(
      aes(x, predicted, color = season), 
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
      title = "Large Finfish Community",
      subtitle = "Median Length and Bottom Temperature Marginal Mean Predictions"
    ))


# # Save
# ggsave(
#   plot = lc_p2_btemp_region_margeffect, 
#   filename = here::here("Figs/large_community/lc_medlen_btemp_region_margeffects.png"))




# Plot marginal effects plots over observed data for:
# Bottom Temperature
land_preds <- as.data.frame(
  ggpredict(
    model = ffish_len_mod2,
    #terms = ~ landings + survey_area + season)
    terms = list("landings" = 10^c(1:10), 
                 "survey_area" = area_levels,
                 season = c("Spring", "Fall")))
  )



# Plotting over landings differences
(lc_p2_land_region_margeffect <- land_preds %>%
    mutate(
      season = factor(facet, levels = c("Spring", "Fall")),
      survey_area = factor(group, levels = area_levels)) %>%
    ggplot() +
    geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = season), alpha = 0.1) +
    geom_line(
      aes(x, predicted, color = season),
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
      title = "Large Finfish Community",
      subtitle = "Median Length and Lendings Marginal Mean Predictions"
    ))


# # Save
# ggsave(
#   plot = lc_p2_btemp_region_margeffect,
#   filename = here::here("Figs/large_community/lc_medlen_l10landings_region_margeffects.png"))





##### a.  Intercept Post-Hoc  ####

# Summary 
summary(ffish_len_mod2)
# area
# season
# area:season


# Region and Season Interaction - Significant
regseas_phoc_len2 <- emmeans(
  object = ffish_len_mod2,
  specs = list(pairwise ~ survey_area * season),
  adjust = "tukey")
regseas_phoc_len2


# Plot of the region & seasonal differences
(lc_p2_regionseason_emmeans <- regseas_phoc_len2$`emmeans of survey_area, season` %>%
  as_tibble() %>%
  ggplot(aes(survey_area, emmean, ymin = lower.CL, ymax = upper.CL)) +
  geom_pointrange(aes(color = season), position = position_dodge(width = 0.25), alpha = 0.8) +
  scale_color_gmri() +
  labs(
    y = "Median Length (cm)",
    x = "Region",
    title = "Large Finfish Community",
    color = "Season",
    subtitle = "Median Length Region x Season Fixed Effects - Post-Hoc"))


# Save
ggsave(
  plot = lc_p2_regionseason_emmeans, 
  filename = here::here("Figs/large_community/lc_medlen_regionseason_emmeans.png"))




##### b.  Trend Marginal Effects  ####

# Significant Trend in Bot temp
emtrends(
  object = ffish_len_mod2, 
  specs =  ~ survey_area,
  var = "bot_temp",
  adjust = "sidak")


# Just temp
# Plot marginal effects plots over observed data for:
# Bottom Temperature
(lc_p2_btemp_margeffect <- as.data.frame(
  ggpredict(ffish_len_mod2, ~ bot_temp) ) %>% 
  ggplot(aes(x, predicted, ymin = conf.low, ymax = conf.high)) +
  geom_ribbon(alpha = 0.1) +
  geom_line() +
  labs(
    y = "Median Length (cm)", 
    x = "Bottom Temperature",
    title = "Large Finfish Community",
    subtitle = "Median Length and Bottom Temperature Marginal Mean Effect"
  ))


# Save
ggsave(
  plot = lc_p2_btemp_margeffect, 
  filename = here::here("Figs/large_community/lc_medlen_btemp_margeffects.png"))




# Significant Trend in landings
emtrends(
  object = ffish_len_mod2, 
  specs =  ~ survey_area,
  var = "landings",
  adjust = "sidak")


# Plot marginal effects plots over observed data for:
# Landings Effect
(lc_p2_land_margeffect <- as.data.frame(
  ggpredict(
    model = ffish_len_mod2, 
    terms = list("landings" = 10^c(1:10)))) %>% 
    ggplot(aes(x, predicted, ymin = conf.low, ymax = conf.high)) +
    geom_ribbon(alpha = 0.1) +
    geom_line() +
    scale_x_log10(labels = scales::label_log(10),
                  limits = c(10^1, 10^10))+
    labs(
      y = "Median Length (cm)", 
      x = "Total Landings (lb.)",
      title = "Large Finfish Community",
      subtitle = "Median Length and Landings Marginal Mean Effect"
    ))


# Save
ggsave(
  plot = lc_p2_land_margeffect, 
  filename = here::here("Figs/large_community/lc_medlen_landings_margeffects.png"))





