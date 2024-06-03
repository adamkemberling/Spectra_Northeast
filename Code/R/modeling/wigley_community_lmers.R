
#### Modeling Length Spectra  ####
library(lme4)
library(lmerTest)
library(emmeans)
library(merTools)
library(tidyquant)
library(tidyverse)
library(gmRi)
library(ggeffects)



#### Load Data  ####
wigley_medlen_df <- read_csv(here::here("Data/model_ready/wigley_community_medsize_mod.csv"))
wigley_lenspectra_df <- read_csv(here::here("Data/model_ready/wigley_community_lenspectra_mod.csv"))
wigley_bmspectra_df <- read_csv(here::here("Data/model_ready/wigley_community_bmspectra_mod.csv"))



# vectors for factor levels
area_levels <- c("GoM", "GB", "SNE", "MAB")
area_levels_long <- c("Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight")



# Drop NA's
len_model_df <- drop_na(wigley_medlen_df, total_weight_lb, bot_temp, med_len_cm) %>% 
  mutate(yr_num = as.numeric(est_year),
         yr_fac = factor(est_year),
         survey_area = factor(survey_area, levels = area_levels),
         season = factor(season, levels = c("Spring", "Fall")))
lenb_model_df <- drop_na(wigley_lenspectra_df, total_weight_lb, bot_temp, b) %>% 
  mutate(yr_num = as.numeric(est_year),
         yr_fac = factor(est_year),
         survey_area = factor(survey_area, levels = area_levels),
         season = factor(season, levels = c("Spring", "Fall")))
wtb_model_df <- drop_na(wigley_bmspectra_df, total_weight_lb, bot_temp, b) %>% 
  mutate(yr_num = as.numeric(est_year),
         yr_fac = factor(est_year),
         survey_area = factor(survey_area, levels = area_levels),
         season = factor(season, levels = c("Spring", "Fall")))



#### Colinearity Avoidance  ####

# Approach 1:
# Proportion/residuals of landings not explained by bot temp

# Approach 2:
# Differencing of tthe two predictors to remove trend


# 1. length model df
landings_temp_lm <- lm(total_weight_lb ~ bot_temp, data = len_model_df)
len_model_df$land_resid <- resid(landings_temp_lm)
len_model_df <- len_model_df %>% 
  mutate(
    bt_diff = lag(bot_temp),
    landings_diff = lag(total_weight_lb)
  )


# 2. length spectra model df
landings_temp_lm <- lm(total_weight_lb ~ bot_temp, data = lenb_model_df)
lenb_model_df$land_resid <- resid(landings_temp_lm)
lenb_model_df <- lenb_model_df %>% 
  mutate(
    bt_diff = lag(bot_temp),
    landings_diff = lag(total_weight_lb)
  )


# 3. bodymass spectra model df
landings_temp_lm <- lm(total_weight_lb ~ bot_temp, data = wtb_model_df)
wtb_model_df$land_resid <- resid(landings_temp_lm)
wtb_model_df <- wtb_model_df %>% 
  mutate(
    bt_diff = lag(bot_temp),
    landings_diff = lag(total_weight_lb)
  )



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

# Covariates 3/ time
len_model_df %>% 
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




####________________####
####  Pass 1: Change Over Time  ####


####____####
####  1. Median Length Model  ####

wig_len_mod <- lmerTest::lmer(
  formula = med_len_cm ~ survey_area * yr_num * season + ((1 | yr_fac)),
  data = len_model_df)



####___ a. Summary  ####

# Model Summary
summary(wig_len_mod)



# ggpredictions
simple_preds <- as.data.frame(
  ggpredict(wig_len_mod, ~ yr_num + season + survey_area) )

# Plot over observed data
simple_preds %>% 
  mutate(
    season = factor(group, levels = c("Spring", "Fall")),
    survey_area = factor(facet, levels = area_levels)) %>% 
  ggplot() +
  geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = season), alpha = 0.1) +
  geom_line(
    aes(x, predicted, color = season, linetype = "Predicted"), 
    linewidth = 1) +
  facet_wrap(~survey_area, scales = "free") +
  geom_point(
    data = len_model_df, 
    aes(yr_num, med_len_cm, color = season),
    alpha = 0.6) +
  geom_line(
    data = len_model_df,
    aes(yr_num, med_len_cm, color = season, linetype = "Observed"),
    alpha = 0.3, linewidth = 0.5) +
  scale_linetype_manual(values = c(3, 1)) +
  scale_color_gmri() +
  labs(y = "Median Length (cm)",
       title = "Wigley Community, median length")

####___ b. Intercept PH  ####

# Regions - significant
region_phoc_len<- emmeans(wig_len_mod, list(pairwise ~ survey_area), adjust = "tukey")
region_phoc_len
plot(region_phoc_len) + 
  coord_flip() +
  labs(y = "Region", 
       title = "Regional Intercept - Post-Hoc", 
       x = "Estimated Marginal Mean (Length)")



####___ c. Slope PH  ####

##. Just survey area and year
sloperegion_phoc_len <- emtrends(
  object = wig_len_mod, 
  specs =  ~ survey_area,
  var = "yr_num",
  adjust = "sidak")


# Plotting Slope 
sloperegion_phoc_len %>% 
  as_tibble() %>% 
  mutate(flag_alpha = ifelse(lower.CL > 0 | upper.CL <0, 1, 0.4)) %>% 
  ggplot(aes(survey_area, yr_num.trend, ymin = lower.CL, ymax = upper.CL)) +
  geom_hline(yintercept = 0, linetype = 1, color = "gray30", linewidth = 1) +
  geom_pointrange(aes(alpha = I(flag_alpha)), color = gmri_cols("blue"), size = 1) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(
    title = "Region Slopes",
    x = NULL,
    y = "Length (cm) ~ Year, EMTrend Coefficients")






####____####
####  2. Length Spectra Model  ####

wig_lenb_mod <- lmerTest::lmer(
  formula = b ~ survey_area * yr_num * season + ((1 | yr_fac)),
  data = lenb_model_df)


####___ a. Summary  ####

# Summary
summary(wig_lenb_mod)


# ggpredictions
simple_preds <- as.data.frame(
  ggpredict(wig_lenb_mod, ~ yr_num + season + survey_area) )

# Plot over observed data
simple_preds %>% 
  mutate(
    season = factor(group, levels = c("Spring", "Fall")),
    survey_area = factor(facet, levels = area_levels)) %>% 
  ggplot() +
  geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = season), alpha = 0.1) +
  geom_line(
    aes(x, predicted, color = season, linetype = "Predicted"), 
    linewidth = 1) +
  facet_wrap(~survey_area, scales = "free") +
  geom_point(
    data = lenb_model_df, 
    aes(yr_num, b, color = season),
    alpha = 0.6) +
  geom_line(
    data = lenb_model_df,
    aes(yr_num, b, color = season, linetype = "Observed"),
    alpha = 0.3, linewidth = 0.5) +
  scale_linetype_manual(values = c(3, 1)) +
  scale_color_gmri() +
  labs(y = "Median Length (cm)",
       title = "Wigley Community, Length Spectra")



####___ b. Intercept PH  ####

# Regions - significant
region_phoc_lenb <- emmeans(wig_lenb_mod, list(pairwise ~ survey_area), adjust = "tukey")
region_phoc_lenb
plot(region_phoc_lenb) + 
  coord_flip() +
  labs(y = "Region", 
       title = "Regional Intercept - Post-Hoc", 
       x = "Estimated Marginal Mean (Length Spectra)")


# Regions - significant
regseas_phoc_lenb <- emmeans(wig_lenb_mod, list(pairwise ~ survey_area *season), adjust = "tukey")
regseas_phoc_lenb
plot(regseas_phoc_lenb) + 
  coord_flip() +
  labs(y = "Region * Season", 
       title = "Regional Intercept - Post-Hoc", 
       x = "Estimated Marginal Mean (Length Spectra)")

####___ c. Slope PH  ####

##. Just survey area and year
sloperegion_phoc_lenb <- emtrends(
  object = wig_lenb_mod, 
  specs  =  ~ survey_area,
  var    = "yr_num",
  adjust = "sidak")


# Plotting Slope 
sloperegion_phoc_lenb %>% 
  as_tibble() %>% 
  mutate(flag_alpha = ifelse(lower.CL > 0 | upper.CL < 0, 1, 0.4)) %>% 
  ggplot(aes(survey_area, yr_num.trend, ymin = lower.CL, ymax = upper.CL)) +
  geom_hline(yintercept = 0, linetype = 1, color = "gray30", linewidth = 1) +
  geom_pointrange(aes(alpha = I(flag_alpha)), color = gmri_cols("blue"), size = 1) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(
    title = "Region Slopes",
    x = NULL,
    y = "Year slope EMTrend Coefficients")



####____####
####  3. Weight Spectra Model  ####


wig_wtb_mod <- lmerTest::lmer(
  formula = b ~ survey_area * yr_num * season + ((1 | yr_fac)),
  data = wtb_model_df)


####___ a. Summary  ####

# Summary
summary(wig_wtb_mod)


# ggpredictions
simple_preds <- as.data.frame(
  ggpredict(wig_wtb_mod, ~ yr_num + season + survey_area) )

# Plot over observed data
simple_preds %>% 
  mutate(
    season = factor(group, levels = c("Spring", "Fall")),
    survey_area = factor(facet, levels = area_levels)) %>% 
  ggplot() +
  geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = season), alpha = 0.1) +
  geom_line(
    aes(x, predicted, color = season, linetype = "Predicted"), 
    linewidth = 1) +
  facet_wrap(~survey_area, scales = "free") +
  geom_point(
    data = wtb_model_df, 
    aes(yr_num, b, color = season),
    alpha = 0.6) +
  geom_line(
    data = wtb_model_df,
    aes(yr_num, b, color = season, linetype = "Observed"),
    alpha = 0.3, linewidth = 0.5) +
  scale_linetype_manual(values = c(3, 1)) +
  scale_color_gmri() +
  labs(y = "Bodymass Spectra slope",
       title = "Wigley community, Bodymass Spectra")



####___ b. Intercept PH  ####

# Regions - significant
region_phoc_wtb <- emmeans(wig_wtb_mod, list(pairwise ~ survey_area), adjust = "tukey")
region_phoc_wtb
plot(region_phoc_wtb) + 
  coord_flip() +
  labs(y = "Region", 
       title = "Regional Intercept - Post-Hoc", 
       x = "Estimated Marginal Mean (Bodymass Spectra)")


# Regions * Seasons - significant
regseas_phoc_wtb <- emmeans(wig_wtb_mod, list(pairwise ~ survey_area *season), adjust = "tukey")
regseas_phoc_wtb
plot(regseas_phoc_wtb) + 
  coord_flip() +
  labs(y = "Region * Season", 
       title = "Regional Intercept - Post-Hoc", 
       x = "Estimated Marginal Mean (Bodymass Spectra)")



####___ c. Slope PH  ####

##. Just survey area and year
sloperegion_phoc_wtb <- emtrends(
  object = wig_wtb_mod, 
  specs =  ~ survey_area,
  var = "yr_num",
  adjust = "sidak")


# Plotting Slope 
sloperegion_phoc_wtb %>% 
  as_tibble() %>% 
  mutate(flag_alpha = ifelse(lower.CL > 0 | upper.CL <0, 1, 0.4)) %>% 
  ggplot(aes(survey_area, yr_num.trend, ymin = lower.CL, ymax = upper.CL)) +
  geom_hline(yintercept = 0, linetype = 1, color = "gray30", linewidth = 1) +
  geom_pointrange(aes(alpha = I(flag_alpha)), color = gmri_cols("blue"), size = 1) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(
    title = "Region Slopes",
    x = NULL,
    y = "Bodymass Spectra ~ Year, EMTrend Coefficients")




####________________####
####________________####
####________________####
####  Pass 2: Change w/ Covariates  ####


####____####
####  1. Median Length Model  ####

wig_len_mod2 <- lmerTest::lmer(
  formula = med_len_cm ~ scale(bot_temp) * survey_area * season + scale(land_resid) * survey_area + (1 | yr_fac),
  data = len_model_df)




####___ a. Summary  ####


# Summary 
summary(wig_len_mod2)


# Plot effects over observed data
bt_preds <- as.data.frame(
  ggpredict(wig_len_mod2, ~ bot_temp + season + survey_area) )

bt_preds %>% 
  mutate(
    season = factor(group, levels = c("Spring", "Fall")),
    survey_area = factor(facet, levels = area_levels)) %>% 
  ggplot() +
  geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = season), alpha = 0.1) +
  geom_line(
    aes(x, predicted, color = season, linetype = "Predicted"), 
    linewidth = 1) +
  facet_wrap(~survey_area, scales = "free") +
  geom_point(
    data = len_model_df, 
    aes(bot_temp, med_len_cm, color = season),
    alpha = 0.6) +
  scale_color_gmri() +
  labs(y = "Median Length (b)", x = "Bottom Temperature")


# Landings Residuals
land_preds <- as.data.frame(
  ggpredict(wig_len_mod2, ~ land_resid + season + survey_area) )

land_preds %>% 
  mutate(
    season = factor(group, levels = c("Spring", "Fall")),
    survey_area = factor(facet, levels = area_levels)) %>% 
  ggplot() +
  geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = season), alpha = 0.1) +
  geom_line(
    aes(x, predicted, color = season, linetype = "Predicted"), 
    linewidth = 1) +
  facet_wrap(~survey_area, scales = "free") +
  geom_point(
    data = len_model_df, 
    aes(land_resid, med_len_cm, color = season),
    alpha = 0.6) +
  scale_color_gmri() +
  labs(y = "Median Length (cm)", x = "resid(landings ~ bot_temp)")




####___ b. Intercept PH  ####

# Regions - significant
region_phoc_len2 <- emmeans(wig_len_mod2, list(pairwise ~ survey_area), adjust = "tukey")
region_phoc_len2
plot(region_phoc_len2) + 
  coord_flip() +
  labs(y = "Region", 
       title = "Regional Intercept - Post-Hoc", 
       x = "Estimated Marginal Mean (median length)")


# Regions * Seasons - significant
regseas_phoc_len2 <- emmeans(wig_len_mod2, list(pairwise ~ survey_area *season), adjust = "tukey")
regseas_phoc_len2
plot(regseas_phoc_len2) + 
  coord_flip() +
  labs(y = "Region * Season", 
       title = "Regional Intercept - Post-Hoc", 
       x = "Estimated Marginal Mean (median length)")


####___ c. Slope PH  ####

# Landings residuals
slope_phoc_len2 <- emtrends(
  object = wig_len_mod2, 
  specs =  ~ survey_area,
  var = "land_resid",
  adjust = "sidak")

slope_phoc_len2 %>% 
  as_tibble() %>% 
  mutate(flag_alpha = ifelse(lower.CL > 0 | upper.CL <0, 1, 0.3)) %>% 
  ggplot(aes(survey_area, land_resid.trend, ymin = lower.CL, ymax = upper.CL)) +
  geom_hline(yintercept = 0, linetype = 1, color = "gray30", linewidth = 1) +
  geom_pointrange(aes(alpha = I(flag_alpha)), color = gmri_cols("blue"), size = 1) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(
    title = "Wigley Median Length Slopes",
    x = NULL,
    y = "resid(landings) Trend")





####____####
####  2. Length Spectra Model  ####
wig_lenb_mod2 <- lmerTest::lmer(
  formula = b ~ scale(bot_temp) * survey_area * season + scale(land_resid) * survey_area + (1 | yr_fac),
  data = lenb_model_df)




####___ a. Summary  ####


# Summary 
summary(wig_lenb_mod2)


# Plot effects over observed data
bt_preds <- as.data.frame(
  ggpredict(wig_lenb_mod2, ~ bot_temp + season + survey_area) )

bt_preds %>% 
  mutate(
    season = factor(group, levels = c("Spring", "Fall")),
    survey_area = factor(facet, levels = area_levels)) %>% 
  ggplot() +
  geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = season), alpha = 0.1) +
  geom_line(
    aes(x, predicted, color = season, linetype = "Predicted"), 
    linewidth = 1) +
  facet_wrap(~survey_area, scales = "free") +
  geom_point(
    data = lenb_model_df, 
    aes(bot_temp, b, color = season),
    alpha = 0.6) +
  scale_color_gmri() +
  labs(y = "Length Spectra Slope (b)", x = "Bottom Temperature")


# Landings Residuals
land_preds <- as.data.frame(
  ggpredict(wig_lenb_mod2, ~ land_resid + season + survey_area) )

land_preds %>% 
  mutate(
    season = factor(group, levels = c("Spring", "Fall")),
    survey_area = factor(facet, levels = area_levels)) %>% 
  ggplot() +
  geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = season), alpha = 0.1) +
  geom_line(
    aes(x, predicted, color = season, linetype = "Predicted"), 
    linewidth = 1) +
  facet_wrap(~survey_area, scales = "free") +
  geom_point(
    data = lenb_model_df, 
    aes(land_resid, b, color = season),
    alpha = 0.6) +
  scale_color_gmri() +
  labs(y = "Length Spectra Slope (b)", x = "resid(landings ~ bot_temp)")



####___ b. Intercept PH  ####


# Regions - significant
region_phoc_lenb2 <- emmeans(wig_lenb_mod2, list(pairwise ~ survey_area), adjust = "tukey")
region_phoc_lenb2
plot(region_phoc_lenb2) + 
  coord_flip() +
  labs(y = "Region", 
       title = "Regional Intercept - Post-Hoc", 
       x = "Estimated Marginal Mean (Length Spectra)")


# Regions * Seasons - significant
regseas_phoc_lenb2 <- emmeans(wig_lenb_mod2, list(pairwise ~ survey_area *season), adjust = "tukey")
regseas_phoc_lenb2
plot(regseas_phoc_len2) + 
  coord_flip() +
  labs(y = "Region * Season", 
       title = "Regional Intercept - Post-Hoc", 
       x = "Estimated Marginal Mean (Length Spectra)")




####___ c. Slope PH  ####


# Landings residuals
slope_phoc_lenb2 <- emtrends(
  object = wig_lenb_mod2, 
  specs =  ~ survey_area,
  var = "land_resid",
  adjust = "sidak")

slope_phoc_lenb2 %>% 
  as_tibble() %>% 
  mutate(flag_alpha = ifelse(lower.CL > 0 | upper.CL <0, 1, 0.3)) %>% 
  ggplot(aes(survey_area, land_resid.trend, ymin = lower.CL, ymax = upper.CL)) +
  geom_hline(yintercept = 0, linetype = 1, color = "gray30", linewidth = 1) +
  geom_pointrange(aes(alpha = I(flag_alpha)), color = gmri_cols("blue"), size = 1) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(
    title = "Wigley Length Spectra Slopes",
    x = NULL,
    y = "resid(landings) Trend")



####____####
####  3. Bodymass Spectra Model  ####

wig_wtb_mod2 <- lmerTest::lmer(
  formula = b ~ scale(bot_temp) * survey_area * season + scale(land_resid) * survey_area + (1 | yr_fac),
  data = wtb_model_df)


# Summary 
summary(wig_wtb_mod2)


# Plot effects over observed data
bt_preds <- as.data.frame(
  ggpredict(wig_wtb_mod2, ~ bot_temp + season + survey_area) )

bt_preds %>% 
  mutate(
    season = factor(group, levels = c("Spring", "Fall")),
    survey_area = factor(facet, levels = area_levels)) %>% 
  ggplot() +
  geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = season), alpha = 0.1) +
  geom_line(
    aes(x, predicted, color = season, linetype = "Predicted"), 
    linewidth = 1) +
  facet_wrap(~survey_area, scales = "free") +
  geom_point(
    data = wtb_model_df, 
    aes(bot_temp, b, color = season),
    alpha = 0.6) +
  scale_color_gmri() +
  labs(y = "Bodymass Spectra Slope (b)", x = "Bottom Temperature")


# Landings Residuals
land_preds <- as.data.frame(
  ggpredict(wig_wtb_mod2, ~ land_resid + season + survey_area) )

land_preds %>% 
  mutate(
    season = factor(group, levels = c("Spring", "Fall")),
    survey_area = factor(facet, levels = area_levels)) %>% 
  ggplot() +
  geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = season), alpha = 0.1) +
  geom_line(
    aes(x, predicted, color = season, linetype = "Predicted"), 
    linewidth = 1) +
  facet_wrap(~survey_area, scales = "free") +
  geom_point(
    data = wtb_model_df, 
    aes(land_resid, b, color = season),
    alpha = 0.6) +
  scale_color_gmri() +
  labs(y = "Bodymass Spectra Slope (b)", x = "resid(landings ~ bot_temp)")


####___ b. Intercept PH  ####


# Regions - significant
region_phoc_wtb2 <- emmeans(wig_wtb_mod2, list(pairwise ~ survey_area), adjust = "tukey")
region_phoc_wtb2
plot(region_phoc_wtb2) + 
  coord_flip() +
  labs(y = "Region", 
       title = "Regional Intercept - Post-Hoc", 
       x = "Estimated Marginal Mean (Bodymass Spectra)")


# Regions * Seasons - significant
regseas_phoc_wtb2 <- emmeans(wig_wtb_mod2, list(pairwise ~ survey_area *season), adjust = "tukey")
regseas_phoc_wtb2
plot(regseas_phoc_wtb2) + 
  coord_flip() +
  labs(y = "Region * Season", 
       title = "Regional Intercept - Post-Hoc", 
       x = "Estimated Marginal Mean (Bodymass Spectra)")



####___ b. Slope PH  ####


# Landings residuals
landings_phoc_wtb2 <- emtrends(
  object = wig_wtb_mod2, 
  specs =  ~ survey_area,
  var = "land_resid",
  adjust = "sidak")

landings_phoc_wtb2 %>% 
  as_tibble() %>% 
  mutate(flag_alpha = ifelse(lower.CL > 0 | upper.CL <0, 1, 0.3)) %>% 
  ggplot(aes(survey_area, land_resid.trend, ymin = lower.CL, ymax = upper.CL)) +
  geom_hline(yintercept = 0, linetype = 1, color = "gray30", linewidth = 1) +
  geom_pointrange(aes(alpha = I(flag_alpha)), color = gmri_cols("blue"), size = 1) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(
    title = "Wigley Bodymass Spectra Slopes",
    x = NULL,
    y = "resid(landings) Trend")



# Bottom Temperatures
btemp_phoc_wtb2 <- emtrends(
  object = wig_wtb_mod2, 
  specs =  ~ survey_area,
  var = "bot_temp",
  adjust = "sidak")

btemp_phoc_wtb2 %>% 
  as_tibble() %>% 
  mutate(flag_alpha = ifelse(lower.CL > 0 | upper.CL <0, 1, 0.3)) %>% 
  ggplot(aes(survey_area, bot_temp.trend, ymin = lower.CL, ymax = upper.CL)) +
  geom_hline(yintercept = 0, linetype = 1, color = "gray30", linewidth = 1) +
  geom_pointrange(aes(alpha = I(flag_alpha)), color = gmri_cols("blue"), size = 1) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(
    title = "Wigley Bodymass Spectra Slopes",
    x = NULL,
    y = "resid(landings) Trend")
akemberling