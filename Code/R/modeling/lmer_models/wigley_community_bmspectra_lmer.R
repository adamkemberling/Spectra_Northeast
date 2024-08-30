
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

theme_set(theme_gmri(rect = element_rect(fill = "white", color = NA)))

# Degree symbol
deg_c <- "\u00b0C"




# vectors for factor levels
area_levels <- c("GoM", "GB", "SNE", "MAB")
area_levels_long <- c("Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight")



#### Load Data  ####
wigley_bmspectra_df <- read_csv(here::here("Data/model_ready/wigley_community_bmspectra_mod.csv")) 



# Drop NA's
wtb_model_df <- drop_na(wigley_bmspectra_df, total_weight_lb, bot_temp, b) %>%
  group_by(survey_area, season) %>%
  #group_by(survey_area) %>%
  mutate(roll_temp = zoo::rollapply(bot_temp, 5, mean, na.rm = T, align = "right",  fill = NA),
         .groups = "drop") %>% 
  mutate(yr_num = as.numeric(est_year),
         yr_fac = factor(est_year),
         survey_area = factor(survey_area, levels = area_levels),
         season = factor(season, levels = c("Spring", "Fall")),
         landings = total_weight_lb,
         yr_seas = str_c(season, est_year))







####______####
##### Key Plots for Paper  ####


#### Bottom Temperature Plot ####

# model for significant trends
temp_mod <- lm(bot_temp ~ est_year * survey_area * season,
               data = wtb_model_df)



# Get predictions, filter out non-zero trends from post-hoc comparison
# Plot the predictions over data
temp_preds <- as.data.frame(
  ggpredict(
    temp_mod, 
    terms = ~ est_year*survey_area*season) ) %>% 
  mutate(
    survey_area = factor(group, levels = area_levels),
    season = factor(facet, levels = c("Spring", "Fall")))




# Drop effect fits that are non-significant  ###
# Just survey area and year
temp_emtrend <- emtrends(
  object = temp_mod,
  specs =  ~ survey_area*season,
  var = "est_year") %>% 
  as_tibble() %>% 
  mutate(zero = 0,
         non_zero = if_else(between(zero, lower.CL, upper.CL), F, T))



# Plot the temperature changes and the significant trends
temp_preds %>% 
  left_join(select(temp_emtrend, season, survey_area, non_zero)) %>% 
  mutate(season = factor(season, levels = c("Spring", "Fall"))) %>% 
  filter(non_zero) %>% 
  ggplot() +
  geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = survey_area), alpha = 0.1) +
  geom_line(
    aes(x, predicted, color = survey_area), 
    linewidth = 1) +
  geom_line(
    data = wtb_model_df,
    aes(est_year, bot_temp, color = survey_area), alpha = 0.4, linewidth = 0.5) +
  facet_grid(survey_area~season, scales = "free") +
  scale_color_gmri() +
  labs(
    title = "Seasonal Bottom Temperature Trends",
    subtitle = "Significant Annual Trend Fits Shown",
    y = "Seasonal Bottom Temperature")







#### For context, what   ####


# b ~ temp
ggplot(wtb_model_df) +
  geom_point(aes(bot_temp, b, color = survey_area)) +
  facet_grid(.~season) +
  scale_color_gmri() 


# b ~ landings
ggplot(wtb_model_df) +
  geom_point(aes(log10(total_weight_lb), b, color = survey_area)) +
  facet_grid(survey_area~season) +
  scale_color_gmri() 


# What do they look like on the same x axis?
wtb_model_df %>% 
  select(survey_area, season, est_year, b, total_weight_lb, bot_temp) %>% 
  pivot_longer(cols = c(b, total_weight_lb, bot_temp), names_to = "var", values_to = "val") %>% 
  ggplot(aes(est_year, val)) +
  geom_point(aes(color = survey_area), alpha = 0.35, size = 0.5) +
  geom_ma(aes(color = survey_area), n = 5, ma_fun = SMA, linetype = 1) +
  facet_grid(var~season, scales = "free") +
  scale_color_gmri() +
  theme(legend.position = "bottom") +
  labs(color = "5-Year Moving Average",
       y = "Value",
       x = "Year")
  







####____ Key Models  ####





#### Model 1:  Trends in Time  ####
# mod1 <- lmer(
#   formula = b ~ est_year*survey_area*season + (1|est_year), 
#   data = wtb_model_df)
# summary(mod1)
mod1_c <- lmer(
  formula = b ~ est_year*survey_area + (1|yr_seas), 
  data = wtb_model_df)
summary(mod1_c)

# Same stories
plot(ggeffects::ggpredict(mod1, terms = ~est_year*survey_area*season), add.data = T )





# Plot the predictions over data
mod1_preds <- as.data.frame(
  ggpredict(mod1, terms = ~ est_year*survey_area*season) ) %>% 
  mutate(
    survey_area = factor(group, levels = area_levels),
    season = factor(facet, levels = c("Spring", "Fall")))


# Drop effect fits that are non-significant  ###
# Just survey area and year
mod1_emtrend <- emtrends(
  object = mod1,
  specs =  ~ survey_area*season,
  var = "est_year") %>% 
  as_tibble() %>% 
  mutate(zero = 0,
         non_zero = if_else(between(zero, lower.CL, upper.CL), F, T))


# Plot over observed data
mod1_preds  %>% 
  left_join(select(mod1_emtrend, survey_area, season, non_zero)) %>% 
  filter(non_zero) %>% 
  ggplot() +
  geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = season), alpha = 0.1) +
  geom_line(
    aes(x, predicted, color = season), 
    linewidth = 1) +
  geom_line(
    data = wtb_model_df,
    aes(yr_num, b, color = season),
    alpha = 0.4,
    linewidth = 0.5) +
  facet_wrap(~survey_area, scales = "free") +
  scale_color_gmri() +
  labs(y = "Body Mass Spectra Slope (b)",
       title = "Wigley Species, Body Mass Spectra",
       x = "Year")




# Plot over observed data
mod1_preds  %>% 
  left_join(select(mod1_emtrend, survey_area, season, non_zero)) %>% 
  filter(non_zero) %>% 
  ggplot() +
  geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = survey_area), alpha = 0.1) +
  geom_line(
    aes(x, predicted, color = survey_area), 
    linewidth = 1) +
  geom_line(
    data = wtb_model_df,
    aes(yr_num, b, color = survey_area),
    alpha = 0.4,
    linewidth = 0.5) +
  facet_grid(survey_area~season, scales = "free") +
  scale_color_gmri() +
  labs(y = "Body Mass Spectra Slope (b)",
       title = "Wigley Species, Body Mass Spectra",
       x = "Year")



#### Model 1b. year_seas ####

mod1 <- lmer(
  formula = b ~ est_year*survey_area*season + (1|yr_seas), 
  data = wtb_model_df)
summary(mod1)


# Same stories
plot(ggeffects::ggpredict(mod1, terms = ~est_year*survey_area*season), add.data = T )





# Plot the predictions over data
mod1_preds <- as.data.frame(
  ggpredict(mod1, terms = ~ est_year*survey_area*season) ) %>% 
  mutate(
    survey_area = factor(group, levels = area_levels),
    season = factor(facet, levels = c("Spring", "Fall")))


# Drop effect fits that are non-significant  ###
# Just survey area and year
mod1_emtrend <- emtrends(
  object = mod1,
  specs =  ~ survey_area*season,
  var = "est_year") %>% 
  as_tibble() %>% 
  mutate(zero = 0,
         non_zero = if_else(between(zero, lower.CL, upper.CL), F, T))

# # Same significance results
# trend_df <- emtrends(
#   object = mod1,  
#   ~survey_area*season, 
#   var = "est_year")
# summary(trend_df, infer= c(T,T))




# Plot over observed data
mod1_preds  %>% 
  left_join(select(mod1_emtrend, survey_area, season, non_zero)) %>% 
  filter(non_zero) %>% 
  ggplot() +
  geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = season), alpha = 0.1) +
  geom_line(
    aes(x, predicted, color = season), 
    linewidth = 1) +
  geom_line(
    data = wtb_model_df,
    aes(yr_num, b, color = season),
    alpha = 0.4,
    linewidth = 0.5) +
  facet_wrap(~survey_area, scales = "free") +
  scale_color_gmri() +
  labs(y = "Body Mass Spectra Slope (b)",
       title = "Wigley Species, Body Mass Spectra",
       x = "Year")




# Plot over observed data
mod1_preds  %>% 
  left_join(select(mod1_emtrend, survey_area, season, non_zero)) %>% 
  filter(non_zero) %>% 
  ggplot() +
  geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = survey_area), alpha = 0.1) +
  geom_line(
    aes(x, predicted, color = survey_area), 
    linewidth = 1) +
  # geom_line(
  #   data = wtb_model_df,
  #   aes(yr_num, b, color = survey_area),
  #   alpha = 0.4,
  #   linewidth = 0.5) +
  geom_point(
    data = wtb_model_df,
    aes(yr_num, b, color = survey_area), alpha = 0.35, size = 0.8) +
  facet_grid(survey_area~season, scales = "free") +
  scale_color_gmri() +
  labs(y = "Body Mass Spectra Slope (b)",
       title = "Wigley Species, Body Mass Spectra",
       x = "Year")






#### Model 2: Temperature & Fishing ####

# actual mod
mod2 <- lmer(
  b ~ survey_area*season*scale(roll_temp) + survey_area*scale(log(total_weight_lb)) + (1|est_year), 
  data = wtb_model_df)


summary(mod2)
# check_model(mod2)



# Quick plot
plot(ggeffects::ggpredict(mod2, terms = ~roll_temp*survey_area*season), add.data = T)




# Clean up the plot:
# Plot the predictions over data
mod2_preds <- as.data.frame(
  ggpredict(
    mod2, 
    terms = ~ roll_temp*survey_area*season))   %>% 
  mutate(
    survey_area = factor(group, levels = area_levels),
    season = factor(facet, levels = c("Spring", "Fall")))



# #### Mod 2: Trend Posthoc  ####
trend_df <- emtrends(
  object = mod2,
  ~survey_area * season,
  var = "roll_temp")
summary(trend_df, infer= c(T,T))


# Drop effect fits that are non-significant  ###
mod2_emtrend <- emtrends(
  object = mod2,
  specs =  ~ survey_area*season,
  var = "roll_temp") %>%
  as_tibble() %>%
  mutate(
    zero = 0,
    non_zero = if_else(between(zero, lower.CL, upper.CL), F, T))




# Limit temperature plotting range
actual_range <- wtb_model_df %>% 
  group_by(season, survey_area) %>% 
  summarise(min_temp = min(bot_temp)-2,
            max_temp = max(bot_temp)+2,
            .groups = "drop")


  

# Plot over observed data

mod2_preds %>% 
  left_join(actual_range) %>% 
  filter((x < min_temp) == F,
         (x > max_temp) == F) %>% 
  left_join(select(mod2_emtrend, survey_area, season, non_zero)) %>%
  filter(non_zero) %>%
  ggplot() +
  geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = survey_area), alpha = 0.1) +
  geom_line(
    aes(x, predicted, color = survey_area), 
    linewidth = 1) +
  geom_point(
    data = wtb_model_df,
    aes(roll_temp, b, color = survey_area),
    alpha = 0.4,
    size = 1) +
  facet_grid(survey_area~season, scales = "free") +
  scale_color_gmri() +
  labs(y = "Body Mass Spectra Slope (b)",
       title = "Wigley Species, Body Mass Spectra",
       x = "5-Year Rolling Average Bottom Temperature")






#### Model 2b: Temperature & Fishing + year_seas ####

# actual mod
mod2 <- lmer(
  b ~ survey_area*scale(roll_temp) + survey_area*scale(log(total_weight_lb)) + (1|yr_seas), 
  data = wtb_model_df)


# summary(mod2)
# check_model(mod2)



# Quick plot
# plot(ggeffects::ggpredict(mod2, terms = ~roll_temp*survey_area*season), add.data = T)
plot(ggeffects::ggpredict(mod2, terms = ~roll_temp*survey_area), add.data = T) +
  scale_color_gmri()



# Clean up the plot:
# Plot the predictions over data
mod2_preds <- as.data.frame(
  ggpredict(
    mod2, 
    terms = ~ roll_temp*survey_area,
    # terms = list(
    #   "roll_temp" = seq(3,20,.2)), 
    #   "survey_area" = factor(area_levels, levels = area_levels)
  )) %>% 
  mutate(
    survey_area = factor(group, levels = area_levels)  )

# Trend Posthoc 
trend_df <- emtrends(
  object = mod2,  
  ~survey_area, 
  var = "roll_temp")
summary(trend_df, infer= c(T,T))


# Drop effect fits that are non-significant  ###
mod2_emtrend <- emtrends(
  object = mod2,
  specs =  ~ survey_area,
  var = "roll_temp") %>%
  as_tibble() %>%
  mutate(
    zero = 0,
    non_zero = if_else(between(zero, lower.CL, upper.CL), F, T))


# Limit temperature plotting range
actual_range <- wtb_model_df %>% 
  group_by(survey_area) %>% 
  summarise(min_temp = min(roll_temp, na.rm = T)-3,
            max_temp = max(roll_temp, na.rm = T)+3,
            .groups = "drop")


# Plot over observed data
mod2_preds %>% 
  left_join(actual_range) %>%
  filter((x < min_temp) == F,
         (x > max_temp) == F) %>%
  left_join(select(mod2_emtrend, survey_area, non_zero)) %>%
  filter(non_zero) %>%
  ggplot() +
  geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = survey_area), alpha = 0.1) +
  geom_line(
    aes(x, predicted, color = survey_area), 
    linewidth = 1) +
  geom_point(
    data = wtb_model_df,
    aes(roll_temp, b, color = survey_area, shape = season),
    alpha = 0.6,
    size = 1.2) +
  facet_grid(survey_area~.) +
  scale_color_gmri() +
  scale_x_continuous(labels = label_number(suffix = deg_c)) +
  labs(y = "Body Mass Spectra Slope (b)",
       title = "Wigley Species, Body Mass Spectra",
       shape = "Season",
       color = "Area",
       x = "5-Year Rolling Bottom Temperature")










####______####
####______####
####______####
####______####
####______####













####_______________________#####


#### Old Stuff  ####

####_______________________#####
#### Colinearity Avoidance  ####

# Approach 1:
# Proportion/residuals of landings not explained by bot temp

# Approach 2:
# Differencing of tthe two predictors to remove trend


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

# Covariateschange over time
wtb_model_df %>% 
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

####____________________####

####__  Pass 1: Change/Time  __####

####  1. Body Mass Spectra Model  ####
wig_wtb_mod <- lmerTest::lmer(
  formula = b ~ survey_area * yr_num * season + (1 | yr_fac),
  data = wtb_model_df)



# Check important predictors
summary(wig_wtb_mod)
# survey area and year, not season


# Diagnostics
check_model(wig_wtb_mod)
check_outliers(wig_wtb_mod)

# Calculate intra-class correlation
icc(wig_wtb_mod)
# 4.1% of total variance is attributable to between year variance






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

# # Regions - significant
# region_phoc_wtb <- emmeans(
#   wig_wtb_mod, 
#   list(pairwise ~ survey_area), 
#   adjust = "tukey",
#   type = "response")
# 
# 
# # Custom Plot
# (sc_p1_regionemmeans <- region_phoc_wtb$`emmeans of survey_area` %>% 
#     as_tibble() %>%
#     ggplot(aes(survey_area, emmean, ymin = lower.CL, ymax = upper.CL)) +
#     geom_pointrange(position = position_dodge(width = 0.25), alpha = 0.8) +
#     scale_color_gmri() +
#     labs(
#       y = "EMMean - Body Mass Spectra Slope (b)",
#       x = "Region",
#       title = "Small Finfish Community",
#       subtitle = "Body Mass Spectra Regional Post-Hoc Comparison"))
# 
# 
# # Save
# ggsave(
#   plot = sc_p1_regionemmeans, 
#   filename = here::here("Figs/small_community/sc_wtb_region_emmeans.png"))

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
#   adjust = "sidak"))
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
  adjust = "sidak"))


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
wig_wtb_mod2 <- lmerTest::lmer(
  formula = b ~ survey_area * season + log10(landings) + scale(bot_temp) + (1 | yr_fac),
  data = wtb_model_df)

# vif check - seems tolerable
plot(performance::check_collinearity(wig_wtb_mod2)) +
  coord_flip() +
  geom_hline(yintercept = 3, linetype = 3, aes(color = "Zuur Threshold"), linewidth = 1)


# Summary 
summary(wig_wtb_mod2)
# log10(landings)
# survey area
# season
# season:area



performance::check_model(wig_wtb_mod2)




# Calculate intra-class correlation
# ratio of the random intercept variance (between year variance)
# to the total variance
RandomEffects <- as.data.frame(VarCorr(wig_wtb_mod2))
ICC_between <- RandomEffects[1,4]/(RandomEffects[1,4]+RandomEffects[2,4]) 
ICC_between
# 1.2%




##### b. Model Predictions  ####

# Not significant in log(weight) model

# Plot marginal effects plots over observed data for:

# Not significant
# # Bottom Temperature
# bt_preds <- as.data.frame(
#   ggpredict(wig_wtb_mod2, ~ bot_temp + survey_area + season) )
# 
# 
# # Plotting over bottom temp differences
# (sc_p2_btemp_region_margeffect <- bt_preds %>%
#     mutate(
#       season = factor(facet, levels = c("Spring", "Fall")),
#       survey_area = factor(group, levels = area_levels)) %>%
#     ggplot() +
#     geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = season), alpha = 0.1) +
#     geom_line(
#       aes(x, predicted, color = season),
#       linewidth = 1) +
#     geom_point(
#       data = wtb_model_df,
#       aes(bot_temp, b, color = season),
#       alpha = 0.6,
#       size = 1) +
#     scale_color_gmri() +
#     facet_wrap(~survey_area) +
#     scale_x_continuous(labels = scales::number_format(suffix = deg_c))+
#     labs(
#       y = "Body Mass Spectra Slope (b)",
#       x = "Bottom Temperature",
#       title = "Small Finfish Community",
#       subtitle = "Body Mass Spectra and Bottom Temperature Marginal Mean Predictions"
#     ))


# # Save
# ggsave(
#   plot = sc_p2_btemp_region_margeffect, 
#   filename = here::here("Figs/small_community/sc_wtb_btemp_regseas_margeffects.png"))




# Plot marginal effects plots over observed data for:
# Bottom Temperature
land_preds <- as.data.frame(
  ggpredict(
    model = wig_wtb_mod2,
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


# # Save
# ggsave(
#   plot = sc_p2_btemp_region_margeffect,
#   filename = here::here("Figs/small_community/sc_wtb_l10landings_regseas_margeffects.png"))





##### a.  Intercept Post-Hoc  ####

# Summary 
summary(wig_wtb_mod2)
# area
# season
# area:season


# Region and Season Interaction - Significant

# emmean coefs
regseas_phoc_wtb2 <- emmeans(
  object = wig_wtb_mod2,
  specs = list(pairwise ~ survey_area * season),
  adjust = "tukey")

regseas_phoc_wtb2


# Response scale
regseas_phoc_wtb2 <- emmeans(
  object = wig_wtb_mod2,
  specs = list(pairwise ~ survey_area * season),
  adjust = "tukey", 
  type = "response")


# Plot of the region & seasonal differences
(sc_p2_regionseason_emmeans <- regseas_phoc_wtb2$`emmeans of survey_area, season` %>%
    as_tibble() %>%
    ggplot(aes(survey_area, emmean, ymin = lower.CL, ymax = upper.CL)) +
    geom_pointrange(aes(color = season), position = position_dodge(width = 0.25), alpha = 0.8) +
    scale_color_gmri() +
    labs(
      y = "Body Mass Spectra Slope (b)",
      x = "Region",
      title = "Small Finfish Community",
      color = "Season",
      subtitle = "Body Mass Spectra Region x Season Fixed Effects - Post-Hoc"))


# Save
ggsave(
  plot = sc_p2_regionseason_emmeans, 
  filename = here::here("Figs/small_community/sc_wtb_regseason_emmeans.png"))




##### b.  Trend Marginal Effects  ####

# # No Significant Trend in Bot temp
# emtrends(
#   object = wig_wtb_mod2, 
#   specs =  ~ survey_area,
#   var = "bot_temp",
#   adjust = "sidak")
# 
# 
# # Just temp
# # Plot marginal effects plots over observed data for:
# # Bottom Temperature
# (sc_p2_btemp_margeffect <- as.data.frame(
#   ggpredict(wig_wtb_mod2, ~ bot_temp) ) %>%
#     ggplot(aes(x, predicted, ymin = conf.low, ymax = conf.high)) +
#     geom_ribbon(alpha = 0.1) +
#     geom_line() +
#     labs(
#       y = "Body Mass Spectra Slope (b)",
#       x = "Bottom Temperature",
#       title = "Small Finfish Community",
#       subtitle = "Body Mass Spectra and Bottom Temperature Marginal Mean Effect"
#     ))
# 
# 
# # Save
# ggsave(
#   plot = sc_p2_btemp_margeffect,
#   filename = here::here("Figs/small_community/sc_wtb_btemp_margeffects.png"))




# Significant Trend in landings
emtrends(
  object = wig_wtb_mod2, 
  specs =  ~ survey_area,
  var = "landings",
  adjust = "sidak")


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




