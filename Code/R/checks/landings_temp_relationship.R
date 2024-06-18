####  Whether or not to bother taking residuals from landings ~ bot_temp

# Note: the original intention of taking the residuals from this model
# instead of raw landings was to remove any shared linear trends with 
# time



#### Packages  ####
library(tidyquant)
library(tidyverse)
library(gmRi)
library(ggeffects)



#### Load Data  ####
ffish_medlen_df <- read_csv(here::here("Data/model_ready/large_community_medlength_mod.csv"))

# vectors for factor levels
area_levels <- c("GoM", "GB", "SNE", "MAB")

# Drop NA's
len_model_df <- drop_na(ffish_medlen_df, total_weight_lb, bot_temp, med_len_cm) %>% 
  mutate(yr_num = as.numeric(est_year),
         yr_fac = factor(est_year),
         survey_area = factor(survey_area, levels = area_levels),
         season = factor(season, levels = c("Spring", "Fall")))


####  Overview Plots  ####


# Median Length Changes
ggplot(len_model_df, aes(yr_num, med_len_cm, color = season)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_ma(aes(linetype = "5-Year Moving Average"),n = 5, ma_fun = SMA) +
  geom_smooth(method = "lm", linewidth = 1, se = T, 
              aes(linetype = "Regression Fit")) +
  facet_wrap(~survey_area) +
  scale_color_gmri() +
  labs(title = "Median Length",
       subtitle = "Broad Finfish Community",
       y = "Length (cm)",
       x = "Year",
       color = "Season")




#### Predictor Colinearity Avoidance  ####

# Approach 1:
# Proportion/residuals of landings not explained by bot temp

# 1. length model df
landings_temp_lm <- lm(total_weight_lb ~ bot_temp, data = len_model_df)
len_model_df$land_resid <- resid(landings_temp_lm)
len_model_df <- len_model_df %>% 
  group_by(survey_area, season) %>%
  arrange(yr_num) %>% 
  mutate(
    bt_diff = lag(bot_temp) - bot_temp,
    land_diff = lag(total_weight_lb) - total_weight_lb) %>% 
  ungroup()

# How different are they?

# almost no correlation R-squared:  0.002762 
summary(landings_temp_lm) 

# basically the same impact data
plot(len_model_df$land_resid, len_model_df$total_weight_lb,
     main = "Model Residuals and Original Values Highly Correlated") 




####  So these look nearly identical
# And would get nearly identical model results

# For scaling within dplyr
scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}



# Covariates / time
#### No adjustments, only scaling
len_model_df %>% 
  select(est_year, survey_area, total_weight_lb, bot_temp) %>% 
  mutate(bot_temp = scale_this(bot_temp),
         total_weight_lb = scale_this(total_weight_lb)) %>% 
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



# New Covariate / Time
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


# No not really, not needed


#####  Log or scale landings  ####





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