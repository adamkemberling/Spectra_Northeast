# Building Expectations for Arrhenius equation / TSR
# Credit to Dr. Bart Difiore for the code-up


library(tidyverse)
library(gmRi)
library(patchwork)
library(scales)



# Set ggplot theme for figures
theme_set(
  theme_classic() + 
    theme(
      # Titles
      plot.title = element_text(hjust = 0, face = "bold", size = 14),
      plot.subtitle = element_text(size = 12),
      plot.caption = element_text(size = 7.2, margin = margin(t = 20), color = "gray40"),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 7.5),
      # Axes
      axis.line.y = element_line(color = "black"),
      axis.ticks.y = element_line(), 
      axis.line.x = element_line(color = "black"),
      axis.ticks.x = element_line(), 
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      rect = element_rect(fill = "transparent", color = "black"),
      # Facets
      strip.text = element_text(color = "white", 
                                face = "bold",
                                size = 11),
      strip.background = element_rect(
        color = "#00736D", 
        fill = "#00736D", 
        size = 1, 
        linetype="solid"))
)


####  Metabolic Theory of Ecology  ####

"LINDMARK 2018: Across species, rates are often
assumed and found to scale as power functions 
of mass with exponents of 3/4 for whole organism rates, 
exponentially with temperature, and with 
independent mass and temperature effects (e.g., in
the Arrhenius fractal supply model applied in the metabolic theory
of ecology (MTE) (Brown et al., 2004; Downs et al., 2008; Gillooly
et al., 2001))."



#### Allometric scaling relationship B ~ M  ####
# Metabolism B as a function of body size (M)
# Overall scaling exponent 0.75 or 3/4
# Estimated in West et al. 1994

# Range of masses: m
m <- seq(2,1000)

# age 0 mass: a0_m
a0_m <- 1

# Metabolism Scales as a function of mody mass
B <- a0_m*m^0.75

# Plot
data.frame("B" = B, "m" = m) %>% 
  ggplot(aes(m, B)) +
  geom_line() +
  labs(
    x = "Mass",
    y = "Metabolism (B)",
    title = "Allometric Scaling: Resting Metabolism ~ Body Mass"
  )




####  Arrhenius equation  ####
# The Arrhenius equation is a formula for the temperature dependence of reaction rates
# Higher temperatures increase reaction rates exponentially
ar <- function(r0 = 10, E = 0.67, k = 8.62*10^-5, temp){
  r0*exp(-1*E/(k*temp))
}

# For a range of temperatures we can predict reaction rates
temp = seq(10,30, length.out = 100)+273.15
R <- ar(temp = temp)

data.frame("R" = R, "Temp" = temp) %>% 
  ggplot(aes(Temp, R)) +
  geom_line() +
  labs(
    title = "Arrhenius equation - Reaction Rate ~ Temperature", 
    x = "Temperature (k)", 
    y = "Reaction Rate")





####  Applications: Temperature Size Rule  ####

# Temperature Size Rule

# First half is power-law relationship of meass dependency (scaling) 
# Second half of this equation: is arrhenius equation
temp_size <- function(
    # normalization constant (sometimes called b0), fit empirically at the level of class (mammal, reptile)
    # the temperature component of the equation the Universal Temperature Dependence of Metabolism (UTD)
    # Changes the intercept, varies within and across taxonomic groups
    r0 = 5, 
    #‘mean activation energy of metabolism’, 
    # its value being estimated empirically from measurements of enzyme kinetics 
    # in vitro
    E = 0.67, 
    # Boltzmanns constant
    k = 8.62*10^-5, 
    temp,  # temperature (k)
    mass   # Mass
    ){
  
  # the equation
  r0*mass^0.75 * exp(-1*E/(k*temp))
}



# Set up range of mass for 3 different temperatures
df <- expand.grid(
  temp = c(10, 20, 30)+273.15, 
  mass = seq(2,10000))

# Estimate size-based metabolism
df$metabolism <- temp_size(temp = df$temp, mass = df$mass)




# Mass-specific metabolism
ggplot(df, aes(x = mass, y = metabolism/mass))+
  geom_line(aes(color = as.factor(temp-273.15)))+
  scale_x_log10(labels = label_log(base = 10))+
  scale_y_log10(labels = label_log(base = 10)) +
  labs(
    color = "Temp C", 
    y = "Mass-Normalized Metabolism\n(metabolism / mass)", 
    x = "Mass",
    title = "TSR - Higher Temperatures Increase the Mass-Specific Metabolism")






####  Applications: Sharpe Schoolfield Equation  ####

# Not original citation:
# https://onlinelibrary.wiley.com/doi/full/10.1046/j.1420-9101.2001.00272.x

"
The temperature dependence of any biological rate has a roughly triangular 
shape (Huey & Kingsolver, 1989; Kingsolver & Woods, 1997). The Sharpe–Schoolfield 
equation yields this familiar and general shape accurately (Fig.1C). 
Wagner et.al. (1984) examined the utility of the Sharpe–Schoolfield equation to 
describe biological rates – mostly development rate of insects – and found the model to give an extremely accurate description of development rates."


# Sharpe Schoolfield equation
sharpe_schoolfield <- function(
    c0   = 0.79, 
    E    = 0.73, 
    Eh   = 1.89, 
    k    = 8.62*10^-5, 
    Tc   = 273.15, 
    Th   = 25+273.15, 
    beta = 0.75, 
    temp, 
    mass){
  
  #
  numerator = c0*mass^beta*exp(E*(1/(k*Tc) - 1/(k*temp)))
  denominator = 1 + exp(Eh*(1/(k*Th) - 1/(k*temp)))
  
  numerator/denominator
  
}



# Make plot on 
df <- expand.grid(temp = c(10+273.15, 30+273.15), mass = seq(2,10))

df$metab <- sharpe_schoolfield(temp = df$temp, mass = df$mass)


p1 <- df %>%
  mutate(temp_cat = ifelse(temp == 10+273.15, "cold", "warm")) %>%
  ggplot(aes(x = mass, y = metab/mass))+
  geom_line(aes(color = temp_cat))+
  scale_color_manual(values = c("blue", "red"))+
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  labs(title = "Mass-Specific Metabolism Higher\nUnder Warmer Temperatures")

# Make another plot showing how thermal preferences scale with size
# using sharpe schoolfield
df <- expand.grid(
  temp = seq(10,30)+273.15, 
  mass = c(2,5,10))  %>%
  mutate(mass_cat = case_when(
    mass == 2 ~ "small", 
    mass == 5 ~ "medium", 
    mass == 10 ~ "large"))

# Estimate metabolism from temp and body mass
df$metab <- sharpe_schoolfield(temp = df$temp, mass = df$mass)

p2 <- df %>%
  mutate(mass_cat = case_when(
    mass == 2 ~ "small", 
    mass == 5 ~ "medium", 
    mass == 10 ~ "large")) %>%
  ggplot(aes(x = temp-273.15, y = metab))+
  geom_line(aes(color = mass_cat))+
  scale_color_manual(values = c("darkgreen", "green", "lightgreen"))+
  theme_classic()  +
  labs(title = "Thermal Sensitivity Higher\nin Larger Individuals")


# Put these together
p1 | p2




####  Applications of Theory:   ####


# Make another plot showing how Mass-specific metabolism relates to metabolism
# using sharpe schoolfield

p3 <- df %>%
  ggplot(aes(x = temp-273.15, y = metab/mass))+
  geom_line(aes(color = mass_cat))+
  scale_color_manual(values = c("darkgreen", "green", "lightgreen"))+
  theme_classic() +
  labs(title = "When Accounting for Mass,\nSomething...")


(p1 / p2 / p3 )




#### Relative Energy Reuirements  ####




df <- expand.grid(
  temp = c(19+273.15), 
  mass = seq(2,10, length.out = 100))


df$I_normal <- 1 * df$mass^0.75
df$I_warm <- 1 * df$mass^0.5
df$metabolism <- 0.5 * df$mass^0.75


df %>%
  pivot_longer(cols = I_normal:metabolism) %>%
  ggplot(aes(x = mass, y = value/mass))+
    geom_line(aes(color = name), linewidth = 1)+
  scale_x_log10()+
  scale_y_log10()+
  scale_color_manual(values = c("blue", "red", "gray50"), labels = c("Ingestion (Normal)", "Ingestion (Warm)", "Standard metabolic demand"))+
  labs(x = "Mass", y = "Energy gains or requirements", color = "")+
  theme_classic()+
  theme(legend.position = c(0.8,0.9))











# # random check of sst before I pass it to a student...
# gom_sst <- oisst_access_timeseries("gmri", "apershing gulf of maine", "cloudstorage")
# nrow(gom_sst)
# length(unique(gom_sst$time))
# 
# # Subset date and temperature and save
# gom_sst_lite <- select(gom_sst, date = time, sea_surface_temp_c = area_wtd_sst) %>% 
#   filter(date < as.Date("2024-09-01"))
# write_csv(gom_sst_lite, here::here("stash", "GMRI_GulfOfMaine_SST_OISSTv2.csv"))
