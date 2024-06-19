#####  Checking Landings and Lags  ####



#### Load existing models and Data  ####

# Median Length - Large Community
lc_len <- read_rds(here::here("Data/models_and_results/lc_medlen.RDS"))

# Length Spectra - Wigley Community
sc_lenb <- read_rds(here::here("Data/models_and_results/sc_lenspectra.RDS"))

# Median Weight - Small Community
sc_wt <- read_rds(here::here("Data/models_and_results/sc_medwt.RDS"))
