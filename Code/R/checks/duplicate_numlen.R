# Checking duplicate length + numlen records for species
# Are they always in here? yea to some extent---

# Load data from local directory
#Path where "survdat" & spp_keys.csv data exists
data_path <- here::here("Data/")

# Check the source
load(str_c(data_path, "/NEFSC_BTS_all_seasons_03032021.RData"))
survdat <- as.data.frame(survey$survdat)

# Stratum Area Key for which stratum correspond to larger regions we use
region_strata <-  c(
  # GOM
  as.character(13:23),
  #GB
  as.character(24:40),
  # SNE
  stringr::str_pad(as.character(1:12), width = 2, pad = "0", side = "left"),
  # MAB
  as.character(61:76))

# Dim check with filters
survdat %>%
  rename_all(tolower) %>%
  mutate(
    id  = stringr::str_c(cruise6, station, stratum),
    strat_num = stringr::str_sub(stratum, 2, 3)) %>%
  filter(
    season %in% c("SPRING", "FALL"),
    !is.na(biomass),
    !is.na(abundance),
    !is.na(length),
    !is.na(numlen),
    numlen > 0,
    strat_num %in% region_strata,
    stratum >= 01010,
    stratum <= 01760,
    stratum != 1310,
    stratum != 1320,
    stratum != 1330,
    stratum != 1350,
    stratum != 1410,
    stratum != 1420,
    stratum != 1490,
    svspp %not in% c(285:299, 305, 306, 307, 316, 323, 910:915, 955:961),
    svspp %not in% c(0, 978, 979, 980, 998)) %>%
  distinct(id, comname, catchsex, length, .keep_all = T) %>%
  dim()
