# Processing Functions
library(sizeSpectra)
library(tidyverse)


# Cleaning function
# loads survdat
# tidies up the names
# produces numlen_adj
# joins in species name key
# Assign region names to strata groups
tidy_nefsc_trawl <- function(data_path){
  
  
  ####  1. Import SURVDAT File  ####
  
  # Convenience change to make names lowercase
  load(str_c(data_path, "/NEFSC_BTS_all_seasons_03032021.RData"))
  survdat <- as.data.frame(survey$survdat)
  
  # Clean+lowercase names up for convenience
  trawldat <- janitor::clean_names(survdat)
  rm(survdat)
  
  
  #### 2.  Column Detection  ####
  # different versions of the data have come with different columns
  # generate these for consistency when comparing across
  
  ####_ a. Missing column flags  ####
  
  # Flags for missing columns that need to be merged in or built
  has_comname  <- "comname"     %in% names(trawldat)
  has_id_col   <- "id"          %in% names(trawldat)
  has_towdate  <- "est_towdate" %in% names(trawldat)
  has_month    <- "est_month"   %in% names(trawldat)
  
  # Flags for renaming or subsetting the data due to presence/absence of columns
  has_year      <- "est_year"      %in% names(trawldat)
  has_catchsex  <- "catchsex"      %in% names(trawldat)
  has_decdeg    <- "decdeg_beglat" %in% names(trawldat)
  has_avg_depth <- "avgdepth"      %in% names(trawldat)
  
  
  
  ####_ b. Missing comname  ####
  
  # Use SVSPP to get common names for species
  if(has_comname == FALSE){
    message("no comnames found, merging records in with spp_keys/sppclass.csv")
    
    # Load sppclass codes and common names
    spp_classes <- readr::read_csv(
      paste0(data_path, "sppclass.csv"),
      col_types = readr::cols())
    spp_classes <- janitor::clean_names(spp_classes)
    spp_classes <- dplyr::mutate(
      .data = spp_classes,
      comname  = stringr::str_to_lower(common_name),
      scientific_name = stringr::str_to_lower(scientific_name))
    spp_classes <- dplyr::distinct(spp_classes, svspp, comname, scientific_name)
    
    # Add the common names over and format for rest of build
    trawldat <- dplyr::mutate(trawldat, svspp = stringr::str_pad(svspp, 3, "left", "0"))
    trawldat <- dplyr::left_join(trawldat, spp_classes, by = "svspp")
    
  }
  
  
  ####_ c. Missing ID  ####
  if(has_id_col == FALSE) {
    message("creating station id from cruise-station-stratum fields")
    # Build ID column
    trawldat <- dplyr::mutate(
      .data = trawldat,
      cruise6 = stringr::str_pad(cruise6, 6, "left", "0"),
      station = stringr::str_pad(station, 3, "left", "0"),
      stratum = stringr::str_pad(stratum, 4, "left", "0"),
      id      = stringr::str_c(cruise6, station, stratum))}
  
  
  ####_ d. Field renaming  ####
  
  # Rename select columns for consistency
  if(has_year == FALSE)      {
    message("renaming year column to est_year")
    trawldat <- dplyr::rename(trawldat, est_year = year) }
  if(has_decdeg == FALSE) {
    message("renaming lat column to decdeg_beglat")
    trawldat <- dplyr::rename(trawldat, decdeg_beglat = lat) }
  if(has_decdeg == FALSE) {
    message("renaming lon column to decdeg_beglon")
    trawldat <- dplyr::rename(trawldat, decdeg_beglon = lon) }
  if(has_avg_depth == FALSE)      {
    message("renaming depth column to avgdepth")
    trawldat <- dplyr::rename(trawldat, avgdepth = depth) }
  
  
  
  ####____ d. build date structure for quick grab of date components
  if(has_towdate == TRUE) {
    message("building month/day columns from est_towdate")
    trawldat <- dplyr::mutate(
      .data = trawldat,
      est_month = stringr::str_sub(est_towdate, 6,7),
      est_month = as.numeric(est_month),
      est_day   = stringr::str_sub(est_towdate, -2, -1),
      est_day   = as.numeric(est_day), .before = season)}
  
  
  
  #### 4. Column Changes  ####
  
  # a. Text Formatting
  trawldat <- dplyr::mutate(
    .data = trawldat,
    comname = tolower(comname),
    id      = format(id, scientific = FALSE),
    svspp   = as.character(svspp),
    svspp   = stringr::str_pad(svspp, 3, "left", "0"),
    season  = stringr::str_to_title(season),
    
    # Format Stratum number,
    # exclude leading and trailing codes for inshore/offshore,
    # used for matching to stratum areas
    strat_num = stringr::str_sub(stratum, 2, 3))
  
  
  # b. Rename to make units more clear
  trawldat <- dplyr::rename(
    .data = trawldat,
    biomass_kg = biomass,
    length_cm  = length)
  
  # c. Replace 0's that must be greater than 0
  trawldat <- dplyr::mutate(
    .data = trawldat,
    biomass_kg = ifelse(biomass_kg == 0 & abundance > 0, 0.0001, biomass_kg),
    abundance  = ifelse(abundance == 0 & biomass_kg > 0, 1, abundance))
  
  
  #### 5. Row Filtering  ####
  
  # Things filtered:
  # 1. Strata
  # 2. Seasons
  # 3. Year limits
  # 4. Vessels
  # 5. Species Exclusion
  
  # Eliminate Canadian Strata and Strata No longer in Use
  trawldat <- dplyr::filter(
    .data = trawldat,
    stratum >= 01010,
    stratum <= 01760,
    stratum != 1310,
    stratum != 1320,
    stratum != 1330,
    stratum != 1350,
    stratum != 1410,
    stratum != 1420,
    stratum != 1490)
  
  
  # Drop NA Biomass and Abundance Records
  trawldat <- dplyr::filter(
    .data = trawldat,
    !is.na(biomass_kg),
    !is.na(abundance))
  
  # Exclude the Skrimps
  trawldat <-  dplyr::filter(
    .data = trawldat,
    svspp %not in% c(285:299, 305, 306, 307, 316, 323, 910:915, 955:961))
  
  # Exclude the unidentified fish
  trawldat <- dplyr::filter(trawldat, svspp %not in% c(0, 978, 979, 980, 998))
  
  
  #### 6. Spatial Filtering - Stratum  ####
  
  # This section merges stratum area info in
  # And drops stratum that are not sampled or in Canada
  # these are used to relate catch/effort to physical areas in km squared
  
  # Stratum Area Key for which stratum correspond to larger regions we use
  strata_key <- list(
    "Georges Bank"          = as.character(13:23),
    "Gulf of Maine"         = as.character(24:40),
    "Southern New England"  = stringr::str_pad(
      as.character(1:12), width = 2, pad = "0", side = "left"),
    "Mid-Atlantic Bight"    = as.character(61:76))
  
  
  # Add the labels to the data
  trawldat <- dplyr::mutate(
    trawldat,
    survey_area =  dplyr::case_when(
      strat_num %in% strata_key$`Georges Bank`         ~ "GB",
      strat_num %in% strata_key$`Gulf of Maine`        ~ "GoM",
      strat_num %in% strata_key$`Southern New England` ~ "SNE",
      strat_num %in% strata_key$`Mid-Atlantic Bight`   ~ "MAB",
      TRUE                                             ~ "stratum not in key"))
  
  
  # Use strata_select to pull the strata we want individually
  strata_select <- c(strata_key$`Georges Bank`,
                     strata_key$`Gulf of Maine`,
                     strata_key$`Southern New England`,
                     strata_key$`Mid-Atlantic Bight`)
  
  
  # Filtering areas using strata_select
  trawldat <- dplyr::filter(trawldat, strat_num %in% strata_select)
  trawldat <- dplyr::mutate(trawldat, stratum = as.character(stratum))
  
  
  #### 7. Adjusting NumLength  ####
  
  # NOTE:
  # numlen is not adjusted to correct for the change in survey vessels and gear
  # these values consequently do not equal abundance, nor biomass which are adjusted
  
  # Because of this and also some instances of bad data,
  # there are cases of more/less measured than initially tallied* in abundance
  # this section ensures that numlen totals out to be the same as abundance
  
  
  # If catchsex is not a column then total abundance is assumed pooled
  if(has_catchsex == TRUE){
    abundance_groups <- c("id", "comname", "catchsex", "abundance")
  } else {
    message("catchsex column not found, ignoring sex for numlen adjustments")
    abundance_groups <- c("id", "comname", "abundance")}
  
  
  # Get the abundance value for each sex
  # arrived at by summing across each length
  abundance_check <- dplyr::group_by(trawldat, !!!rlang::syms(abundance_groups))
  abundance_check <- dplyr::summarise(
    .data = abundance_check,
    abund_actual = sum(numlen),
    n_len_class  = dplyr::n_distinct(length_cm),
    .groups      = "drop")
  
  
  # Get the ratio between the original abundance column
  # and the sum of numlen we just grabbed
  conv_factor <- dplyr::distinct(trawldat, !!!rlang::syms(abundance_groups), length_cm)
  conv_factor <- dplyr::inner_join(conv_factor, abundance_check, by = abundance_groups)
  conv_factor <- dplyr::mutate(conv_factor, convers = abundance / abund_actual)
  
  # Merge back and convert the numlen field
  # original numlen * conversion factor = numlength adjusted
  survdat_processed <- dplyr::left_join(trawldat, conv_factor, by = c(abundance_groups, "length_cm"))
  survdat_processed <- dplyr::mutate(survdat_processed, numlen_adj = numlen * convers, .after = numlen)
  survdat_processed <- dplyr::select(survdat_processed, -c(abund_actual, convers))
  
  # remove conversion factors from environment
  rm(abundance_check, conv_factor, strata_key, strata_select)
  
  
  
  #### 8. Distinct Station & Species Length Info   ####
  
  # For each station we need unique combinations of
  # station_id, species, catchsex, length_cm, adjusted_numlen
  # to capture what and how many of each length fish is caught
  
  # Record of unique station catches:
  # One row for every species * sex * length_cm, combination in the data
  trawl_lens <- dplyr::filter(
    .data = survdat_processed,
    is.na(length_cm) == FALSE,
    is.na(numlen) == FALSE,
    numlen_adj > 0)
  
  
  # Do we want to just keep all the station info here as well?
  # question to answer is whether any other columns repeat,
  # or if these are the only ones
  trawl_clean <- dplyr::distinct(
    .data = trawl_lens,
    id, 
    svspp, 
    comname, 
    catchsex, 
    abundance, 
    n_len_class,
    length_cm, 
    numlen, 
    numlen_adj, 
    biomass_kg, .keep_all = TRUE)
  
  
  
  
  # Return the dataframe
  # Contains 1 Row for each length class of every species caught
  return(trawl_clean)
  
}





####___________________________####

add_wigley_lw <- function(tidy_trawl, data_path){
  
  
  #### 1. Match Species to LW Coefficients  ####

  ####  Resource Paths
  lw_key_path <- paste0(data_path, "wigley_06_lwreg.csv")
  
  # This table is a combined table of wigley and fishbase L-W coefficients
  wigley_coefs <- readr::read_csv(lw_key_path, col_types = readr::cols()) %>% 
    janitor::clean_names()
  wigley_coefs <- dplyr::mutate(
    wigley_coefs,
    svspp = stringr::str_pad(svspp, 3, "left", "0"),
    season = stringr::str_to_title(season),
    common_name = tolower(common_name),
    scientific_name = tolower(scientific_name))
  
  # Pull important columns, rename if necessary
  wigley_coefs <- dplyr::select(
    wigley_coefs,
    season, svspp, comname = common_name, scientific_name, catchsex, b, ln_a = lna)
  
  
  
  # Mismatched svspp code(s) - Fix scup in trawl data
  wigley_lookup <- function(x){
    unique(wigley_coefs$svspp[which(wigley_coefs$comname == x)])}
  
  tidy_trawl <- dplyr::mutate(
    .data = tidy_trawl,
    svspp = dplyr::if_else(comname == "scup", wigley_lookup("scup"), svspp))
  
  
  
  
  # Merge on comname, season, and catchsex
  # Join just by svspp to account for common name changes or misspellings
  trawl_weights <- dplyr::select(tidy_trawl, -comname) %>% 
     dplyr::inner_join(wigley_coefs)
  trawl_weights <- dplyr::arrange(trawl_weights, est_year, season)
  
  # Estimate Weights from Lengths
  trawl_weights <- dplyr::mutate(
    trawl_weights,
    b             = as.numeric(b),
    a             = exp(ln_a), 
    llen          = log(length_cm),
    ind_log_wt    = ln_a + (b * llen),
    ind_weight_kg = exp(ind_log_wt),                # weight of an individual in size class
    sum_weight_kg = ind_weight_kg * numlen_adj)     # Individual weight * adjusted numlen
  
  trawl_weights <- tidyr::drop_na(trawl_weights, ind_weight_kg)
  trawl_weights <- dplyr::select(trawl_weights, -ind_log_wt, -llen)
  
  ####  2. Use Coefficients to Estimate length-based biomass  ####
  
  # calculate total biomass again using weights from key
  # make a key for the length weight coefficient sources
  survdat_weights <- dplyr::arrange(trawl_weights, est_year, season, comname, length_cm)
  survdat_weights <- dplyr::mutate(survdat_weights, lw_group = stringr::str_c(comname, season, catchsex))
  
  # Return the data
  return(survdat_weights)
  
  
}




#___________________________####


# ----- Make it a function  ------------


# Code to do mle for one group of data:



#' @title {MLE-Bins Size Spectra Estimation}
#'
#' @param ss_input Dataframe containing a column of abundance and a column of sizes
#' @param grouping_vars string identifiers of columns to group_by prior to analysis
#' @param abundance_vals string indicating column of abundances
#' @param size_vals string indicating column with individual length/weight data
#' @param isd_xmin lower limit for size distribution fitting
#' @param isd_xmax upper limit for size distribution fitting
#' @param global_min T/F Enforce a minimum size across all groups
#' @param global_max T/F Enforce a maximum size across all groups
#'
#' @return
#' @export
#'
#' @examples
group_binspecies_mle <-  function(
    ss_input, 
    grouping_vars, 
    abundance_vals = "numlen_adj",
    size_vals = "length_cm",
    isd_xmin = NULL,
    isd_xmax = NULL,
    global_min = TRUE,
    global_max = TRUE,
    bin_width = 1){
  
  # 1. toggle which abundance/weight columns to use with switch
  .abund_col  <- sym(abundance_vals)
  .size_col   <- sym(size_vals)
  .group_cols <- grouping_vars
  .agg_cols   <- c(.group_cols, "comname", size_vals)
  
  
  # 2. Select only the columns we need:
  # Aggregate numbers by size 
  # within the groups we are measuring:
  ss_input_summs <- ss_input %>% 
    group_by(!!!syms(.agg_cols)) %>% 
    summarise(
      Number = sum(!!.abund_col),
      .groups = "drop")
  
  
  # Create a group column that we can map() through
  # Drop columns we don't need, rename to match edwards code 
  # edwards calls this df databinforlike in eightMethods()
  mle_input <- ss_input_summs  %>% 
    unite(col = "group_var", {{.group_cols}}, sep = "-", remove = FALSE, na.rm = FALSE)  %>% 
    select(
      group_var, 
      Number, 
      wmin = !!.size_col)  %>% 
    # Set up the min/max for the bins
    mutate(wmax = wmin + bin_width)
  
  #------ Optional - Set Global Constants
  
  # Power-law limits:
  # set left bound & right bounds, grams
  if(global_min == TRUE){
    if(is.null(isd_xmin)){ isd_xmin <- min(mle_input$wmin, na.rm = T)}
  }
  if(global_max == TRUE){
    if(is.null(isd_xmax)){ isd_xmax <- max(mle_input$wmax, na.rm = T)}
  }
  
  
  #------ Loop through Groups 
  group_results_df <- mle_input %>% 
    split(.$group_var) %>% 
    map_df(
      function(ss_input_i){
        
        # 3. Tune subgroup the power-law limits:
        # set left bound & right bounds, grams/cm
        min_i <- min(ss_input_i$wmin, na.rm = T)
        max_i <- max(ss_input_i$wmax, na.rm = T)
        if(global_min == FALSE){
          if(is.null(isd_xmin)){ isd_xmin <- min_i}
        }
        if(global_max == FALSE){
          if(is.null(isd_xmax)){ isd_xmax <- max_i}
        }
        
        # Filter the range of sizes we want to include
        ss_input_i <- ss_input_i %>% 
          filter(wmin >= isd_xmin,
                 wmax <= isd_xmax)
        
        # Total individuals in subgroup
        n_i <- sum(ceiling(ss_input_i$Number) )
        
        # Do the exponent estimation for group i
        group_est <- calcLike(
          negLL.fn          = negLL.PLB.bins.species,
          p                 = -1.9,
          suppress.warnings = TRUE,
          dataBinForLike    = ss_input_i,
          n                 = n_i,
          xmin              = isd_xmin,
          xmax              = isd_xmax, 
          vecDiff = 2)
        
        # Put it into a dataframe to rejoin neatly
        mle_group_results <- data.frame(
          xmin_fit = isd_xmin,
          xmax_fit = isd_xmax,
          xmin_actual = min_i,
          xmax_actual = max_i,
          n = n_i,
          b = group_est$MLE,
          confMin = group_est$conf[1],
          confMax = group_est$conf[2])
        
        
        # Process C and standard error
        mle_group_results <- mle_group_results %>% 
          mutate(
            stdErr = (abs(confMin - b) + abs(confMax - b)) / (2 * 1.96),
            C = (b != -1 ) * (b + 1) / ( xmax_fit^(b + 1) - xmin_fit^(b + 1) ) + (b == -1) * 1 / ( log(xmax_fit) - log(xmin_fit)))
        return(mle_group_results)}, 
      
      # Name the column with group ID
      .id = "group_var") %>% 
    
    # Decompose the group ID back into original columns
    separate(
      group_var, 
      sep = "-", 
      into = grouping_vars, 
      remove = F)
  
  
  # Spit it out
  return(group_results_df)
  
}




#' @title {MLE-Bins Size Spectra Estimation}
#'
#' @param ss_input Dataframe containing a column of abundance and a column of sizes
#' @param grouping_vars string identifiers of columns to group_by prior to analysis
#' @param abundance_vals string indicating column of abundances
#' @param size_vals string indicating column with individual length/weight data
#' @param isd_xmin lower limit for size distribution fitting
#' @param isd_xmax upper limit for size distribution fitting
#' @param global_min T/F Enforce a minimum size across all groups
#' @param global_max T/F Enforce a maximum size across all groups
#'
#' @return
#' @export
#'
#' @examples
group_binspecies_bodymass_spectra <-  function(
    ss_input, 
    grouping_vars, 
    abundance_vals = "numlen_adj",
    size_vals = "length_cm",
    isd_xmin = NULL,
    isd_xmax = NULL,
    global_min = TRUE,
    global_max = TRUE,
    bin_width = 1){
  
  # 1. toggle which abundance/weight columns to use with switch
  .abund_col  <- sym(abundance_vals)
  .size_col   <- sym(size_vals)
  .group_cols <- grouping_vars
  .agg_cols   <- c(.group_cols, "comname", size_vals)
  
  
  # 2. Select only the columns we need:
  # Aggregate numbers by size 
  # within the groups we are measuring:
  ss_input_summs <- ss_input %>% 
    group_by(!!!syms(.agg_cols)) %>% 
    summarise(
      Number = sum(!!.abund_col),
      .groups = "drop")
  
  
  # Create a group column that we can map() through
  # Drop columns we don't need, rename to match edwards code 
  # edwards calls this df databinforlike in eightMethods()
  mle_input <- ss_input_summs  %>% 
    unite(col = "group_var", {{.group_cols}}, sep = "-", remove = FALSE, na.rm = FALSE)  %>% 
    select(
      group_var, 
      Number, 
      wmin = !!.size_col)  %>% 
    # Set up the min/max for the bins
    mutate(wmax = wmin + bin_width)
  
  #------ Optional - Set Global Constants
  
  # Power-law limits:
  # set left bound & right bounds, grams
  if(global_min == TRUE){
    if(is.null(isd_xmin)){ isd_xmin <- min(mle_input$wmin, na.rm = T)}
  }
  if(global_max == TRUE){
    if(is.null(isd_xmax)){ isd_xmax <- max(mle_input$wmax, na.rm = T)}
  }
  
  
  #------ Loop through Groups 
  group_results_df <- mle_input %>% 
    split(.$group_var) %>% 
    map_df(
      function(ss_input_i){
        
        # 3. Tune subgroup the power-law limits:
        # set left bound & right bounds, grams/cm
        min_i <- min(ss_input_i$wmin, na.rm = T)
        max_i <- max(ss_input_i$wmax, na.rm = T)
        if(global_min == FALSE){
          if(is.null(isd_xmin)){ isd_xmin <- min_i}
        }
        if(global_max == FALSE){
          if(is.null(isd_xmax)){ isd_xmax <- max_i}
        }
        
        # Filter the range of sizes we want to include
        ss_input_i <- ss_input_i %>% 
          filter(wmin >= isd_xmin,
                 wmax <= isd_xmax)
        
        # Total individuals in subgroup
        n_i <- sum(ceiling(ss_input_i$Number) )
        
        # Do the exponent estimation for group i
        group_est <- calcLike(
          negLL.fn          = negLL.PLB.bins.species,
          p                 = -1.9,
          suppress.warnings = TRUE,
          dataBinForLike    = ss_input_i,
          n                 = n_i,
          xmin              = isd_xmin,
          xmax              = isd_xmax, 
          vecDiff = 2)
        
        # Put it into a dataframe to rejoin neatly
        mle_group_results <- data.frame(
          xmin_fit = isd_xmin,
          xmax_fit = isd_xmax,
          xmin_actual = min_i,
          xmax_actual = max_i,
          n = n_i,
          b = group_est$MLE,
          confMin = group_est$conf[1],
          confMax = group_est$conf[2])
        
        
        # Process C and standard error
        mle_group_results <- mle_group_results %>% 
          mutate(
            stdErr = (abs(confMin - b) + abs(confMax - b)) / (2 * 1.96),
            C = (b != -1 ) * (b + 1) / ( xmax_fit^(b + 1) - xmin_fit^(b + 1) ) + (b == -1) * 1 / ( log(xmax_fit) - log(xmin_fit)))
        return(mle_group_results)}, 
      
      # Name the column with group ID
      .id = "group_var") %>% 
    
    # Decompose the group ID back into original columns
    separate(
      group_var, 
      sep = "-", 
      into = grouping_vars, 
      remove = F)
  
  
  # Spit it out
  return(group_results_df)
  
}




