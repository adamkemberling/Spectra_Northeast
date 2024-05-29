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
  
  # # c. Replace 0's that must be greater than 0
  # trawldat <- dplyr::mutate(
  #   .data = trawldat,
  #   biomass_kg = ifelse(biomass_kg == 0 & abundance > 0, 0.0001, biomass_kg),
  #   abundance  = ifelse(abundance == 0 & biomass_kg > 0, 1, abundance))
  
  
  # c. EDIT: Remove these^ because downstream problems
  no_biomass <- trawldat$biomass_kg == 0 & trawldat$abundance > 0
  no_abundance <- trawldat$abundance == 0 & trawldat$biomass_kg > 0
  trawldat %>% 
    filter(!no_biomass,
           !no_abundance)
  
  
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


# ----- Processing sizeSpectra for Groups  ------------


# Tweaking the likelihood function:
# It seems to be solving to -1 in a strange way, not sure if its the fail state of 
# too small a vecdiff or, NA's not filtered

# ##' Calculate the negative log-likelihood of `b` for the PLB model, given
# ##'  species-specific binned data (MLEbins method)
# ##'
# ##' Calculate the negative log-likelihood of *b* for the PLB model,
# ##'  given binned data where the bins can be different for each species, namely
# ##'  the MLEbins method derived as equations (S.18) and (S.26) in MEPS paper.
# ##'  Returns the negative log-likelihood.
# ##'  Will be called by `nlm()` or similar, but `xmin` and `xmax` will just be estimated
# ##'  as the min of lowest bin and max of the largest bin (i.e. their MLEs),
# ##'  no need to do numerically. See Supplementary Material of MEPS paper for derivation, and
# ##'  the vignettes for example use.
# ##'
# ##' @param b value of `b` for which to calculate the negative log-likelihood
# ##' @param dataBinForLike table data frame (tbl_df) where each row is the count in a bin
# ##'  of a species, and columns (and corresponding mathematical notation in MEPS
# ##'   Supplementary Material) are:
# ##'  * `SpecCode`: code for each species, `s`
# ##'  * `wmin`: lower bound of the bin, `w_\{sj\}` where `j` is the bin number
# ##'  * `wmax`: upper bound of the bin, `w_\{s, j+1\}`
# ##'  * `Number`: count in that bin for that species, `d_\{sj\}`
# ##'  For each species the first and last bins must be non-empty, i.e.
# ##'   `w_\{s1\}, w_\{s,J_s +1\} > 0`.
# ##' @param n total number of counts `n = \sum_\{sj\} d_\{sj\}` over all `s` and `j`
# ##' @param xmin maximum likelihood estimate for `xmin`, `xmin = min_\{sj\}
# ##'   w_\{s, 1\}`
# ##' @param xmax maximum likelihood estimate for `xmax`, `xmax = max_\{sj\}
# ##'   w_\{s, J_s+1\}`
# ##' @return  negative log-likelihood of the parameters given the data
# ##' @author Andrew Edwards
# ##' @export
# negLL.PLB.bins.species.new = function(
#     b, dataBinForLike, n, xmin, xmax){
#   # Would be useful to put into a pre-processing function:
#   #    if(xmin <= 0 | xmin >= xmax | length(d) != J | length(w) != J+1 |
#   #       d[1] == 0 | d[J] == 0 | min(d) < 0)
#   #       stop("Parameters out of bounds in negLL.PLB.bins.species")
#   # if(b != -1)
#   # {  
#   #   # From MEPS equation (S.18), first calculate each component in the
#   #   # summations and then sum them:
#   #   temp2 = dplyr::mutate(
#   #     dataBinForLike,
#   #     comp2 = Number * log( abs( wmax^(b+1) - wmin^(b+1) ) ) )
#   #   comp2Sum = sum(temp2$comp2)
#   #   
#   #   logLL = - n * log( abs( xmax^(b+1) - xmin^(b+1) ) ) + comp2Sum
#   #   neglogLL = - logLL      # Negative log-likelihood
#   # } else
#   # {
#     # Not fully tested, 
#     # but should work:
#     temp2 = dplyr::mutate(
#       dataBinForLike,
#       comp2 = Number * log( log(wmax) - log(wmin) ) )
#     comp2Sum = sum(temp2$comp2)
#     
#     logLL = - n * log( abs( xmax^(b+1) - xmin^(b+1) ) ) + comp2Sum
#     neglogLL = - logLL      # Negative log-likelihood
#   
#   return(neglogLL)
# }





#' @title {MLE-Bins Size Spectra Estimation - Length or Weight}
#'
#' @param ss_input Dataframe containing a column of abundance and a column of sizes
#' @param grouping_vars string identifiers of columns to group_by prior to analysis
#' @param abundance_vals string indicating column of abundances
#' @param size_vals string indicating column with individual length/weight data
#' @param use_weight FALSE
#' @param isd_xmin lower limit for size distribution fitting
#' @param isd_xmax upper limit for size distribution fitting
#' @param global_min T/F Enforce a minimum size across all groups
#' @param global_max T/F Enforce a maximum size across all groups
#' @param vdiff value defining  range over which to test the negative log-likelihood
#'   to construct the confidence interval; range is `MLE` \eqn{\pm} `vecDiff`. Default is 0.5 and a symmetric
#'   range is tested for fitting size spectra, since for movement data
#'   sets in Table 2 of Edwards (2011; 92(6):1247-1257) the intervals were
#'   symmetric, so symmetric seems a good start.
#'
#' @return
#' @export
#'
#' @examples
group_binspecies_spectra <-  function(
    ss_input, 
    grouping_vars, 
    abundance_vals = "numlen_adj",
    length_vals = "length_cm",
    use_weight = FALSE,
    isd_xmin = NULL,
    isd_xmax = NULL,
    global_min = TRUE,
    global_max = TRUE,
    bin_width = 1,
    vdiff = 1){
  
  
  # 1. toggle which abundance/weight columns to use with switch
  .abund_col  <- sym(abundance_vals)
  .size_col   <- sym(length_vals)
  .group_cols <- grouping_vars
  .agg_cols   <- c(.group_cols, "comname", length_vals)
  
  
  
  # Add min/max weight for size bins if using weight
  if(use_weight == TRUE){
    # Make sure they come through the aggregation step
    .agg_cols <- c(.agg_cols, "wmin_g", "wmax_g")
  }
          
  
  # 2. Aggregate numbers by size 
  # within the groups we are measuring:
  ss_input_summs <- ss_input %>% 
    group_by(!!!syms(.agg_cols)) %>% 
    summarise(
      Number = sum(!!.abund_col),
      .groups = "drop")
  
  
  # Create a group column that we can map() through
  # edwards calls this df databinforlike in eightMethods()
  mle_input <- ss_input_summs  %>% 
    unite(
      col = "group_var", 
      {{.group_cols}}, 
      sep = "-", 
      remove = TRUE, 
      na.rm = TRUE)  
  
  
  # 2. Rename important columns to use either length or weight bins
  # the calclike function expects: 
  # "wmin" and "wmax" for the bin min/max
  # SpecCode for species ID
  # Number for number of individuals within the bins
  
  # For lengths
  if(use_weight == FALSE){
    
    # Set up the min/max for the bins
    mle_input <- mle_input %>% 
      rename(wmin = !!.size_col)  %>% 
      mutate(wmax = wmin + bin_width) %>% 
      drop_na(wmin, wmax, Number)
  }
  
  # For weights
  if(use_weight == TRUE){
    
    # Set up the min/max for the bins
    mle_input <- mle_input %>% 
      rename(wmin = wmin_g,
             wmax = wmax_g) %>% 
      drop_na(wmin, wmax, Number)
  }

  
  #------ Optional - Set Global Constants
  
  # Power-law limits:
  # set left bound & right bounds, grams
  if(global_min == TRUE){
    if(is.null(isd_xmin)){ isd_xmin <- min(mle_input$wmin, na.rm = T)}
  }
  if(global_max == TRUE){
    if(is.null(isd_xmax)){ isd_xmax <- max(mle_input$wmax, na.rm = T)}
  }
  
  #------ End optional
  
  
  # Subset to just the columns needed to speed it up a touch
  mle_input <- mle_input %>% 
    dplyr::select(
      group_var,
      SpecCode = comname,
      wmin,
      wmax,
      Number)
  
  

  
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
                 wmin <= isd_xmax) %>% 
         mutate(Number = ceiling(Number))
        
        # Total individuals in subgroup
        #n_i <- sum(ceiling(ss_input_i$Number) )
        n_i <- sum(ss_input_i$Number)
        
        # return(
        #   list(
        #     n = n_i,
        #     xmin = isd_xmin,
        #     xmax = isd_xmax,
        #     vdiff = vdiff
        #   )
        # )
        
        # Do the exponent estimation for group i
        group_est <- calcLike(
          #negLL.fn          = negLL.PLB.bins.species.new,
          negLL.fn          = negLL.PLB.bins.species,
          suppress.warnings = TRUE,
          dataBinForLike    = ss_input_i,
          p                 = 1.9,
          n                 = n_i,
          xmin              = isd_xmin,
          xmax              = isd_xmax, 
          vecDiff           = vdiff)
        
        
        # Put it into a dataframe to rejoin neatly
        mle_group_results <- data.frame(
          xmin_fit = isd_xmin,
          xmax_fit = isd_xmax,
          xmin_actual = min_i,
          xmax_actual = max_i,
          n = n_i,
          b = group_est$MLE,
          confMin = group_est$conf[1],
          confMax = group_est$conf[2],
          spectra_type = ifelse(use_weight, "bodymass", "length"))
        
        
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
      remove = T)
  
  
  # Spit it out
  return(group_results_df)
  
}



###___________________________####


#### Processing Mean+Median Sizes  ####



# Get weighted mean lengths and weights using 
# tow-level and total stratified abundances
# Note: uses numlen because individuals were ID'd and measured
group_size_metrics <- function(
    size_data, 
    .group_cols = "Year", 
    .abund_col = "numlen_adj",
    .length_col = "length_cm",
    has_weights = FALSE,
    .weight_col = "ind_weight_kg"){
  
  
  # 1. Build group_level from desired group columns
  group_size_data <- size_data %>% 
    drop_na({{.abund_col}}, {{.length_col}}) %>% 
    unite(
      col = "group_var", 
      {{.group_cols}}, 
      sep = "-", 
      remove = TRUE, 
      na.rm = TRUE)
  
  
  # Run Min/Max/Avg. Size for the group
  group_results <- group_size_data %>% 
    split(.$group_var) %>% 
    imap_dfr(function(group_data, group_name){
      
      # Length (measured for all)
      mean_len    <- weighted.mean(
        group_data[, .length_col], 
        group_data[, .abund_col], 
        na.rm = T)
      med_len     <- matrixStats::weightedMedian(
        group_data[, .length_col], 
        group_data[, .abund_col], 
        na.rm = T)
      min_len     <- min(group_data[, .length_col], na.rm = T)
      max_len     <- max(group_data[, .length_col], na.rm = T)
      
      
      
      # Weights (estimated from length w/ wigley data)
      if(has_weights == TRUE){
        group_data <- drop_na(group_data, {{.weight_col}})
        
        mean_weight <- weighted.mean(
          group_data[,.weight_col], 
          group_data[, .abund_col], 
          na.rm = T)
        med_weight  <- matrixStats::weightedMedian(
          group_data[, .weight_col], 
          group_data[, .abund_col], 
          na.rm = T)
        min_weight  <- min(group_data[, .weight_col], na.rm = T)
        max_weight  <- max(group_data[, .weight_col], na.rm = T)
      }
      
      # Abundance totals
      total_abund <- sum(group_data[, .abund_col], na.rm = T)
      
      # Total number of species
      num_species <- group_data %>% 
        filter(is.na(comname) == FALSE) %>% 
        distinct(comname) %>% 
        nrow()
      
      # Put in table
      table_out <- data.frame(
        "group_var"        = group_name,
        "n_species"        = num_species,  
        "numlen_adj_total" = total_abund,
        "mean_len_cm"      = mean_len,
        "med_len_cm"       = med_len,
        "min_len_cm"       = min_len,
        "max_len_cm"       = max_len)
      
      # Add weight information
      if(has_weights == TRUE){
        table_out <- bind_cols(
          table_out,
          data.frame(
            "mean_wt_kg"       = mean_weight,
            "med_wt_kg"        = med_weight,
            "min_wt_kg"        = min_weight,
            "max_wt_kg"        = max_weight)) %>% 
          
          # Replace Inf with NA
          mutate(across(
            .cols = numlen_adj_total:max_wt_kg, 
            .fns = ~ifelse(is.infinite(abs(.x)), NA, .x)))
        
      }
        
        
        
      
      # return the table
      return(table_out)
    }) %>% 
    # Decompose the group ID back into original columns
    separate(
      group_var, 
      sep = "-", 
      into = .group_cols, 
      remove = T)
  
  
  
  # Return the results
  return(group_results)
  
  
}







# rstan compiler options: Only helpful/needed for isdbayes
# https://github.com/stan-dev/rstan/wiki/Configuring-C---Toolchain-for-Mac
# dotR <- file.path(Sys.getenv("HOME"), ".R")
# if (!file.exists(dotR)) dir.create(dotR)
# M <- file.path(dotR, "Makevars")
# if (!file.exists(M)) file.create(M)
# arch <- ifelse(R.version$arch == "aarch64", "arm64", "x86_64")
# cat(paste("\nCXX17FLAGS += -O3 -mtune=native -arch", arch, "-ftemplate-depth-256"),
#     file = M, sep = "\n", append = FALSE)
