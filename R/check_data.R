# Copyright 2024 Fisheries Queensland

# This file is part of DDUST.
# DDUST is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DDUST is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with DDUST. If not, see <https://www.gnu.org/licenses/>.

#'  Check input data for DDUST model
#'
#' @param data list of data
#' @param silent TRUE show warnings and messages
#' 
#' @description
#' 'check_data' returns an amended list of model input data. See details.
#' 
#' @details
#' The minimum requirements for data include
#' - proportion_spawning: a vector of length 12 indicating the proportion of spawning that occurs each month
#' - weight_at_recruitment: a vector of length 2 indicating the weight of a pre-recruit and the weight of a recruit
#' - weight_inf: the average maximum weight
#' - last_year_catch: the last year of catch (int)
#' - first_year_catch: the first year of catch (int)
#' - first_year_rec_devs: the first year of recruitment deviations (int)
#' - ctch: a vector of catch
#' - cpue: an array of standardised catch per unit effort with dimensions with rows for fleets
#' - cpue_sd: an array of standardised catch per unit effort standard deviation with dimensions with rows for fleets
#' - month_sequence: a vector of months with 1 indicating the warmest month
#' - calculate_rho: an integer indicating if rho is calculated inside the model using weight_at_recruitment and weight_inf
#' - use_recruit_penalty: an integer indicating if the recruitment penalty is used
#' - prior_mean_mu: a value for the mean of the mu normal prior
#' - prior_mean_k: a value for the mean of the k normal prior
#' - prior_mean_xi: a value for the mean of the xi normal prior
#' 
#' The following data will be amended with default values:
#' - rec_dev_type = 'off': treatment of recruitment deviations ('off', 'barycentric', 'redundant')
#' - Number_months_per_timestep = 1: number of months per timestep (1, 2, 3, 4, 6 or 12)
#' - use_cpue_sd = 0: an integer to indicate if cpue_sd will be used
#' - num_years_warmup = 100: an integer to indicate number of years to run equilibirum warm up
#' - recruit_penalty_strength = 0.001: a value to scale the strength of the recruitment penalty
#' - recruit_penalty_exponent = 2: a value to scale the strength of the recruitment penalty
#' - absolute_biomass = (0,0,...,0): a vector of observed absolute biomass data
#' - absolute_biomass_sd  = (1,1,...,1): a vector of observed absolute biomass standard deviation data
#' - minimum_annual_harvest = 0.05: a ratio indicating the harvest rate that induces the harvest penalty
#' - do_projections = 0: an integer to indicate if projections are run
#' - projection_years = 50: the number of years to run projections (int)
#' - projection_harvest = (0,0,...,0): a vector of harvest used in harvest projections
#' - target_relative_biomass  = 0.6: the relative biomass used in target projections
#' - F_initial = 0.05: the initial fishing mortality used for MSY and FMSY calculations
#' - prior_sd_mu = 1: a value for the standard deviation of the mu normal prior
#' - prior_sd_k = 1: a value for the standard deviation of the k normal prior
#' - prior_sd_xi = 1: a value for the standard deviation of the xi normal prior
#' 
#' @return Amended list of model input data
#' @export
#'
#' @examples
#' data <- check_data(data)
check_data <- function(data, silent = FALSE){
  # Data provided to R function (no default): ----
  # DATA_VECTOR(proportion_spawning);
  if (is.null(data$proportion_spawning)) {
    stop("No 'data$proportion_spawning' value defined. \n
           Please provide a vector of length 12 indicating the proportion of spawning that occurs each month. \n
           Vector will aggregate according to data$Number_months_per_timestep.")
  }
    if (!length(data$proportion_spawning) == 12) {
    stop("Incorrect length of 'data$proportion_spawning' vector defined. \n
           Please provide a vector of length 12 indicating the proportion of spawning that occurs each month. \n
           Vector will aggregate according to data$Number_months_per_timestep.")
  }
  # DATA_VECTOR(weight_at_recruitment);
  if (is.null(data$weight_at_recruitment)) {
    stop("No 'data$weight_at_recruitment' value defined. \n
             Please provide a vector of length 2 indicating the weight of a pre-recruit and the weight of a recruit.")
  }
  if (!length(data$weight_at_recruitment) == 2) {
    stop("Incorrect length of 'data$weight_at_recruitment' vector defined. \n
             Please provide a vector of length 2 indicating the weight of a pre-recruit and the weight of a recruit.")
  }
  # DATA_SCALAR(weight_inf);
  if (is.null(data$weight_inf)) {
    stop("No 'data$weight_inf' value defined. \n
             Please provide the average maximum weight.")
  }
  # DATA_INTEGER(last_year_catch);
  if (is.null(data$last_year_catch)) {
    stop("No 'data$last_year_catch' value defined. \n
             Please provide the last year of catch.")
  }
  # DATA_INTEGER(first_year_catch);
  if (is.null(data$first_year_catch)) {
    stop("No 'data$first_year_catch' value defined. \n
             Please provide the first year of catch.")
  }
  # DATA_INTEGER(first_year_rec_devs);
  if (is.null(data$first_year_rec_devs)) {
    stop("No 'data$first_year_rec_devs' value defined. \n
             Please provide the first year of recruitment deviations.")
  }
  # DATA_VECTOR(ctch);
  if (is.null(data$ctch)) {
    stop("No 'data$ctch' value defined. \n
             Please provide a vector of catch.")
  }
  # DATA_ARRAY(cpue);
  if (is.null(data$cpue)) {
    stop("No 'data$cpue' value defined. \n
             Please provide an array of standardised catch per unit effort.")
  }
  # DATA_ARRAY(cpue_sd);
  if (is.null(data$cpue_sd)) {
    stop("No 'data$cpue_sd' value defined. \n
             Please provide an array of standardised catch per unit effort standard deviation.")
  }
  if (!any(dim(data$cpue)==dim(data$cpue_sd))){
    warning("Dimension of cpue and cpue_sd should match.")
  }
  # DATA_VECTOR(month_sequence); # to be refined/removed in a later issue.
  if (is.null(data$month_sequence)) {
    stop("No 'data$month_sequence' value defined. \n
             Please provide a vector of months with 1 indicating the warmest month.")
  }
  # DATA_SCALAR(rho_input);
  if (is.null(data$rho_input) & data$calculate_rho == 0) {
    stop("Please provide 'data$rho_input' value when data$calculate_rho == 0.")
  }
  if (!is.null(data$rho_input) & data$calculate_rho == 1) {
    if (!silent){message(crayon::magenta("Ignoring data$rho_input. rho calculated inside the model using weight_at_recruitment and weight_inf."))}
  }
  if (is.null(data$rho_input) & data$calculate_rho == 1) {
    if (!silent){message(crayon::magenta("No 'data$rho_input' value defined. Using '1' by default but replaced by calculation in model."))}
    data$rho_input <- 1
  }
  # DATA_INTEGER(calculate_rho);
  if (is.null(data$calculate_rho)) {
    stop("No 'data$calculate_rho' value defined. \n
            Please provide a value indicating whether you want rho to be \n
           calculated inside the model using weight_at_recruitment and weight_inf.")
  }
  # DATA_INTEGER(use_recruit_penalty);
  if (is.null(data$use_recruit_penalty)) {
    stop("No 'data$use_recruit_penalty' value defined. \n
            Please provide a value indicating whether you want use the recruitment penalty.")
  }
  # DATA_SCALAR(prior_mean_mu);
  if (is.null(data$prior_mean_mu)) {
    stop("No 'data$prior_mean_mu' value defined. \n
            Please provide a value for the mean of the mu normal prior.")
  }
  # DATA_SCALAR(prior_mean_k);
  if (is.null(data$prior_mean_k)) {
    stop("No 'data$prior_mean_k' value defined. \n
            Please provide a value for the mean of the kappa normal prior.")
  }
  # DATA_SCALAR(prior_mean_xi);
  if (is.null(data$prior_mean_xi)) {
    stop("No 'data$prior_mean_xi' value defined. \n
            Please provide a value for the mean of the xi normal prior.")
  }
  # DATA_STRING(rec_dev_type);
  if (is.null(data$rec_dev_type)) {
    data$rec_dev_type <- 'off'
    if (!silent){message(crayon::magenta("No 'data$rec_dev_type' value defined. Using 'off' by default."))}
  }
  # DATA_INTEGER(Number_months_per_timestep);
  if (is.null(data$Number_months_per_timestep)) {
    if (!silent){message(crayon::magenta("No 'data$Number_months_per_timestep' value defined. Using '1' by default."))}
    data$Number_months_per_timestep <- "1"
  } else if (!data$Number_months_per_timestep %in% c(1,2,3,4,6,12)){
    stop("Invalid 'data$Number_months_per_timestep' value defined. Valid values are: 1, 2, 3, 4, 6, and 12.")
  }
  # DATA_INTEGER(use_cpue_sd);
  if (is.null(data$use_cpue_sd)) {
    if (!silent){message(crayon::magenta("No 'data$use_cpue_sd' value defined. Using '0' by default. Values in data$cpue_sd will not be used."))}
    data$use_cpue_sd <- 0
  } else if (!data$use_cpue_sd %in% c(0, 1)){
    stop("Invalid 'data$use_cpue_sd' value defined. Valid values are: 0 and 1.")
  }
  # DATA_INTEGER(num_years_warmup);
  if (is.null(data$num_years_warmup)) {
    if (!silent){message(crayon::magenta("No 'data$num_years_warmup' value defined. Using '100' by default."))}
    data$num_years_warmup <- 100
  }
  # DATA_SCALAR(recruit_penalty_strength);
  if (is.null(data$recruit_penalty_strength)) {
    if (!silent){message(crayon::magenta("No 'data$recruit_penalty_strength' value defined. Using '0.001' by default."))}
    data$recruit_penalty_strength <- 0.001
  }
  # DATA_SCALAR(recruit_penalty_exponent);
  if (is.null(data$recruit_penalty_exponent)) {
    if (!silent){message(crayon::magenta("No 'data$recruit_penalty_exponent' value defined. Using '2' by default."))}
    data$recruit_penalty_exponent <- 2
  }
  # DATA_VECTOR(absolute_biomass);
  if (is.null(data$absolute_biomass)) {
    if (!silent){message(crayon::magenta("No 'data$absolute_biomass' value defined. Vector filled with zeros, which the model will ignore."))}
    data$absolute_biomass <- 0*data$cpue
  }
    if (!is.matrix(data$absolute_biomass)) {
    if (!silent){message(crayon::magenta("Suggest making bsolute_biomass a matrix."))}
  }
  # DATA_SCALAR(absolute_biomass_sd);
  if (is.null(data$absolute_biomass_sd)) {
    if (!silent){message(crayon::magenta("No 'data$absolute_biomass_sd' value defined. Using '1' by default."))}
    data$absolute_biomass_sd <- rep(1,length(data$cpue))
  }
  # DATA_SCALAR(minimum_annual_harvest);
  if (is.null(data$minimum_annual_harvest)) {
    if (!silent){message(crayon::magenta("No 'data$minimum_annual_harvest' value defined. Using '0.05' by default."))}
    data$minimum_annual_harvest <- 0.05
  }
  # DATA_MATRIX(CoordBasis);
  if (is.null(data$CoordBasis)) {
    if (!silent){message(crayon::magenta("No 'data$CoordBasis' value defined. Creating Matrix..."))}
    number_rec_devs <- data$last_year_catch - data$first_year_rec_devs
    CoordBasis <- matrix(0,nrow=number_rec_devs, ncol=number_rec_devs+1)
    for (i in 1:number_rec_devs) {
      hh <- sqrt(0.5 * i/ (i + 1))
      CoordBasis[i, 1:i] <- -hh/ i
      CoordBasis[i, i + 1] <- hh
    }
    CoordBasis <- CoordBasis/hh
    data$CoordBasis <- t(CoordBasis)
  }
  # DATA_INTEGER(do_projections);
  if (is.null(data$do_projections)) {
    data$do_projections <- 0
    if (!silent){message(crayon::magenta("No 'data$do_projections' value defined. Using '0' by default."))}
  }
  # DATA_SCALAR(projection_years)
  if (is.null(data$projection_years)) {
    if (!silent){message(crayon::magenta("No 'data$projection_years' value defined. Using '50' by default."))}
    data$projection_years <- 50
  }
  # DATA_SCALAR(projection_harvest)
  if (is.null(data$projection_harvest)) {
    if (!silent){message(crayon::magenta("No 'data$projection_harvest' value defined. Using '0' by default."))}
    data$projection_harvest <- rep(0,data$projection_years*12/data$Number_months_per_timestep)
  }
  # DATA_SCALAR(target_relative_biomass)
  if (is.null(data$target_relative_biomass)) {
    if (!silent){message(crayon::magenta("No 'data$target_relative_biomass' value defined. Using '0.6' by default."))}
    data$target_relative_biomass <- 0.6
  }
  # DATA_SCALAR(F_initial)
  if (is.null(data$F_initial)) {
    if (!silent){message(crayon::magenta("No 'data$F_initial' value defined. Using '0.05' by default."))}
    data$F_initial <- 0.05
  }
  # DATA_SCALAR(prior_sd_mu)
  if (is.null(data$prior_sd_mu)) {
    if (!silent){message(crayon::magenta("No 'data$prior_sd_mu' value defined. Using '1 by default."))}
    data$prior_sd_mu <- 1
  }
  # DATA_SCALAR(prior_sd_k)
  if (is.null(data$prior_sd_k)) {
    if (!silent){message(crayon::magenta("No 'data$prior_sd_k' value defined. Using '1 by default."))}
    data$prior_sd_k <- 1
  }
  # DATA_SCALAR(prior_sd_xi)
  if (is.null(data$prior_sd_xi)) {
    if (!silent){message(crayon::magenta("No 'data$prior_sd_xi' value defined. Using '1 by default."))}
    data$prior_sd_xi <- 1
  }
  # DATA_STRUCT(economic_data)
  if (!is.list(data$economic_data)) {
    if (!silent){message(crayon::magenta("Economic_data not provided as list. Setting default values."))}
    data$economic_data <- list(cK = 1)
  }

  for (value in c('cL', 'cM', 'cK', 'cF', 'cO', 'v', 'W', 'K', 'prop', 'B', 'dmean', 'o', 'd')) {
    default_values <- list(cL = 0.29, cM = 0.41, cK = 384, cF = 619, cO = 44.26, 
                           v = rep(100,12/data$Number_months_per_timestep), W = 42646, 
                           K = 589719, prop = 0.63, B = 143, dmean = 68, o = 0.05, d = 0.02)

    if (is.null(data$economic_data[value][[1]])) {
      if (!silent){message(crayon::magenta(paste0("No 'data$economic_data$",value,"' defined. Using ", default_values[value], " by default.")))}
        data$economic_data[value] = default_values[value] 
    }
  }

  if (!length(data$economic_data$v) == 12/data$Number_months_per_timestep){
    stop(paste0("Length of data$economic_data$v should be ",12/data$Number_months_per_timestep,". Provide the average price/value for every timestep (",12/data$Number_months_per_timestep," timestep/s)."))
  }

  # Economic_data vector size

  # Data vector size
  num_timesteps <- (data$last_year_catch - data$first_year_catch + 1)*12/data$Number_months_per_timestep
  lengths <-  c(cpue = ncol(data$cpue),
                # cpue_sd = ncol(data$cpue_sd),
                ctch = length(data$ctch),
                absolute_biomass = ncol(data$absolute_biomass),
                absolute_biomass_sd = length(data$absolute_biomass_sd))
  if (!length(unique(lengths))==1 || !unique(lengths==num_timesteps)) {
    if (!silent){stop(paste0(
      "Vectors must have", num_timesteps, "columns.",
      "'data$cpue' has ", lengths['cpue'], " columns.
      'data$ctch' has ", lengths['ctch'], " columns.
      'data$absolute_biomass' has ", lengths['absolute_biomass'], " columns.
      'data$absolute_biomass_sd' has ", lengths['absolute_biomass_sd'], " columns"))}
  }

   # Projection harvest vector size
  if (!length(data$projection_harvest) == data$projection_years*12/data$Number_months_per_timestep){
    if (!silent){stop(paste0("data$projection_harvest must have length data$projection_years (", 
                              data$projection_years, ")"))}
  }
 
  # rec dev lengths
  num_recdev_years <- data$last_year_catch - data$first_year_rec_devs + 1
  if (!all(dim(data$CoordBasis) == c(num_recdev_years,num_recdev_years-1))) {stop(paste0('data$CoordBasis must have dimensions ',num_recdev_years,' by ', num_recdev_years-1))}
  if (!silent){message(crayon::magenta(paste0("Expecting 'parameters$log_R_star' to have length ",num_recdev_years,".")))}
  if (!silent){message(crayon::magenta(paste0("Expecting 'parameters$zeta' to have length ",num_recdev_years-1,".")))}
  number_rec_devs <- data$last_year_catch - data$first_year_rec_devs

  return(data)
}
