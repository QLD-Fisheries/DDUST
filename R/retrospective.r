# Copyright 2024 Fisheries Queensland

# This file is part of DDUST.
# DDUST is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DDUST is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with DDUST. If not, see <https://www.gnu.org/licenses/>.

#' Retrospective analysis
#' 
#' @param data list of data. See ?check_data
#' @param parameters list of parameters
#' @param map list of fixed parameters
#' @param peel number of years to peel (int)
#' @param silent TRUE to silence console output
#' @param laplace TRUE to use laplace approximation for recruitment deviation parameters
#' @param lower upper bounds on the variables passed to nlminb()
#' @param upper lower bounds on the variables passed to nlminb()
#'
#' @return List of MLE outputs. `dd_mle` is a list of model inputs and outputs from nlminb(). 
#' @export
#'
#' @examples
#'
#' dd_out <- retrospective(data, parameters, map)
#'
retrospective <- function(data, parameters, map, peel = 5, silent = FALSE, laplace = TRUE, lower = NULL, upper = NULL){

    peel_length <- 12/data$Number_months_per_timestep
    dd_out <- list()

    for (i in 0:peel){
        if (!silent){message(crayon::blue(paste0("Running model with ", i, " years removed...")))}
        data$last_year_catch <- data$last_year_catch - i
        data$ctch <- data$ctch[1:(length(data$ctch)-i*peel_length)]
        data$cpue <- data$cpue[,1:(length(data$cpue)-i*peel_length)]  |> as.matrix()  |> t()
        data$cpue_sd <- data$cpue_sd[,1:(length(data$cpue_sd)-i*peel_length)] |> as.matrix()  |> t()
        data$do_projections <- 0
        data$absolute_biomass <- data$absolute_biomass[,1:(length(data$absolute_biomass)-i*peel_length)] |> as.matrix()  |> t()
        data$absolute_biomass_sd <- data$absolute_biomass_sd[1:(length(data$absolute_biomass_sd)-i*peel_length)]
        data$CoordBasis <- NULL

        parameters$zeta <- parameters$zeta[1:(length(parameters$zeta)-i)]
        parameters$log_R_star <- parameters$log_R_star[1:(length(parameters$log_R_star)-i)]

        map$zeta <- map$zeta[1:(length(map$zeta)-i)]
        map$log_R_star <- map$log_R_star[1:(length(map$log_R_star)-i)]

        dd_out[[i+1]] <- run_DDUST(data, parameters, map, silent = TRUE, laplace = laplace, lower = lower, upper = upper)
        if (!silent) message(crayon::white("\U2714 Success \n"))
    }
    if (!silent) message(crayon::blue("Mohn's rho: ", mohnsrho(dd_out)))
    return(dd_out)
}
