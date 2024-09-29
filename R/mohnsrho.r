# Copyright 2024 Fisheries Queensland

# This file is part of DDUST.
# DDUST is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DDUST is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with DDUST. If not, see <https://www.gnu.org/licenses/>.

#' Calculate Mohn's rho
#' 
#' @description 
#' From 'The retrospective problem in sequential population analysis: An investigation using cod fishery and simulated data.' Mohn, R. (1999)
#' 
#' @param dd_out List of dd_out. Output from retrospective(). First list element must be model with full data timeseries.
#'
#' @return A value

mohnsrho <- function(dd_out){
    rho <- 0
    for (i in 2:length(dd_out)) {
        index <- length(dd_out[[i]]$dd_mle$B_annual)
        rho <- rho + (dd_out[[i]]$dd_mle$B_annual[index] - dd_out[[1]]$dd_mle$B_annual[index])/dd_out[[1]]$dd_mle$B_annual[index]
    }
    return(rho/(length(dd_out)-1))
}