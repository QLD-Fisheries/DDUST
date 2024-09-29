# Copyright 2024 Fisheries Queensland

# This file is part of DDUST.
# DDUST is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DDUST is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with DDUST. If not, see <https://www.gnu.org/licenses/>.

#' Check parameters for DDUST model
#'
#' @param parameters list of parameters
#' @param rec_dev_type Type of recruitment deviations: "barycentric", "redundant" or "none"
#' @param silent TRUE to receive more information
#'
#' @return updated list of parameters
#' @export
#'
#' @examples
#' parameters <- check_parameters(parameters, rec_dev_type = 'redundant')
check_parameters <- function(parameters, rec_dev_type, silent = FALSE){

  if (is.null(parameters$xi)){
    if (is.null(parameters$h)){
      if (!silent){message(crayon::magenta("No 'parameters$xi' or 'parameters$h' value defined. Using h = 0.7 by default."))}
      parameters$xi <- log( (1-5*0.7) / (0.7-1) )
    } else {
      parameters$xi <- log((1-5*parameters$h) / (parameters$h-1))
    }

  }
  if (rec_dev_type == 'barycentric'){
    if (is.null(parameters$zeta)){
      stop("No 'parameters$zeta' value defined. \n
             Please provide the initial value for the recruitment deviations.")
    }
    if (!is.null(parameters$log_R_star)){
      if (!silent){message(crayon::magenta('Barycentric recruitment deviations specified, data$log_R_star will be ignored.'))}
    } else {
      parameters$log_R_star <- 0*parameters$zeta
    }
  }
  if (rec_dev_type == 'redundant'){
    if (is.null(parameters$log_R_star)) {
      stop("No 'parameters$log_R_star' value defined. \n
           Did you mean to use recruitment deviations 'redundant'? \n
           Please provide the initial value for the recruitment deviations using parameters$log_R_star or set data$rec_dev_type = 'barycentric'.")
    }
    if (!is.null(parameters$zeta)){
      if (!silent){message(crayon::magenta('Redundant recruitment deviations specified, data$zeta will be ignored.'))}
    } else {
      parameters$zeta <- rnorm(length(parameters$log_R_star)-1,0,sqrt(exp(parameters$lsigmaR_sq)))
    }
  }

  return(parameters)
}
