# Copyright 2024 Fisheries Queensland

# This file is part of DDUST.
# DDUST is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DDUST is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with DDUST. If not, see <https://www.gnu.org/licenses/>.

#' Make the DDUST TMB model using the data, parameters and map supplied
#'
#' @param data list of data
#' @param parameters list of parameters
#' @param map list of fixed parameters. Used in TMB::MakeADFun()
#' @param laplace TRUE to use laplace approximation for recruitment deviation parameters
#' @param hessian calculate hessian: TRUE or FALSE
#'
#' @return model object
#' @export
#'
#' @examples
#' data <- check_data(data)
#' parameters <- check_parameters(parameters, rec_dev_type = 'redundant')
#'
#' make_DD_model(data, parameters, map)
make_DD_model <- function(data, parameters, map, laplace = TRUE, hessian = TRUE){

  if (laplace) {
    if (data$rec_dev_type == 'redundant'){
    model <- TMB::MakeADFun(data = c(model = "DDUST",data),
                       parameters,
                       random = "log_R_star",
                       DLL = "DDUST_TMBExports",
                       map = map,
                       hessian = hessian,
                       checkParameterOrder=TRUE,
                       silent = FALSE) 
    } else if (data$rec_dev_type == 'barycentric'){
      model <- TMB::MakeADFun(data = c(model = "DDUST",data),
                       parameters,
                       random = "zeta",
                       DLL = "DDUST_TMBExports",
                       map = map,
                       hessian = hessian,
                       checkParameterOrder=TRUE,
                       silent = FALSE) 
    } else {
      stop('To use no recruitment deviations, turn off laplace approximation with `laplace = FALSE`.')
    }
  } else {
    model <- MakeADFun(data = c(model = "DDUST",data),
                       parameters,
                       DLL = "DDUST_TMBExports",
                       map = map,
                       hessian = hessian,
                       checkParameterOrder=TRUE,
                       silent = FALSE)
  }
  return(model)
}
