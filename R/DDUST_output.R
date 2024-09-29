# Copyright 2024 Fisheries Queensland

# This file is part of DDUST.
# DDUST is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DDUST is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with DDUST. If not, see <https://www.gnu.org/licenses/>.

#' Condense sdreport, fit, data, parameters and map into one full report
#'
#' @param rep sdreport(model)
#' @param fit optim(model)
#' @param data data used in MakeADFun
#' @param parameters parameters used in MakeADFun
#' @param map map used in MakeADFun
#' @param model model created with MakeADFun
#'
#' @return A dataframe
#' @export
#'
DDUST_output <- function(rep, fit, data, parameters, map, model){

  report <- rep$value
  report_sd <- rep$sd
  output_names = unique(names(report))
  X = list()
  for (i in 1:length(output_names)){
    temp = which(names(report) == output_names[i])
    X[[output_names[i]]] = report[names(report) == output_names[i]]
    X[[paste0(output_names[i],"_sd")]] = report_sd[temp]
  }

  output <- X
  output$convergence = fit$convergence
  output$hessian = rep$pdHess
  output$data = data
  output$parameters = parameters
  output$map = map
  output$model = model
  return(output)
}
