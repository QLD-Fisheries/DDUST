# Copyright 2024 Fisheries Queensland

# This file is part of DDUST.
# DDUST is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DDUST is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with DDUST. If not, see <https://www.gnu.org/licenses/>.

#' Convert report output from TMB::sdreport() to dataframe
#'
#' @param report A list with eight elements created by TMB::sdreport()
#'
#' @return A dataframe
#' @export
#'
report2df <- function(report){

  output_names = unique(names(report))
  X <- list()
  for (i in 1:length(output_names)){
    X[[output_names[i]]] <- report[names(report) == output_names[i]]
  }
  return(X)

}
