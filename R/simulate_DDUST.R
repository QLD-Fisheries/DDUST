# Copyright 2024 Fisheries Queensland

# This file is part of DDUST.
# DDUST is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DDUST is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with DDUST. If not, see <https://www.gnu.org/licenses/>.

#' Simulate DDUST model with MCMC parameters to extract derived quantities
#'
#' @param dd_mle List of MLE outputs
#' @param dd_mcmc List of MCMC outputs
#'
#' @return List of simulated quantities
#' @export
#'
#' @examples
#' \dontrun{
#' dd_out <- run_DDUST(data, parameters, map, MCMC = TRUE)
#' dd_sim <- simulate_DDUST(dd_out$dd_mle, dd_out$dd_mcmc)
#' }
#'
simulate_DDUST <- function(dd_mle, dd_mcmc){

  if (!is.list(dd_mle) || !all(sapply(dd_mle,is.list))){
    dd_mle <- list(dd_mle)
  }
  if (!is.list(dd_mcmc) || !all(sapply(dd_mcmc, inherits, "stanfit"))){
    dd_mcmc <- list(dd_mcmc)
  }

  simulation_list <- list()
  for (scenario in 1:length(dd_mle)){
    fit <- as.matrix(dd_mcmc[[scenario]])
    simulation_list[[scenario]] <- sapply(1:nrow(fit), function(i){list(dd_mle[[scenario]]$model$simulate(fit[i,1:ncol(fit)-1]))})
  }
  return(simulation_list)
}
