# Copyright 2024 Fisheries Queensland

# This file is part of DDUST.
# DDUST is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DDUST is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with DDUST. If not, see <https://www.gnu.org/licenses/>.

#' Quickly run DDUST
#' 
#' 1. Model data is checked with check_data().
#' 2. Model parameters is checked with check_parameters()
#' 3. DDUST model is formulated with make_DD_model(). Laplace approximation is used by default.
#' 4. Model is optimised with nlminb().
#' 5. Standard deviation report is calculated with sdreport().
#' 6. Outputs are collated with DDUST_output().
#'
#' @param data list of data. See ?check_data
#' @param parameters list of parameters
#' @param map list of fixed parameters
#' @param silent TRUE to silence console output
#' @param laplace FALSE to not use laplace approximation for recruitment deviation parameters
#' @param MCMC TRUE to run MCMC after MLE
#' @param chains number of chains, passed to tmbstan()
#' @param iter number of iterations, passed to tmbstan()
#' @param maxit maximum number of iterations passed to nlminb()
#' @param init initial values for the sampler, passed to tmbstan()
#' @param lower upper bounds on the variables passed to nlminb()
#' @param upper lower bounds on the variables passed to nlminb()
#' @param ... other parameters to be passed to tmbstan()
#'
#' @return List of MLE and MCMC outputs. `dd_mle` is a list of model inputs and outputs from optim(). `dd_mcmc` is a list of MCMC outputs.
#' @export
#'
#' @examples
#'
#' # Quickly run DDUST
#' dd_out <- run_DDUST(data, parameters, map)
#'
#' # Step-by-step run DDUST
#' \dontrun{
#' # Check data and parameters
#' data <- check_data(data, silent = TRUE)
#' parameters <- check_parameters(parameters, rec_dev_type = data$rec_dev_type, silent = TRUE)
#'
#' # Make TMB model
#' model <- TMB::MakeADFun(data = c(model = "DDUST",data),
#'                         parameters,
#'                         random = "log_R_star",
#'                         DLL = "DDUST_TMBExports",
#'                         map = map,
#'                         hessian = TRUE,
#'                         checkParameterOrder=FALSE)
#'
#' # Optimise
#' fit <- optim(model$par, model$fn, model$gr)
#'
#' # TMB report
#' rep <- TMB::sdreport(model)
#'
#' # Collate objects
#' dd_mle <- DDUST_output(rep, fit, data, parameters, map, model)
#'}
run_DDUST <- function(data, parameters, map, silent = FALSE, laplace = FALSE, maxit = 500000, MCMC = FALSE, chains = 5, iter = 1000, init = 'par', lower = NULL, upper = NULL, ...){

  if (!silent) message(crayon::blue("Checking data... \n"))
  data <- check_data(data, silent = silent)
  if (!silent) message(crayon::white("\U2714 Success \n"))
  if (!silent) message(crayon::blue("Checking parameters...\n"))
  parameters <- check_parameters(parameters, rec_dev_type = data$rec_dev_type, silent = silent)
  if (!silent) message(crayon::white("\U2714 Success \n"))
  if (is.null(map$zeta) & is.null(map$log_R_star)) {
    stop('Do not estimate zeta and log_R_star in the same model. Fix one of the parameters using map$zeta or map$log_R_star.')
  }
  if (!silent) message(crayon::blue("Making DDUST model... \n"))
  model <- make_DD_model(data, parameters, map, laplace)
  if (!silent){
    message(crayon::white("\U2714 Success \n"))
    message(crayon::blue("Optimising... \n"))
    message(crayon::magenta(paste0("Estimated parameters: ")))
    message(crayon::magenta(paste0(unique(names(model$par))," \n")))
    message(crayon::magenta(paste0("Fixed parameters: ")))
    message(crayon::magenta(paste0(unique(setdiff(names(parameters),names(model$par)))," \n")))
    message(crayon::magenta(paste0("Parameters using laplace approximation: ")))
    message(crayon::magenta(paste0(unique(names(model$env$par[model$env$random]))," \n")))
  }
  sink("NULL")
  fit <- nlminb(model$par, model$fn, model$gr, lower = lower, upper = upper)
  sink()
  if (!silent) message(crayon::white("\U2714 Success \n"))
  if (!silent) message(crayon::blue("Calculating standard deviation report... \n"))
  sink("NULL")
  if (!laplace){
    rep <- TMB::sdreport(model, hessian.fixed = model$he())
  } else {
    rep <- TMB::sdreport(model)
  }
  sink()
  if (!silent) message(crayon::white("\U2714 Success \n"))
  dd_mle <- DDUST_output(rep, fit, data, parameters, map, model)

  if (!MCMC){return(list(dd_mle = dd_mle,
                         dd_mcmc = NULL))}

  if (MCMC){
    options(mc.cores = parallel::detectCores())
    dd_mcmc <- tmbstan::tmbstan(dd_mle$model, chains = chains, iter = iter, init = init, ...)
    return(list(dd_mle = dd_mle,
                dd_mcmc = dd_mcmc))
  }
}
