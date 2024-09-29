# Copyright 2024 Fisheries Queensland

# This file is part of DDUST.
# DDUST is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DDUST is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with DDUST. If not, see <https://www.gnu.org/licenses/>.

# This script performs tests for different combinations of options used in DDUST. 
# Not all tests have been written yet.
# The dimensions for testing are:
# a)	Timestep (monthly, bimonthly, trimonthly, quadmonthly, semi-annual and annual)
# b)	Recruitment deviation type (barycentric, as is, off)
# c)	Laplace approximation (yes, no)
# d)	Late start recruitment deviations (yes, no)

# The convention for naming tests should specify which combination of these options is being tested:
# For example,
#         Timestep monthly, Recdevs barycentric, laplace
#         Timestep annual, Recdevs redundant, no laplace, late start

# Bare-bones test data is provided in fixtures/test.rda

# Projections OFF ----

## Monthly timestep ----

### rec_dev_type 'off' ----

test_that("Test 2", { # Monthly timestep, recdevs off, laplace, projections off"
  load("fixtures/test.rda")
  expect_error(run_DDUST(test_data, test_pars, test_map, laplace = TRUE, silent = TRUE),
  info = 'Successfully detects incorrect use of laplace.')
})

test_that("Test 4", { # Monthly timestep, recdevs off, no laplace, projections off
  load("fixtures/test.rda")
  data <- test_data
  parameters <- test_pars
  map <- test_map
  data$use_recruit_penalty <- 1

  data <- check_data(data, silent = TRUE)
  parameters <- check_parameters(parameters, rec_dev_type = data$rec_dev_type, silent = TRUE)
  model <- make_DD_model(data, parameters, map, laplace=FALSE)
  model$env$tracepar <- TRUE
  sink("NULL")
  fit <- optim(model$par, model$fn, model$gr, method = "Brent", lower = 9, upper = 20, control=list(trace=1))
  sink()
  sink("NULL")
  rep <- TMB::sdreport(model)
  sink()
  dd_mle <- DDUST_output(rep, fit, data, parameters, map, model)
  dd_out <- list(dd_mle = dd_mle)

  expect_false(any(sapply(dd_out$dd_mle, function(list_element) {
    return(is.numeric(list_element) && any(is.nan(list_element)))
  })), info = "The function result contains NaN values")
})

### rec_dev_type 'barycentric' ----

test_that("Test 6", { # Monthly timestep, recdevs barycentric, laplace, projections off
  load("fixtures/test.rda")
  test_data$rec_dev_type <- "barycentric"
  test_map$zeta <- NULL
  test_data$ctch <- test_data$ctch/12
  test_pars$lsigmaR_sq <- log(0.1^2)
  
  dd_out <- run_DDUST(test_data, test_pars, test_map, silent = TRUE, upper = 11)

  expect_false(any(sapply(dd_out$dd_mle, function(list_element) {
    return(is.numeric(list_element) && any(is.nan(list_element)))
  })), info = "The function result contains NaN values")
})

### rec_dev_type 'redundant' ----

## Annual timestep ----

### rec_dev_type 'off' ----

test_that("Test 15", { # Annual timestep, recdevs off, no laplace, late start
  data <- DDUST::data
  parameters <- DDUST::parameters
  map <- DDUST::map

  data$rec_dev_type <- "off"
  data$first_year_rec_devs <- 1970
  data$CoordBasis <- NULL
  n <- data$last_year_catch - data$first_year_rec_devs
  parameters$log_R_star <- rep(19, n + 1)
  parameters$zeta <- rnorm(n, 0, 0.2)
  map$log_R_star <- rep(factor(NA), n + 1)
  map$zeta <- rep(factor(NA), n)
  map$lsigmaR_sq <- factor(NA)

  data <- check_data(data, silent = FALSE)
  parameters <- check_parameters(parameters, "barycentric", silent = FALSE)

  dd_out <- run_DDUST(data, parameters, map, laplace = FALSE, silent = FALSE)

  expect_false(any(sapply(dd_out$dd_mle, function(list_element) {
    return(is.numeric(list_element) && any(is.nan(list_element)))
  })), info = "The function result contains NaN values")
})

### rec_dev_type 'barycentric' ----

test_that("Test 17", { # Annual timestep, recdevs barycentric, laplace, late start, projections off
  data <- DDUST::data
  parameters <- DDUST::parameters
  map <- DDUST::map

  data$rec_dev_type <- "barycentric"
  data$first_year_rec_devs <- 1970
  data$CoordBasis <- NULL
  n <- data$last_year_catch - data$first_year_rec_devs
  parameters$log_R_star <- rep(19, n + 1)
  parameters$zeta <- rnorm(n, 0, 0.2)
  map$log_R_star <- rep(factor(NA), n + 1)
  map$zeta <- NULL
  map$lsigmaR_sq <- factor(NA)
  map$xi <- factor(NA)

  dd_out <- run_DDUST(data, parameters, map, silent = TRUE)

  expect_false(any(sapply(dd_out$dd_mle, function(list_element) {
    return(is.numeric(list_element) && any(is.nan(list_element)))
  })), info = "The function result contains NaN values")
})

### rec_dev_type 'redundant' ----

test_that("Test 21", { # Annual timestep, recdevs redundant, laplace, late start, projections off
  data <- DDUST::data
  parameters <- DDUST::parameters
  map <- DDUST::map

  data$rec_dev_type <- "redundant"
  data$first_year_rec_devs <- 1970
  data$CoordBasis <- NULL

  n <- data$last_year_catch - data$first_year_rec_devs

  parameters$log_R_star <- rep(19, n + 1)
  parameters$zeta <- rnorm(n, 0, 0.2)

  map$log_R_star <- NULL
  map$zeta <- rep(factor(NA), n)
  map$lsigmaR_sq <- factor(NA)

  data <- check_data(data, silent = TRUE)
  parameters <- check_parameters(parameters, data$rec_dev_type, silent = TRUE)
  
  dd_out <- run_DDUST(data, parameters, map, silent = TRUE)

  expect_false(any(sapply(dd_out$dd_mle, function(list_element) {
    return(is.numeric(list_element) && any(is.nan(list_element)))
  })), info = "The function result contains NaN values")
})


# Projections ON ----

## Monthly timestep ----

### rec_dev_type 'off' ----

test_that("Test 28", { # Monthly timestep, recdevs off, no laplace, projections on
  load("fixtures/test.rda")
  data <- test_data
  parameters <- test_pars
  map <- test_map
  data$use_recruit_penalty <- 1
  data$do_projections <- 1

  data <- check_data(data, silent = TRUE)
  parameters <- check_parameters(parameters, rec_dev_type = data$rec_dev_type, silent = TRUE)
  model <- make_DD_model(data, parameters, map, laplace=FALSE)
  model$env$tracepar <- TRUE
  sink("NULL")
  fit <- optim(model$par, model$fn, model$gr, method = "Brent", lower = 9, upper = 20, control=list(trace=1))
  sink()
  sink("NULL")
  rep <- TMB::sdreport(model)
  sink()
  dd_mle <- DDUST_output(rep, fit, data, parameters, map, model)
  dd_out <- list(dd_mle = dd_mle)

  expect_false(any(sapply(dd_out$dd_mle, function(list_element) {
    return(is.numeric(list_element) && any(is.nan(list_element)))
  })), info = "The function result contains NaN values")
})

### rec_dev_type 'barycentric' ----

### rec_dev_type 'redundant' ----

## Annual timestep ----

### rec_dev_type 'off' ----

### rec_dev_type 'barycentric' ----

test_that("Test 42", { # Annual timestep, recdevs barycentric, laplace, projections on
  tol <- 0.00001
  data <- DDUST::data
  parameters <- DDUST::parameters
  map <- DDUST::map
  data$do_projections <- 1

  # Use projection_harvest
  data$projection_harvest = rep(10000,data$projection_years*12/data$Number_months_per_timestep)

  # Run model
  dd_out <- run_DDUST(data, parameters, map, MCMC = FALSE, silent = TRUE)
  dd_out$dd_mle$msy
  # No NAN values
  expect_false(any(sapply(dd_out$dd_mle, function(list_element) {
    return(is.numeric(list_element) && any(is.nan(list_element)))
  })))

  # Parameter results
  expect_equal(dd_out$dd_mle$h, c(h = 0.227021), tolerance = tol)
  expect_equal(dd_out$dd_mle$R0, c(R0 = 108461.7), tolerance = tol)

  # Biomass results
  expect_equal(tail(dd_out$dd_mle$B_annual_ratio, 1), c(B_annual_ratio = 0.8355458), tolerance = tol)

  # Yield results
  expect_equal(dd_out$dd_mle$msy, c(msy = 616062), tolerance = tol)
  expect_equal(dd_out$dd_mle$Ftarg, c(Ftarg = 0.08036292), tolerance = 0.0001)

  # Priors
  expect_equal(dd_out$dd_mle$prior_k, c(prior_k = 6.125), tolerance = tol)
  expect_equal(dd_out$dd_mle$prior_mu, c(prior_mu = 6.125), tolerance = tol)
  expect_equal(dd_out$dd_mle$prior_xi, c(prior_xi = 4.040785), tolerance = tol)

  # Likelihood
  expect_equal(dd_out$dd_mle$RecDevLL, c(RecDevLL = -184.1358), tolerance = tol)
  expect_equal(dd_out$dd_mle$cpueLL, c(cpueLL = -113.828), tolerance = tol)
  expect_equal(dd_out$dd_mle$biomassLL, c(biomassLL = 0), tolerance = tol)
  expect_equal(dd_out$dd_mle$penLL1, c(penLL1 = 0), tolerance = tol)
  expect_equal(dd_out$dd_mle$penLL2, c(penLL2 = 50.71449), tolerance = tol)
  expect_equal(dd_out$dd_mle$LL, c(LL = -230.9585), tolerance = tol)

  # Projection
  expect_equal(tail(dd_out$dd_mle$HarvestProjection.Biomass,1), 
              c(HarvestProjection.Biomass = 181513.8), tolerance = tol)

})

test_that("Test 44", { # Annual timestep, recdevs barycentric, no laplace, projections on
  data <- DDUST::data
  parameters <- DDUST::parameters
  map <- DDUST::map
  data$do_projections <- 1

  data$rec_dev_type <- "barycentric"
  map$zeta <- NULL
  map$lsigmaR_sq <- factor(NA)
  map$log_R_star <- rep(factor(NA), length(parameters$log_R_star))
  map$xi <- factor(NA)
  parameters$Rinit <- 11
  parameters$xi <- log(3)

  # dd_out <- run_DDUST(data, parameters, map, silent = TRUE)
  data <- check_data(data, silent = TRUE)
  parameters <- check_parameters(parameters, rec_dev_type = data$rec_dev_type, silent = TRUE)
  model <- make_DD_model(data, parameters, map, laplace = FALSE)
  sink('NULL')
  fit <- nlminb(model$par, model$fn, model$gr)
  rep <- TMB::sdreport(model)
  sink()
  dd_mle <- DDUST_output(rep, fit, data, parameters, map, model)

  tol <- 0.00001

  # Parameter results
  expect_equal(dd_mle$h, c(h = 0.5), tolerance = tol)
  expect_equal(dd_mle$R0, c(R0 = 50903), tolerance = tol)

  # Biomass results
  expect_equal(tail(dd_mle$B_annual_ratio, 1), c(B_annual_ratio = 1.0199878), tolerance = tol)

  # Yield results
  expect_equal(dd_mle$msy, c(msy = 1535493), tolerance = tol)
  expect_equal(dd_mle$Ftarg, c(Ftarg = 0.7722452), tolerance = tol)

  # Priors
  expect_equal(dd_mle$prior_k, c(prior_k = 6.125), tolerance = tol)
  expect_equal(dd_mle$prior_mu, c(prior_mu = 6.125), tolerance = tol)
  expect_equal(dd_mle$prior_xi, c(prior_xi = 0), tolerance = tol)

  # Likelihood
  expect_equal(dd_mle$RecDevLL, c(RecDevLL = -3.608457), tolerance = tol)
  expect_equal(dd_mle$cpueLL, c(cpueLL = -115.3222), tolerance = tol)
  expect_equal(dd_mle$biomassLL, c(biomassLL = 0), tolerance = tol)
  expect_equal(dd_mle$penLL1, c(penLL1 = 0), tolerance = tol)
  expect_equal(dd_mle$penLL2, c(penLL2 = 30.53004), tolerance = tol)
  expect_equal(dd_mle$LL, c(LL = -76.1506), tolerance = tol)
})

### rec_dev_type 'redundant' ----
