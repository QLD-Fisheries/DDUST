
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DDUST <img src="man/figures/logo.png" align="right" height="138" alt="" />

<!-- badges: start -->
<!-- badges: end -->

A delay difference R package that allows the user to specify the time
step used for delays and incorporates seasonal variation in recruitment,
spawning, and catchability. The delay difference with user specified
time step (DDUST) model allows for monthly, bimonthly, trimonthly,
quadmonthly, semi-annual and annual biomass dynamics.

For any questions or inquiries, please email
<fisheriesassessment@daf.qld.gov.au>.

## Installation

You can install the development version of DDUST from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("QLD-Fisheries/DDUST")
```

## Example

Once DDUST is installed, the user will create a list of data and
parameters. The easiest way to use DDUST is with the function
`run_DDUST()` which will 1) check the inputs objects with `check_data()`
and `check_parameters()`, 2) use the data and parameters to build a TMB
model with `TMB::makeADreport()` or `make_DD_model()`, 3) optimise using
`stats::optim()`, 4) extract parameter uncertainty with
`TMB::sdreport()` and 5) collate all relevant objects using
`DDUST_output()`. For more customisation, the user can run each step of
`run_DDUST()` separately or replace with their own preferences.

The DDUST package also comes with example data and example parameters to
help with set up.

``` r
# Basic use:
dd_out <- run_DDUST(data,parameters,map)
```

``` r
# Customised use:

# Check data and parameters
data <- check_data(data, silent = TRUE)
parameters <- check_parameters(parameters, rec_dev_type = data$rec_dev_type, silent = TRUE)

# Make TMB model
model <- TMB::MakeADFun(data = c(model = "DDUST",data),
                        parameters,
                        random = "log_R_star",
                        DLL = "DDUST_TMBExports",
                        map = map,
                        hessian = hessian,
                        checkParameterOrder=FALSE)

# Optimise
fit <- optim(model$par, model$fn, model$gr, method = 'L-BFGS-B', control = list(maxit = 50000))

# TMB report
rep <- sdreport(model)

# Collate objects
dd_mle <- DDUST_output(rep, fit, data, parameters, map, model)

# MCMC with tmbstan
dd_mcmc <- tmbstan::tmbstan(dd_mle$model)

# Collect results
dd_out <- list(dd_mle = dd_mle,
               dd_mcmc = dd_mcmc)
```

The best way to plot results is to use the Fisheries Queensland plotting
package SSAND. You can install the development version of SSAND from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("QLD-Fisheries/SSAND")
```

Some plots are for comparing multiple models. This example shows a
spaghetti plot for six (6) scenarios:

``` r
library(SSAND)
SSAND::spaghettiplot(spaghettiplot_prep_DD(dd_mle))
```

<img src="man/figures/README-example-1.png" width="100%" />
