# bare-bones test data ----
## Monthly model
test_data <- list()
test_data$Number_months_per_timestep <- 1
test_data$proportion_spawning <- rep(1/12,12)
test_data$weight_at_recruitment <- c(1,2)
test_data$weight_inf <- c(5)
test_data$last_year_catch <- 2023
test_data$first_year_catch <- 1950
test_data$first_year_rec_devs <- 1950
test_data$month_sequence <- as.numeric(1:12)
test_data$calculate_rho <- 0
test_data$rho_input <- 0.606
test_data$use_recruit_penalty <- 0
test_data$prior_mean_xi <- log(3)
test_data$prior_mean_mu <- 0
test_data$prior_mean_k <- 5
test_data$ctch <- c(seq(0,1000,length.out=round(888*0.5)),seq(1000,5,length.out=round(888*0.5)))
set.seed(123)
test_data$cpue <- rnorm(888,30,5)-c(seq(0,15,length.out=round(888*0.75)),seq(15,5,length.out=round(888*0.25)))
test_data$cpue_sd <- t(matrix(rep(0,888)))
test_data$ctch <- test_data$ctch + rnorm(length(test_data$ctch),0, test_data$ctch/3)
test_data$ctch <- rowSums(matrix(test_data$ctch,ncol=test_data$Number_months_per_timestep, byrow=T))
test_data$cpue <- t(matrix(rowSums(matrix(test_data$cpue,ncol=test_data$Number_months_per_timestep, byrow=T))/test_data$Number_months_per_timestep))

test_pars <- list()
test_pars$Rinit <- 15
test_pars$xi <- log(0.5)
test_pars$M <- 1.2
test_pars$k <- 2
test_pars$mu <- 0
test_pars$q1 <- 0.15
test_pars$q2 <- 0.01
test_pars$lsigmaR_sq <- log(0.36^2)
test_pars$lsigmaI_sq <- log(0.2^2)
set.seed(123)
test_pars$zeta <- rnorm(2023 - 1950, 0, 0.36)
test_pars$log_R_star <- as.numeric(rep(19,74))

test_map <- list(M=factor(NA),
            q1 = factor(NA),
            q2 = factor(NA),
            k = factor(NA),
            mu = factor(NA),
            lsigmaI_sq=factor(NA),
            lsigmaR_sq=factor(NA),
            xi=factor(NA),
            log_R_star = rep(factor(NA),length(test_pars$log_R_star)),
            zeta = rep(factor(NA),length(test_pars$zeta)))

# dd_out <- run_DDUST(test_data, test_pars, test_map, laplace = FALSE)

save(test_data,test_pars,test_map,file='tests/testthat/fixtures/test.rda')
