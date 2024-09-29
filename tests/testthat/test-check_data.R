# Copyright 2024 Fisheries Queensland

# This file is part of DDUST.
# DDUST is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DDUST is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with DDUST. If not, see <https://www.gnu.org/licenses/>.

test_that("recognise ill-formatted data", {
  expect_error(check_data(NULL, silent = TRUE))
})


test_that("Projection harvest vector too short", {
  data <- DDUST::data
  parameters <- DDUST::parameters
  map <- DDUST::map
  data$do_projections <- 1
  data$projection_harvest <- rep(0,10)

  expect_error(check_data(data))
})

test_that("Wrong CoordBasis dimensions", {
  data <- DDUST::data
  parameters <- DDUST::parameters
  map <- DDUST::map
  data$CoordBasis <- matrix(0,nrow=2,ncol=2)
  expect_error(check_data(data))
})

