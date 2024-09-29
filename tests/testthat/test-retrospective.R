# Copyright 2024 Fisheries Queensland

# This file is part of DDUST.
# DDUST is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DDUST is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with DDUST. If not, see <https://www.gnu.org/licenses/>.

test_that("multiplication works", {
    data <- DDUST::data
    parameters <- DDUST::parameters
    map <- DDUST::map
    data$cpue_sd <- data$cpue_sd[,1:length(data$cpue)] |> as.matrix()  |> t()

    dd_out <- retrospective(data, parameters, map, peel = 5)

    expect_false(any(sapply(dd_out$dd_mle, function(list_element) {
    return(is.numeric(list_element) && any(is.nan(list_element)))
    })), info = "The function result contains NaN values")
})
