## Copyright (C) 2018-present, Jiali Lin
## All rights reserved.

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

testthat::context("Unit tests for sample_data.R")

test_that("Precision matrix must be positive definite", {
  set.seed(99)
  p <- 10
  K <- 4
  n <- 400
  dat <- GenerateData(p, K, n)
  A = dat$A
  expect_true(all(eigen(A)$values > 0))
})