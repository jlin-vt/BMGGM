## Copyright (C) 2018-present, Jiali Lin All rights reserved.

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or (at
## your option) any later version.

## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.

#' Fix to ensure positive definiteness from Danaher et al.
#' Divide each off-diagonal element by sum of absolute values of off-diagonal elements in its row.
#' 
#' @param A a raw matrix.
#' @param denom_factor scale parameter.
#' 
#' @return B a symmetric and positive definite matrix.
#' 
#' @export
FixMatrix <- function(A, denom_factor) {
  p <- nrow(A)
  for (cur_row in 1:p) {
    cur_sum <- sum(abs(A[cur_row, ])) - 1
    if (cur_sum != 0) {
      A[cur_row, ] <- A[cur_row, ]/(denom_factor * cur_sum)
    }
    # Make sure diagonal entries are still 1
    A[cur_row, cur_row] <- 1
  }
  # Final matrix is average of matrix with its transpose
  B <- (A + t(A))/2
  return(B)
}