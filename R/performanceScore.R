#' Continuous ranked probability score.
#'
#' Process the CRPS score for a given probability prediction compared to an 
#' observed value. The non-exceedence freequency is processed as follows: 
#' \code{F(x) = (r-a)/(n+b)}
#'
#' @param x Vector of length n (number of members / analogues).
#' @param x0 Observed value.
#' @param a,b Constants that have a law-dependent optimum from which the samples 
#'   are derived.
#'
#' @return Value of the CRPS score.
#'
#' @examples
#' obs <- 25.8
#' precip <- c(35.3, 31.5, 14.9, 30.0, 3.0, 4.2, 15.1, 1.9, 1.2, 2.1)
#' result <- atmoswingRToolbox::crps(precip, obs, a=0.375, b=0.25)
#' 
#' @export
#' 
crps <- function(x, x0, a=0.44, b=0.12) {
  
  x.s <- sort(x)
  n <- length(x)
  r <- 1:n
  Fx <- (r - a) / (n + b)
  res <- 0
  
  # Add rectangle on right side if observed value is on the right
  if (x0 > x.s[n]) {
    res <- res + x0 - x.s[n]
  }
  
  # Add rectangle on the left side if observed value is on the left
  if (x0 < x.s[1]) {
    res <- res + x.s[1] - x0
  }
  
  if (n > 1) {
    for (i in 1:(n - 1)) {
      # If value on left side of observed value
      if (x.s[i] < x0) {
        # Case 1: next value also on left side of observed value
        if (x.s[i + 1] <= x0) {
          res <- res + (Fx[i] ^ 2 + Fx[i+1] ^ 2) * (x.s[i+1] - x.s[i]) / 2
        }
        # Case 2: x0 before next value
        else {
          F0 <- (Fx[i+1] - Fx[i]) * (x0 - x.s[i]) / (x.s[i+1] - x.s[i]) + Fx[i]
          res <- res + (F0 ^ 2 + Fx[i] ^ 2) * (x0 - x.s[i]) / 2
          res <- res + ((F0 - 1) ^ 2 + (Fx[i+1] - 1) ^ 2) * (x.s[i+1] - x0) / 2
        }
      }
      # If value on the right side of observed value
      else {
        res <- res + ((Fx[i] - 1) ^ 2 + (Fx[i+1] - 1) ^ 2) * (x.s[i+1] - x.s[i]) / 2
      }
    }
  }
  
  return(res)
} 


#' Process the CRPS for an increasing number of analogues.
#'
#' Process the CRPS for an increasing number of analogues for every day of the 
#' target period.
#'
#' @param A Results of AtmoSwing as parsed by atmoswingRToolbox::parseNcOutputs.
#' @param filter.size Length of "running window", has to be odd.
#'
#' @return Matrices with CRPS scores for an increasing number of analogues.
#'
#' @examples
#' data <- atmoswingRToolbox::parseNcOutputs('path/to/dir', 22, 'calibration')
#' res <- atmoswingRToolbox::crpsNbAnalogs(data, 9)
#' 
#' @export
#' 
crpsNbAnalogs <- function(A, filter.size = 9) {
  
  # Calculation on a single row
  crpsPerRow <- function(analog.values, target.value) {
    crps <- array(NA, dim = length(analog.values))
    for (i in 1:length(analog.values)) {
      crps[i] <- atmoswingRToolbox::CRPS(analog.values[1:i], target.value)
    }
    return(crps)
  }
  
  # Apply on the whole matrix
  crps.nb.analogs <- t(mapply(
    crpsPerRow,
    split(A$analog.values.norm, row(A$analog.values.norm)),
    A$target.values.norm
  ))
  
  # Moving average
  ma <- function(x, n = 9) {
    filter(x, rep(1 / n, n), sides = 2)
  }
  
  # Smooth with a running average
  crps.nb.analogs.smoothed <- t(apply(crps.nb.analogs, 1, ma, n = filter.size))
  
  # Pack results
  res <- list(crps = crps.nb.analogs,
              crps.smoothed = crps.nb.analogs.smoothed)
  
  return(res)
} 
