#' Continuous ranked probability score
#'
#' Process the CRPS score for a given probability prediction compared to an observed value. 
#' The non-exceedence freequency is processed as follows: \code{"F(x) = (r-a)/(n+b)"}
#'
#' @param x vector of length n (number of members / analogues)
#' @param x0 observed value
#' @param a constant that has a law-dependent optimum from which the samples are derived
#' @param b constant that has a law-dependent optimum from which the samples are derived
#'
#' @return Value of the CRPS score
#'
#' @examples
#' data <- read.table('testdaten.txt',header=F)
#' x <- sqrt(data$V5/70.82)
#' x0 <- sqrt(24.4/70.82)
#' result <- CRPS(x, x0, a=0.375, b=0.25)
#'
#' @export

CRPS <- function(x, x0, a=0.44, b=0.12) {
  x.s <- sort(x)
  n <- length(x)
  r <- 1:n
  Fx <- (r - a) / (n + b)
  res <- 0
  # add rectangle on right side if observed value is on the right
  if (x0 > x.s[n]) {
    res <- res + x0 - x.s[n]
  }
  # add rectangle on the left side if observed value is on the left
  if (x0 < x.s[1]) {
    res <- res + x.s[1] - x0
  }
  if (n > 1) {
    for (i in 1:(n - 1)) {
      # if value on left side of observed value
      if (x.s[i] < x0) {
        # case 1: next value also on left side of observed value
        if (x.s[i + 1] <= x0) {
          res <- res + (Fx[i] ^ 2 + Fx[i + 1] ^ 2) * (x.s[i + 1] - x.s[i]) / 2
        }
        # case 2: x0 before next value
        else {
          F0 <- (Fx[i + 1] - Fx[i]) * (x0 - x.s[i]) / (x.s[i + 1] - x.s[i]) + Fx[i]
          res <- res + (F0 ^ 2 + Fx[i] ^ 2) * (x0 - x.s[i]) / 2
          res <- res + ((F0 - 1) ^ 2 + (Fx[i + 1] - 1) ^ 2) * (x.s[i + 1] - x0) / 2
        }
      }
      # if value on the right side of observed value
      else {
        res <- res + ((Fx[i] - 1) ^ 2 + (Fx[i + 1] - 1) ^ 2) * (x.s[i + 1] - x.s[i]) /
          2
      }
    }
  }
  return(res)
} 
