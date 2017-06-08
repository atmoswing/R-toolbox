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
#' @param w Weight for every analogue day
#'
#' @return Value of the CRPS score.
#'
#' @examples
#' obs <- 25.8
#' precip <- c(35.3, 31.5, 14.9, 30.0, 3.0, 4.2, 15.1, 1.9, 1.2, 2.1)
#' result <- atmoswing::crps(precip, obs, a=0.375, b=0.25)
#' 
#' @export
#' 
crps <- function(x, x0, a=0.44, b=0.12, w=NA) {
  
  n <- length(x)
  
  if (is.na(w)) {
    x.s <- sort(x)
    r <- 1:n
    Fx <- (r - a) / (n + b)
  } else {
    assert_that(length(Fx) == n)
    sorted <- sort.int(x, index.return=TRUE)
    x.s <- sorted$x
    w <- w[sorted$ix]
    fx <- w / sum(w)
    Fx <- cumsum(fx) - fx/2
  }
  
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

#' Process the CRPS for a given number of analogues.
#'
#' Process the CRPS for a given number of analogues for every day of the 
#' target period.
#'
#' @param A Results of AtmoSwing as parsed by atmoswing::parseNcOutputs.
#' @param nb.analogs Number of analogs to consider (all of them if ignored or 0)
#'
#' @return Vector of the CRPS value for every day of the target period.
#'
#' @examples
#' \dontrun{
#' data <- atmoswing::parseNcOutputs('optim/1/results', 1, 'validation')
#' res <- atmoswing::crpsVector(data, 30)
#' }
#' 
#' @export
#' 
crpsVector <- function(A, nb.analogs = 0) {
  
  # Check provided number of analogues
  if (nb.analogs == 0) {
    nb.analogs <- ncol(A$analog.values.norm)
  } else {
    if (nb.analogs > ncol(A$analog.values.norm)) {
      nb.analogs <- ncol(A$analog.values.norm)
      warning("The provided number of analogues is higher than the number of analogues available.")
    }
  }
  
  # Calculation on a single row
  crpsPerRow <- function(analog.values, target.value, nb.analogs) {
    return (atmoswing::crps(analog.values[1:nb.analogs], target.value))
  }
  
  # Apply on the whole matrix
  crps.vect <- t(mapply(
    crpsPerRow,
    split(A$analog.values.norm, row(A$analog.values.norm)),
    A$target.values.norm,
    nb.analogs
  ))
  
  return(crps.vect)
}

#' Process the CRPS weighted by the analogy criteria for the target period..
#'
#' Process the CRPS weighted by the analogy criteria for a given number of 
#' analogues for every day of the target period.
#'
#' @param A Results of AtmoSwing as parsed by atmoswing::parseNcOutputs.
#' @param nb.analogs Number of analogs to consider (all of them if ignored or 0)
#'
#' @return Vector of the CRPS value for every day of the target period.
#'
#' @examples
#' \dontrun{
#' data <- atmoswing::parseNcOutputs('optim/1/results', 1, 'validation')
#' res <- atmoswing::crpsVector(data, 30)
#' }
#' 
#' @export
#' 
crpsVectorWeightedByCriteria <- function(A, nb.analogs = 0) {
  
  # Check provided number of analogues
  if (nb.analogs == 0) {
    nb.analogs <- ncol(A$analog.values.norm)
  } else {
    if (nb.analogs > ncol(A$analog.values.norm)) {
      nb.analogs <- ncol(A$analog.values.norm)
      warning("The provided number of analogues is higher than the number of analogues available.")
    }
  }
  
  # Calculation on a single row
  crpsPerRow <- function(analog.values, target.value, analog.criteria, nb.analogs) {
    weights <- abs(analog.criteria[1:nb.analogs] - analog.criteria[nb.analogs])
    return (atmoswing::crps(analog.values[1:nb.analogs], target.value, a=0.44, b=0.12, w=weights))
  }
  
  # Apply on the whole matrix
  crps.vect <- t(mapply(
    crpsPerRow,
    split(A$analog.values.norm, row(A$analog.values.norm)),
    A$target.values.norm,
    split(A$analog.criteria, row(A$analog.criteria)),
    nb.analogs
  ))
  
  return(crps.vect)
}

#' Get the contribution to the CRPS per predictand classes.
#'
#' Get the contribution to the CRPS for different predictand classes.
#'
#' @param crps.vect A vector of the CRPS value for every target day.
#' @param obs.vect A vector of the observed value (same length as crps.vect).
#'
#' @return Contribution to the CRPS (%).
#'
#' @examples
#' \dontrun{
#' data <- atmoswing::parseNcOutputs('optim/1/results', 1, 'validation')
#' crps.vect <- atmoswing::crpsVector(data, 10)
#' res <- atmoswing::crpsContribClasses(crps.vect, data$target.values.raw)
#' barplot(res$contrib, names.arg=res$labels, las=2, cex.names=0.7)
#' title('Contribution (%) to the CRPS (not cumulative)')
#' }
#' 
#' @export
#' 
crpsContribClasses <- function(crps.vect, obs.vect) {
  
  contrib <- matrix(NA, ncol = 7)
  
  contrib[1] <- sum(crps.vect[obs.vect == 0])
  contrib[2] <- sum(crps.vect[obs.vect > 0
                              & obs.vect <= stats::quantile(obs.vect, probs = 0.8)])
  contrib[3] <- sum(crps.vect[obs.vect > stats::quantile(obs.vect, probs = 0.8)
                              & obs.vect <= stats::quantile(obs.vect, probs = 0.9)])
  contrib[4] <- sum(crps.vect[obs.vect > stats::quantile(obs.vect, probs = 0.9)
                              & obs.vect <= stats::quantile(obs.vect, probs = 0.95)])
  contrib[5] <- sum(crps.vect[obs.vect > stats::quantile(obs.vect, probs = 0.95)
                              & obs.vect <= stats::quantile(obs.vect, probs = 0.97)])
  contrib[6] <- sum(crps.vect[obs.vect > stats::quantile(obs.vect, probs = 0.97)
                              & obs.vect <= stats::quantile(obs.vect, probs = 0.99)])
  contrib[7] <- sum(crps.vect[obs.vect > stats::quantile(obs.vect, probs = 0.99)])
  
  contrib <- 100 * contrib / sum(contrib)
  
  labels <- c("P = 0", "0 < P <= q.8", "q.8 < P <= q.9", "q.9 < P <= q.95",
              "q.95 < P <= q.97", "q.97 < P <= q.99", "P > q.99")
  
  # Pack results
  res <- list(contrib = contrib,
              labels = labels)
                             
  return(res)
}