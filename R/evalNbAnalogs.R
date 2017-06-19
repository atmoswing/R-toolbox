#' Process the CRPS for an increasing number of analogues.
#'
#' Process the CRPS for an increasing number of analogues for every day of the 
#' target period.
#'
#' @param A Results of AtmoSwing as parsed by atmoswing::parseAllNcOutputs.
#' @param filter.size Length of "running window", has to be odd.
#'
#' @return Matrices with CRPS scores for an increasing number of analogues.
#'
#' @examples
#' \dontrun{
#' data <- atmoswing::parseAllNcOutputs('optim/1/results', 1, 'validation')
#' res <- atmoswing::crpsNbAnalogs(data, 9)
#' }
#' 
#' @export
#' 
crpsNbAnalogs <- function(A, filter.size = 9) {
  
  # Calculation on a single row
  crpsPerRow <- function(analog.values, target.value) {
    crps <- array(NA, dim = length(analog.values))
    for (i in 1:length(analog.values)) {
      crps[i] <- atmoswing::crps(analog.values[1:i], target.value)
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
    stats::filter(x, rep(1 / n, n), sides = 2)
  }
  
  # Smooth with a running average
  crps.nb.analogs.smoothed <- t(apply(crps.nb.analogs, 1, ma, n = filter.size))
  
  # Pack results
  res <- list(crps = crps.nb.analogs,
              crps.smoothed = crps.nb.analogs.smoothed)
  
  return(res)
}


#' Identifies the optimal number of analogues on the crps matrix.
#'
#' Identifies the optimal number of analogues on the crps matrix resulting from
#' atmoswing::crpsNbAnalogs.
#'
#' @param crps.matrix Results of atmoswing::crpsNbAnalogs.
#'
#' @return Vector of the optimal number of analogues for each target date.
#'
#' @examples
#' \dontrun{
#' data <- atmoswing::parseAllNcOutputs('optim/1/results', 1, 'validation')
#' res.crps <- atmoswing::crpsNbAnalogs(data)
#' res.nb <- atmoswing::optimalNbAnalog(res.crps$crps.smoothed)
#' }
#' 
#' @export
#' 
optimalNbAnalog <- function(crps.matrix) {
  
  nb.min.first <- apply(crps.matrix, 1, which.min)
  nb.min.med <- apply(crps.matrix, 1, function(x) stats::median(which(x==min(x, na.rm = TRUE))))
  nb.min.last <- apply(crps.matrix, 1, function(x) max(which(x==min(x, na.rm = TRUE))))

  # Pack results
  res <- list(nb.min.first = nb.min.first,
              nb.min.med = nb.min.med,
              nb.min.last = nb.min.last)
  
  return(res)
}


#' Assign the optimal number of analogues to every analogue date
#'
#' Assign the optimal number of analogues to every analogue date. This only
#' works on the calibration period, where the optimal number of analogues can
#' be assessed for the non independent dates.
#'
#' @param A Results of AtmoSwing as parsed by atmoswing::parseAllNcOutputs.
#' @param target.a.nb Results of atmoswing::optimalNbAnalog
#'
#' @return Matrix of the optimal number of analogues for each analogue date.
#'
#' @examples
#' \dontrun{
#' data <- atmoswing::parseAllNcOutputs('optim/1/results', 1, 'validation')
#' res.crps <- atmoswing::crpsNbAnalogs(data)
#' target.a.nb <- atmoswing::optimalNbAnalog(res.crps$crps.smoothed)
#' res <- atmoswing::optimalNbAnalogOfAnalogs(data, target.a.nb$nb.min.med)
#' }
#' 
#' @export
#' 
optimalNbAnalogOfAnalogs <- function(A, target.a.nb) {

  analog.a.nb <- matrix(NA, nrow = nrow(A$analog.dates.MJD), ncol = ncol(A$analog.dates.MJD))
  
  for (i in 1:length(A$target.dates.MJD)){
    targ.date <- A$target.dates.MJD[i]
    targ.nb.analog <- target.a.nb[i]
    
    analog.a.nb[A$analog.dates.MJD==targ.date] <- targ.nb.analog
  }
  
  return(analog.a.nb)
}
