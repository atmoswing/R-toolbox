#' Process the CRPS for an increasing number of analogues.
#'
#' Process the CRPS for an increasing number of analogues for every day of the 
#' target period.
#'
#' @param A Results of AtmoSwing as parsed by atmoswing::parseNcOutputs.
#' @param filter.size Length of "running window", has to be odd.
#'
#' @return Matrices with CRPS scores for an increasing number of analogues.
#'
#' @examples
#' \dontrun{
#' data <- atmoswing::parseNcOutputs('optimizer-outputs/1/results', 1, 'validation')
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
#' @return Vector of the optimal number of analogues for each target date..
#'
#' @examples
#' \dontrun{
#' data <- atmoswing::parseNcOutputs('optimizer-outputs/1/results', 1, 'validation')
#' res.crps <- atmoswing::crpsNbAnalogs(data, 9)
#' res.nb <- atmoswing::optimalNbAnalog(res.crps$crps.smoothed, 9)
#' }
#' 
#' @export
#' 
optimalNbAnalog <- function(crps.matrix) {
  
  nb.min.first <- apply(crps.matrix, 1, which.min)
  nb.min.med <- apply(crps.matrix, 1, function(x) median(which(x==min(x, na.rm = TRUE))))
  nb.min.last <- apply(crps.matrix, 1, function(x) max(which(x==min(x, na.rm = TRUE))))

  # Pack results
  res <- list(nb.min.first = nb.min.first,
              nb.min.med = nb.min.med,
              nb.min.last = nb.min.last)
  
  return(res)
}
