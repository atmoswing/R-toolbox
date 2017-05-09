library(assertthat)

#' Parse NetCDF files resulting from AtmoSwing optimizer
#'
#' Extract results (for both analogues and the target situations: dates, predictand values, 
#' prediction score) from the NetCDF files resulting from AtmoSwing optimizer
#'
#' @param directory directory containing the outputs from AtmoSwing (containing the "calibration" 
#'     or "validation" directories)
#' @param station.id ID of the station time series
#' @param period either "calibration" or "validation"
#'
#' @return Results of the analogue method
#'
#' @examples
#' data <- atmoswingRToolbox::parseNcOutputs('path/to/dir', 22, 'calibration')
#'
#' @export

parseNcOutputs <- function(directory, station.id, period) {
  assert_that((period=="calibration" || period=="validation"), msg = "period must be calibration' or 'validation'")
  assert_that(is.dir(directory), msg = paste(directory, "is not a directory"))
  
  # Look for the files
  path.values <- paste(directory, '/', period, '/AnalogsValues_id_', station.id, '_step_0.nc', sep='')
  assert_that(file.exists(path.values), msg = paste(path.values, "not found"))
  path.dates <- paste(directory, '/', period, '/AnalogsDates_id_', station.id, '_step_0.nc', sep='')
  assert_that(file.exists(path.dates), msg = paste(path.dates, "not found"))
  path.scores <- paste(directory, '/', period, '/AnalogsForecastScores_id_', station.id, '_step_0.nc', sep='')
  assert_that(file.exists(path.scores), msg = paste(path.scores, "not found"))
  
  # Open all files
  AV.nc = ncdf4::nc_open(path.values)
  AD.nc = ncdf4::nc_open(path.dates)
  AS.nc = ncdf4::nc_open(path.scores)
  
  # Extract data
  
  AM <- list(
    analog.values.norm = t(ncdf4::ncvar_get(AV.nc, "analog_values_norm")),
    analog.values.raw = t(ncdf4::ncvar_get(AV.nc, "analog_values_gross")),
    analog.criteria = t(ncdf4::ncvar_get(AV.nc, "analog_criteria")),
    target.dates = ncdf4::ncvar_get(AV.nc, "target_dates"),
    target.values.norm = ncdf4::ncvar_get(AV.nc, "target_values_norm"),
    target.values.raw = ncdf4::ncvar_get(AV.nc, "target_values_gross"),
    analog.dates = t(ncdf4::ncvar_get(AD.nc, "analog_dates")),
    predict.score = t(ncdf4::ncvar_get(AS.nc, "forecast_scores"))
  )
  
  # Close all files
  ncdf4::nc_close(AV.nc)
  ncdf4::nc_close(AD.nc)
  ncdf4::nc_close(AS.nc)
  
  return(AM)
} 