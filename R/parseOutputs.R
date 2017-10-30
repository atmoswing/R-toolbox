#' Parse NetCDF files resulting from AtmoSwing optimizer.
#'
#' Extract results (for both analogues and the target situations: dates, 
#' predictand values, prediction score) from the NetCDF files resulting from 
#' AtmoSwing optimizer.
#'
#' @param directory Directory containing the outputs from AtmoSwing (containing 
#'   the "calibration" or "validation" directories).
#' @param station.id ID of the station time series.
#' @param period Either "calibration" or "validation".
#' @param level Analogy level.
#'
#' @return Results of the analogue method.
#'
#' @examples
#' \dontrun{
#' data <- atmoswing::parseAllNcOutputs('optimizer-outputs/1/results', 1, 'validation')
#' }
#' 
#' @export
#' 
parseAllNcOutputs <- function(directory, station.id, period, level = 1) {
  
  assertthat::assert_that((period=='calibration' || period=='validation'), 
                          msg = 'period must be "calibration" or "validation"')
  assertthat::assert_that(assertthat::is.dir(directory), 
                          msg = paste(directory, 'is not a directory (wd:', 
                                      getwd(), ')'))
  
  # Look for the files
  path.values <- paste(directory, '/', period, '/AnalogsValues_id_', 
                       station.id, '_step_', level-1, '.nc', sep='')
  assertthat::assert_that(file.exists(path.values), 
                          msg = paste(path.values, 'not found'))
  path.dates <- paste(directory, '/', period, '/AnalogsDates_id_', 
                      station.id, '_step_', level-1, '.nc', sep='')
  assertthat::assert_that(file.exists(path.dates), 
                          msg = paste(path.dates, 'not found'))
  path.scores <- paste(directory, '/', period, '/AnalogsForecastScores_id_', 
                       station.id, '_step_', level-1, '.nc', sep='')
  assertthat::assert_that(file.exists(path.scores), 
                          msg = paste(path.scores, 'not found'))
  
  # Open all files
  AV.nc = ncdf4::nc_open(path.values)
  AD.nc = ncdf4::nc_open(path.dates)
  AS.nc = ncdf4::nc_open(path.scores)
  
  # Extract data
  AM <- list(
    analog.dates.MJD = t(ncdf4::ncvar_get(AD.nc, 'analog_dates')),
    analog.criteria = t(ncdf4::ncvar_get(AV.nc, 'analog_criteria')),
    analog.values.norm = t(ncdf4::ncvar_get(AV.nc, 'analog_values_norm')),
    analog.values.raw = t(ncdf4::ncvar_get(AV.nc, 'analog_values_gross')),
    target.dates.MJD = ncdf4::ncvar_get(AV.nc, 'target_dates'),
    target.dates.UTC = as.Date(astroFns::dmjd2ut(
      ncdf4::ncvar_get(AV.nc, 'target_dates'), tz= 'UTC' ), format='%Y.%m.%d'),
    target.values.norm = ncdf4::ncvar_get(AV.nc, 'target_values_norm'),
    target.values.raw = ncdf4::ncvar_get(AV.nc, 'target_values_gross'),
    predict.score = t(ncdf4::ncvar_get(AS.nc, 'forecast_scores'))
  )
  
  # Close all files
  ncdf4::nc_close(AV.nc)
  ncdf4::nc_close(AD.nc)
  ncdf4::nc_close(AS.nc)
  
  return(AM)
} 

#' Parse NetCDF files resulting from AtmoSwing optimizer.
#'
#' Extract results (for both analogues and the target situations: dates) from 
#' the NetCDF files resulting from AtmoSwing optimizer.
#'
#' @param directory Directory containing the outputs from AtmoSwing (containing 
#'   the "calibration" or "validation" directories).
#' @param station.id ID of the station time series.
#' @param period Either "calibration" or "validation".
#' @param level Analogy level.
#'
#' @return Results of the analogue method.
#'
#' @examples
#' \dontrun{
#' data <- atmoswing::parseDatesNcOutputs('optimizer-outputs/1/results', 1, 'validation', 2)
#' }
#' 
#' @export
#' 
parseDatesNcOutputs <- function(directory, station.id, period, level = 1) {
  
  assertthat::assert_that((period=='calibration' || period=='validation'), 
                          msg = 'period must be "calibration" or "validation"')
  assertthat::assert_that(assertthat::is.dir(directory), 
                          msg = paste(directory, 'is not a directory (wd:', 
                                      getwd(), ')'))
  
  # Look for the files
  path.dates <- paste(directory, '/', period, '/AnalogsDates_id_', 
                      station.id, '_step_', level-1, '.nc', sep='')
  assertthat::assert_that(file.exists(path.dates), 
                          msg = paste(path.dates, 'not found'))
  
  # Open all files
  AD.nc = ncdf4::nc_open(path.dates)
  
  # Extract data
  AM <- list(
    analog.dates.MJD = t(ncdf4::ncvar_get(AD.nc, 'analog_dates')),
    target.dates.MJD = ncdf4::ncvar_get(AD.nc, 'target_dates'),
    target.dates.UTC = as.Date(astroFns::dmjd2ut(
      ncdf4::ncvar_get(AD.nc, 'target_dates'), tz= 'UTC' ), format='%Y.%m.%d')
  )
  
  # Close all files
  ncdf4::nc_close(AD.nc)
  
  return(AM)
} 

#' Parse NetCDF files resulting from AtmoSwing optimizer.
#'
#' Extract results (for both analogues and the target situations: predictand 
#' values) from the NetCDF files resulting from AtmoSwing optimizer.
#'
#' @param directory Directory containing the outputs from AtmoSwing (containing 
#'   the "calibration" or "validation" directories).
#' @param station.id ID of the station time series.
#' @param period Either "calibration" or "validation".
#' @param level Analogy level.
#'
#' @return Results of the analogue method.
#'
#' @examples
#' \dontrun{
#' data <- atmoswing::parseValuesNcOutputs('optimizer-outputs/1/results', 1, 'validation')
#' }
#' 
#' @export
#' 
parseValuesNcOutputs <- function(directory, station.id, period, level = 1) {
  
  assertthat::assert_that((period=='calibration' || period=='validation'), 
                          msg = 'period must be "calibration" or "validation"')
  assertthat::assert_that(assertthat::is.dir(directory), 
                          msg = paste(directory, 'is not a directory (wd:', 
                                      getwd(), ')'))
  
  # Look for the files
  path.values <- paste(directory, '/', period, '/AnalogsValues_id_', 
                       station.id, '_step_', level-1, '.nc', sep='')
  assertthat::assert_that(file.exists(path.values), 
                          msg = paste(path.values, 'not found'))
  
  # Open all files
  AV.nc = ncdf4::nc_open(path.values)
  
  # Extract data
  AM <- list(
    analog.values.raw = t(ncdf4::ncvar_get(AV.nc, 'analog_values_gross')),
    target.values.raw = ncdf4::ncvar_get(AV.nc, 'target_values_gross'),
    target.dates.MJD = ncdf4::ncvar_get(AV.nc, 'target_dates'),
    target.dates.UTC = as.Date(astroFns::dmjd2ut(
      ncdf4::ncvar_get(AV.nc, 'target_dates'), tz= 'UTC' ), format='%Y.%m.%d')
  )
  
  # Close all files
  ncdf4::nc_close(AV.nc)
  
  return(AM)
} 

#' Parse NetCDF score files resulting from AtmoSwing optimizer.
#'
#' Extract results (for the target situations: predictand values, prediction 
#' score) from the NetCDF files resulting from AtmoSwing optimizer.
#'
#' @param directory Directory containing the outputs from AtmoSwing (containing 
#'   the "calibration" or "validation" directories).
#' @param station.id ID of the station time series.
#' @param period Either "calibration" or "validation".
#' @param level Analogy level.
#'
#' @return Results of the analogue method.
#'
#' @examples
#' \dontrun{
#' data <- atmoswing::parseScoresNcOutputs('optimizer-outputs/1/results', 1, 'validation')
#' }
#' 
#' @export
#' 
parseScoresNcOutputs <- function(directory, station.id, period, level = 1) {
  
  assertthat::assert_that((period=='calibration' || period=='validation'), 
                          msg = 'period must be "calibration" or "validation"')
  assertthat::assert_that(assertthat::is.dir(directory), 
                          msg = paste(directory, 'is not a directory (wd:', 
                                      getwd(), ')'))
  
  # Look for the files
  path.values <- paste(directory, '/', period, '/AnalogsValues_id_', 
                       station.id, '_step_', level-1, '.nc', sep='')
  assertthat::assert_that(file.exists(path.values), 
                          msg = paste(path.values, 'not found'))
  path.scores <- paste(directory, '/', period, '/AnalogsForecastScores_id_', 
                       station.id, '_step_', level-1, '.nc', sep='')
  assertthat::assert_that(file.exists(path.scores), 
                          msg = paste(path.scores, 'not found'))
  
  # Open all files
  AV.nc = ncdf4::nc_open(path.values)
  AS.nc = ncdf4::nc_open(path.scores)
  
  # Extract data
  AM <- list(
    target.dates.MJD = ncdf4::ncvar_get(AV.nc, 'target_dates'),
    target.dates.UTC = as.Date(astroFns::dmjd2ut(
      ncdf4::ncvar_get(AV.nc, 'target_dates'), tz= 'UTC' ), format='%Y.%m.%d'),
    target.values.norm = ncdf4::ncvar_get(AV.nc, 'target_values_norm'),
    target.values.raw = ncdf4::ncvar_get(AV.nc, 'target_values_gross'),
    predict.score = t(ncdf4::ncvar_get(AS.nc, 'forecast_scores'))
  )
  
  # Close all files
  ncdf4::nc_close(AV.nc)
  ncdf4::nc_close(AS.nc)
  
  return(AM)
} 