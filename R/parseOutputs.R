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
  path.values <- paste(directory, '/', period, '/AnalogValues_id_', 
                       station.id, '_step_', level-1, '.nc', sep='')
  assertthat::assert_that(file.exists(path.values), 
                          msg = paste(path.values, 'not found'))
  
  path.dates <- paste(directory, '/', period, '/AnalogDates_id_', 
                      station.id, '_step_', level-1, '.nc', sep='')
  assertthat::assert_that(file.exists(path.dates), 
                          msg = paste(path.dates, 'not found'))
  
  path.scores <- paste(directory, '/', period, '/Scores_id_', 
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
  path.dates <- paste(directory, '/', period, '/AnalogDates_id_', 
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
#' Extract results (for both analogues and the target situations: dates) from 
#' the NetCDF files resulting from AtmoSwing optimizer.
#'
#' @param filePath Path to the file containing the analog dates
#'
#' @return Results of the analogue method.
#'
#' @examples
#' \dontrun{
#' data <- atmoswing::parseDatesNcFile('path/to/results/validation/AnalogDates_id_291_step_0.nc')
#' }
#' 
#' @export
#' 
parseDatesNcFile <- function(filePath) {
  
  assertthat::assert_that(file.exists(filePath), 
                          msg = paste(filePath, 'not found.'))
  
  # Open all files
  AD.nc = ncdf4::nc_open(filePath)
  
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
  path.values <- paste(directory, '/', period, '/AnalogValues_id_', 
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

#' Parse NetCDF files resulting from AtmoSwing optimizer.
#'
#' Extract results (for both analogues and the target situations: predictand 
#' values) from the NetCDF files resulting from AtmoSwing optimizer.
#'
#' @param filePath Path to the file containing the analog values
#'
#' @return Results of the analogue method.
#'
#' @examples
#' \dontrun{
#' data <- atmoswing::parseDatesNcFile('path/to/results/validation/AnalogValues_id_291_step_0.nc')
#' }
#' 
#' @export
#' 
parseValuesNcFile <- function(filePath) {
  
  assertthat::assert_that(file.exists(filePath), 
                          msg = paste(filePath, 'not found.'))
  
  # Open all files
  AV.nc = ncdf4::nc_open(filePath)
  
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
  path.values <- paste(directory, '/', period, '/AnalogValues_id_', 
                       station.id, '_step_', level-1, '.nc', sep='')
  assertthat::assert_that(file.exists(path.values), 
                          msg = paste(path.values, 'not found'))
  
  path.scores <- paste(directory, '/', period, '/Scores_id_', 
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

#' Build a dataframe to store results.
#'
#' Build an empty dataframe to store the stations properties.
#'
#'
#' @return The empty dataframe.
#'
#' @examples
#' \dontrun{
#' stations <- atmoswing::createStationsDataframe()
#' }
#' 
#' @export
#' 
createStationsDataframe <- function() {
  
  stations <- data.frame(
    id = ncvar_get(stationsFile, varid = "station_ids"),
    x = ncvar_get(stationsFile, varid = "station_x_coords"),
    y = ncvar_get(stationsFile, varid = "station_y_coords"),
    h = ncvar_get(stationsFile, varid = "station_heights"),
    p10 = ncvar_get(stationsFile, varid = "daily_precipitations_for_return_periods")[4,]
  )
  
  return(stations)
} 

#' Parse parameters resulting from AtmoSwing optimizer.
#'
#' Extract resulting parameters from the text files resulting from AtmoSwing 
#' optimizer.
#'
#' @param directory Root directory of multiple runs.
#' @param datasets List of datasets (must be used as folder names - e.g. /JRA-55/)
#' @param methods List of methods (must be used as folder names - e.g. /4Z/)
#'
#' @return Results of the analogue method.
#'
#' @examples
#' \dontrun{
#' datasets <- c('CFSR', 'ERA-20C', 'JRA-55')
#' methods <- c('2Z', '4Z', '4Z-2MI')
#' data <- atmoswing::parseAllResultsText('path/to/runs', datasets, methods)
#' }
#' 
#' @export
#' 
parseAllResultsText <- function(directory, datasets, methods) {
  
  stations <- atmoswing::createStationsDataframe()
  
  # List run files
  files <- list.files(c(directory, ""), pattern = "_station_(.*)_best_parameters.txt", 
                      full.names = TRUE, recursive = TRUE)
  
  # Parse results for all datasets
  pb <- txtProgressBar(max = length(datasets)*length(methods))
  for (dataset in datasets) {
    filesSlctDat <- files[ grep(paste("/", dataset, "/", sep = ""), files) ]
    for (method in methods) {
      filesSlct <- filesSlctDat[ grep(paste("/", method, "/", sep = ""), files) ]
      fieldNameCalib <- paste(dataset, "_", method, "_calib", sep = "")
      fieldNameValid <- paste(dataset, "_", method, "_valid", sep = "")
      fieldAnb <- paste(dataset, "_", method, "_anb", sep = "")
      fieldXmin <- paste(dataset, "_", method, "_xmin", sep = "")
      fieldYmin <- paste(dataset, "_", method, "_ymin", sep = "")
      fieldXstep <- paste(dataset, "_", method, "_xstep", sep = "")
      fieldYstep <- paste(dataset, "_", method, "_ystep", sep = "")
      fieldXw <- paste(dataset, "_", method, "_xw", sep = "")
      fieldYw <- paste(dataset, "_", method, "_yw", sep = "")
      
      # Parse files
      for (file in filesSlct) {
        dat <- read.delim(file, header = FALSE, skip = 1)
        posStation <- which(dat == "Station")
        stationId <- dat[[posStation+1]]
        
        # Analog number
        pos <- which(dat == "Anb")
        for (i in 1:length(pos)) {
          fullFieldName <- paste(fieldAnb, "_", i, sep = "")
          
          if(!fullFieldName %in% colnames(stations)) {
            stations[fullFieldName] <- NA
          }
          stations[[fullFieldName]][which(stations$id == stationId)] <- dat[[pos[[i]]+1]]
        }
        
        # Xmin
        pos <- which(dat == "Xmin")
        for (i in 1:length(pos)) {
          fullFieldName <- paste(fieldXmin, "_", i, sep = "")
          
          if(!fullFieldName %in% colnames(stations)) {
            stations[fullFieldName] <- NA
          }
          stations[[fullFieldName]][which(stations$id == stationId)] <- dat[[pos[[i]]+1]]
        }
        
        # Ymin
        pos <- which(dat == "Ymin")
        for (i in 1:length(pos)) {
          fullFieldName <- paste(fieldYmin, "_", i, sep = "")
          
          if(!fullFieldName %in% colnames(stations)) {
            stations[fullFieldName] <- NA
          }
          stations[[fullFieldName]][which(stations$id == stationId)] <- dat[[pos[[i]]+1]]
        }
        
        # Xstep
        pos <- which(dat == "Xstep")
        for (i in 1:length(pos)) {
          fullFieldName <- paste(fieldXstep, "_", i, sep = "")
          
          if(!fullFieldName %in% colnames(stations)) {
            stations[fullFieldName] <- NA
          }
          stations[[fullFieldName]][which(stations$id == stationId)] <- dat[[pos[[i]]+1]]
        }
        
        # Ystep
        pos <- which(dat == "Ystep")
        for (i in 1:length(pos)) {
          fullFieldName <- paste(fieldYstep, "_", i, sep = "")
          
          if(!fullFieldName %in% colnames(stations)) {
            stations[fullFieldName] <- NA
          }
          stations[[fullFieldName]][which(stations$id == stationId)] <- dat[[pos[[i]]+1]]
        }
        
        # Xpts
        pos <- which(dat == "Xptsnb")
        for (i in 1:length(pos)) {
          fullFieldName <- paste(fieldXw, "_", i, sep = "")
          
          if(!fullFieldName %in% colnames(stations)) {
            stations[fullFieldName] <- NA
          }
          stations[[fullFieldName]][which(stations$id == stationId)] <- (dat[[pos[[i]]+1]] - 1) * dat[[pos[[i]]+9]]
        }
        
        # Ypts
        pos <- which(dat == "Yptsnb")
        for (i in 1:length(pos)) {
          fullFieldName <- paste(fieldYw, "_", i, sep = "")
          
          if(!fullFieldName %in% colnames(stations)) {
            stations[fullFieldName] <- NA
          }
          stations[[fullFieldName]][which(stations$id == stationId)] <- (dat[[pos[[i]]+1]] - 1) * dat[[pos[[i]]+3]]
        }
        
        # Score
        posScore <- which(dat == "Calib")
        if(!fieldNameCalib %in% colnames(stations)) {
          stations[fieldNameCalib] <- NA
          stations[fieldNameValid] <- NA
        }
        stations[[fieldNameCalib]][which(stations$id == stationId)] <- dat[[posScore+1]]
        stations[[fieldNameValid]][which(stations$id == stationId)] <- dat[[posScore+3]]
      }
      
      setTxtProgressBar(pb, value = getTxtProgressBar(pb)+1)
    }
  }
  
  close(pb)
  
  stations
}