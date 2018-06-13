#' Process the CRPSS of multiple methods / datasets
#'
#' Process the CRPSS of multiple methods / datasets
#'
#' @param directory Root directory of multiple runs.
#' @param predictandDB Path to the predictand DB.
#' @param datasets List of datasets (must be used as folder names - e.g. /JRA-55/)
#' @param methods List of methods (must be used as folder names - e.g. /4Z/)
#' @param period Either "calibration" or "validation".
#'
#' @return Dataframe with CRPSS scores for every method / dataset.
#'
#' @examples
#' \dontrun{
#' datasets <- c('CFSR', 'ERA-20C', 'JRA-55')
#' methods <- c('2Z', '4Z', '4Z-2MI')
#' stations <- atmoswing::getCrpss('path/to/runs', datasets, methods, 'validation')
#' }
#' 
#' @export
#' 
getCrpss <- function(directory, predictandDB, datasets, methods, period) {
  
  stations <- atmoswing::createStationsDataframe(predictandDB)
  
  # List run dirs
  dirs <- list.dirs(c(directory, ""), full.names = TRUE, recursive = TRUE)
  dirs <- dirs[ grep("results", dirs) ]
  dirs <- dirs[ grep("results/", dirs, invert = TRUE) ]
  
  # Parse results for all datasets
  pb <- txtProgressBar(max = length(datasets)*length(methods))
  for (dataset in datasets) {
    dirsSlctDat <- dirs[ grep(paste("/", dataset, "/", sep = ""), dirs) ]
    for (method in methods) {
      dirsSlct <- dirsSlctDat[ grep(paste("/", method, "/", sep = ""), dirsSlctDat) ]
      
      fieldName <- paste(dataset, "_", method, "_crpss", sep = "")
      stations[fieldName] <- NA
      
      # Parse files
      for (dirPath in dirsSlct) {
        for (i_st in 1:nrow(stations)) {
          
          pattern <- paste('AnalogValues_id_', i_st, '_step_.{1}', sep='')
          resFile <- list.files(paste(dirPath, period, sep = "/"), 
                                pattern = pattern, full.names = TRUE, 
                                recursive = TRUE)
          assertthat::assert_that(length(resFile) == 1)
          
          res <- atmoswing::parseValuesNcFile(resFile)
          
          obs <- res$target.values.raw
          ens <- res$analog.values.raw
          
          clim <- t(matrix( obs, length(obs), length(obs)))
          
          fieldName <- paste(dataset, "_", method, "_crpss", sep = "")
          stations[[fieldName]][which(stations$id == i_st)] <- 
            SkillScore(EnsCrps(ens, obs), EnsCrps(clim, obs), 
                       handle.na = "use.pairwise.complete")[1]
          
        }
      }
      
      setTxtProgressBar(pb, value = getTxtProgressBar(pb)+1)
    }
  }
  
  close(pb)
  
  stations
}


#' Process the CRPSS conditionally to predictand thresholds
#'
#' Process the CRPSS of multiple methods / datasets conditionally to predictand thresholds
#'
#' @param directory Root directory of multiple runs.
#' @param predictandDB Path to the predictand DB.
#' @param datasets List of datasets (must be used as folder names - e.g. /JRA-55/)
#' @param methods List of methods (must be used as folder names - e.g. /4Z/)
#' @param period Either "calibration" or "validation".
#' @param startYear First year of the total period.
#' @param endYear Last year of the total period.
#'
#' @return Dataframe with CRPSS scores for every method / dataset.
#'
#' @examples
#' \dontrun{
#' datasets <- c('CFSR', 'ERA-20C', 'JRA-55')
#' methods <- c('2Z', '4Z', '4Z-2MI')
#' stations <- atmoswing::getCrpssThresholds('path/to/runs', 'path/to/db', datasets, methods, 'validation', 1981, 2010)
#' }
#' 
#' @export
#' 
getCrpssThresholds <- function(directory, predictandDB, datasets, methods, period, startYear, endYear) {
  
  stations <- atmoswing::createStationsDataframe(predictandDB)
  
  predictandDB.nc <- ncdf4::nc_open(predictandDB)
  
  thresNames <- c("eq0", "sup0", "0.80", "0.90", "0.95", "0.99")
  
  # Precip thresholds
  precipTime <- ncvar_get(predictandDB.nc, varid = "time")
  precipVal <-  ncvar_get(predictandDB.nc, varid = "data")
  
  ncdf4::nc_close(predictandDB.nc)
  
  timeSlct <- precipTime>=astroFns::ut2dmjd(startYear, 1, 1) & 
    precipTime<=astroFns::ut2dmjd(endYear, 12, 31)
  precipVal <- precipVal[, timeSlct]
  
  q.80 <- apply(precipVal, 1, stats::quantile, probs = 0.8, na.rm = TRUE)
  q.90 <- apply(precipVal, 1, stats::quantile, probs = 0.9, na.rm = TRUE)
  q.95 <- apply(precipVal, 1, stats::quantile, probs = 0.95, na.rm = TRUE)
  q.99 <- apply(precipVal, 1, stats::quantile, probs = 0.99, na.rm = TRUE)
  
  stations <- cbind(q.80, q.90, q.95, q.99)
  
  # List run dirs
  dirs <- list.dirs(c(directory, ""), full.names = TRUE, recursive = TRUE)
  dirs <- dirs[ grep("results", dirs) ]
  dirs <- dirs[ grep("results/", dirs, invert = TRUE) ]
  
  # Parse results for all datasets
  pb <- txtProgressBar(max = length(datasets)*length(methods))
  for (dataset in datasets) {
    dirsSlctDat <- dirs[ grep(paste("/", dataset, "/", sep = ""), dirs) ]
    for (method in methods) {
      dirsSlct <- dirsSlctDat[ grep(paste("/", method, "/", sep = ""), dirsSlctDat) ]
      
      for(thresName in thresNames) {
        fieldName <- paste(dataset, "_", method, "_", thresName, sep = "")
        stations[fieldName] <- NA
      }
      
      # Parse files
      for (dirPath in dirsSlct) {
        for (i_st in 1:nrow(stations)) {
          
          pattern <- paste('AnalogValues_id_', i_st, '_step_.{1}', sep='')
          resFile <- list.files(paste(dirPath, period, sep = "/"), 
                                pattern = pattern, full.names = TRUE, 
                                recursive = TRUE)
          assertthat::assert_that(length(resFile) == 1)
          
          res <- atmoswing::parseValuesNcFile(resFile)
          
          obs <- res$target.values.raw
          ens <- res$analog.values.raw
          
          clim <- t(matrix( obs, length(obs), length(obs)))
          
          fieldName <- paste(dataset, "_", method, "_eq0", sep = "")
          stations[[fieldName]][which(stations$id == i_st)] <- 
            SkillScore(EnsCrps(ens[res$target.values.raw == 0,], 
                               obs[res$target.values.raw == 0]), 
                       EnsCrps(clim[res$target.values.raw == 0,], 
                               obs[res$target.values.raw == 0]), 
                       handle.na = "use.pairwise.complete")[1]
          
          fieldName <- paste(dataset, "_", method, "_sup0", sep = "")
          stations[[fieldName]][which(stations$id == i_st)] <- 
            SkillScore(EnsCrps(ens[res$target.values.raw > 0,], 
                               obs[res$target.values.raw > 0]), 
                       EnsCrps(clim[res$target.values.raw > 0,], 
                               obs[res$target.values.raw > 0]), 
                       handle.na = "use.pairwise.complete")[1]
          
          fieldName <- paste(dataset, "_", method, "_0.80", sep = "")
          stations[[fieldName]][which(stations$id == i_st)] <- 
            SkillScore(EnsCrps(ens[res$target.values.raw >= stations$q.80[i_st],], 
                               obs[res$target.values.raw >= stations$q.80[i_st]]), 
                       EnsCrps(clim[res$target.values.raw >= stations$q.80[i_st],], 
                               obs[res$target.values.raw >= stations$q.80[i_st]]), 
                       handle.na = "use.pairwise.complete")[1]
          
          fieldName <- paste(dataset, "_", method, "_0.90", sep = "")
          stations[[fieldName]][which(stations$id == i_st)] <- 
            SkillScore(EnsCrps(ens[res$target.values.raw >= stations$q.90[i_st],], 
                               obs[res$target.values.raw >= stations$q.90[i_st]]), 
                       EnsCrps(clim[res$target.values.raw >= stations$q.90[i_st],], 
                               obs[res$target.values.raw >= stations$q.90[i_st]]), 
                       handle.na = "use.pairwise.complete")[1]
          
          fieldName <- paste(dataset, "_", method, "_0.95", sep = "")
          stations[[fieldName]][which(stations$id == i_st)] <- 
            SkillScore(EnsCrps(ens[res$target.values.raw >= stations$q.95[i_st],], 
                               obs[res$target.values.raw >= stations$q.95[i_st]]), 
                       EnsCrps(clim[res$target.values.raw >= stations$q.95[i_st],], 
                               obs[res$target.values.raw >= stations$q.95[i_st]]), 
                       handle.na = "use.pairwise.complete")[1]
          
          fieldName <- paste(dataset, "_", method, "_0.99", sep = "")
          stations[[fieldName]][which(stations$id == i_st)] <- 
            SkillScore(EnsCrps(ens[res$target.values.raw >= stations$q.99[i_st],], 
                               obs[res$target.values.raw >= stations$q.99[i_st]]), 
                       EnsCrps(clim[res$target.values.raw >= stations$q.99[i_st],], 
                               obs[res$target.values.raw >= stations$q.99[i_st]]), 
                       handle.na = "use.pairwise.complete")[1]
          
        }
      }
      
      setTxtProgressBar(pb, value = getTxtProgressBar(pb)+1)
    }
  }
  
  close(pb)
  
  stations
}


#' Process the correlation coefficient
#'
#' Process the correlation coefficient of multiple methods / datasets
#'
#' @param directory Root directory of multiple runs.
#' @param predictandDB Path to the predictand DB.
#' @param datasets List of datasets (must be used as folder names - e.g. /JRA-55/)
#' @param methods List of methods (must be used as folder names - e.g. /4Z/)
#' @param period Either "calibration" or "validation".
#'
#' @return Dataframe with correlation coefficients for every method / dataset.
#'
#' @examples
#' \dontrun{
#' datasets <- c('CFSR', 'ERA-20C', 'JRA-55')
#' methods <- c('2Z', '4Z', '4Z-2MI')
#' stations <- atmoswing::getCorrel('path/to/runs', datasets, methods, 'validation')
#' }
#' 
#' @export
#' 
getCorrel <- function(directory, predictandDB, datasets, methods, period) {
  
  stations <- atmoswing::createStationsDataframe(predictandDB)
  
  # List run dirs
  dirs <- list.dirs(c(directory, ""), full.names = TRUE, recursive = TRUE)
  dirs <- dirs[ grep("results", dirs) ]
  dirs <- dirs[ grep("results/", dirs, invert = TRUE) ]
  
  # Parse results for all datasets
  pb <- txtProgressBar(max = length(datasets)*length(methods))
  for (dataset in datasets) {
    dirsSlctDat <- dirs[ grep(paste("/", dataset, "/", sep = ""), dirs) ]
    for (method in methods) {
      dirsSlct <- dirsSlctDat[ grep(paste("/", method, "/", sep = ""), dirsSlctDat) ]
      
      fieldName <- paste(dataset, "_", method, "_corr_1st", sep = "")
      stations[fieldName] <- NA
      fieldName <- paste(dataset, "_", method, "_corr_med", sep = "")
      stations[fieldName] <- NA
      fieldName <- paste(dataset, "_", method, "_corr_mean", sep = "")
      stations[fieldName] <- NA
      
      # Parse files
      for (dirPath in dirsSlct) {
        for (i_st in 1:nrow(stations)) {
          
          pattern <- paste('AnalogValues_id_', i_st, '_step_.{1}', sep='')
          resFile <- list.files(paste(dirPath, period, sep = "/"), 
                                pattern = pattern, full.names = TRUE, 
                                recursive = TRUE)
          assertthat::assert_that(length(resFile) == 1)
          
          res <- atmoswing::parseValuesNcFile(resFile)
          
          med <- apply(res$analog.values.raw, 1, stats::median, na.rm = TRUE)
          mea <- apply(res$analog.values.raw, 1, mean, na.rm = TRUE)
          
          fieldName <- paste(dataset, "_", method, "_corr_1st", sep = "")
          stations[[fieldName]][which(stations$id == i_st)] <- 
            cor(res$target.values.raw, y = res$analog.values.raw[,1], use = "complete.obs")
          
          fieldName <- paste(dataset, "_", method, "_corr_med", sep = "")
          stations[[fieldName]][which(stations$id == i_st)] <- 
            cor(res$target.values.raw, y = med, use = "complete.obs")
          
          fieldName <- paste(dataset, "_", method, "_corr_mean", sep = "")
          stations[[fieldName]][which(stations$id == i_st)] <- 
            cor(res$target.values.raw, y = mea, use = "complete.obs")
          
        }
      }
      
      setTxtProgressBar(pb, value = getTxtProgressBar(pb)+1)
    }
  }
  
  close(pb)
  
  stations
}


#' Process the monthly correlation coefficients
#'
#' Process the monthly correlation coefficients of multiple methods / datasets
#'
#' @param directory Root directory of multiple runs.
#' @param predictandDB Path to the predictand DB.
#' @param datasets List of datasets (must be used as folder names - e.g. /JRA-55/)
#' @param methods List of methods (must be used as folder names - e.g. /4Z/)
#' @param period Either "calibration" or "validation".
#'
#' @return Dataframe with monthly correlation coefficients for every method / dataset.
#'
#' @examples
#' \dontrun{
#' datasets <- c('CFSR', 'ERA-20C', 'JRA-55')
#' methods <- c('2Z', '4Z', '4Z-2MI')
#' stations <- atmoswing::getCorrelMonthly('path/to/runs', datasets, methods, 'validation')
#' }
#' 
#' @export
#' 
getCorrelMonthly <- function(directory, predictandDB, datasets, methods, period) {
  
  stations <- atmoswing::createStationsDataframe(predictandDB)
  
  # List run dirs
  dirs <- list.dirs(c(directory, ""), full.names = TRUE, recursive = TRUE)
  dirs <- dirs[ grep("results", dirs) ]
  dirs <- dirs[ grep("results/", dirs, invert = TRUE) ]
  
  # Parse results for all datasets
  pb <- txtProgressBar(max = length(datasets)*length(methods))
  for (dataset in datasets) {
    dirsSlctDat <- dirs[ grep(paste("/", dataset, "/", sep = ""), dirs) ]
    for (method in methods) {
      dirsSlct <- dirsSlctDat[ grep(paste("/", method, "/", sep = ""), dirsSlctDat) ]
      
      for(currMonth in 1:12) {
        fieldName <- paste(dataset, "_", method, "_corr_mean_", currMonth, sep = "")
        stations[fieldName] <- NA
      }
      
      # Parse files
      for (dirPath in dirsSlct) {
        for (i_st in 1:nrow(stations)) {
          
          pattern <- paste('AnalogValues_id_', i_st, '_step_.{1}', sep='')
          resFile <- list.files(paste(dirPath, period, sep = "/"), 
                                pattern = pattern, full.names = TRUE, 
                                recursive = TRUE)
          assertthat::assert_that(length(resFile) == 1)
          
          res <- atmoswing::parseValuesNcFile(resFile)
          
          mea <- apply(res$analog.values.raw, 1, mean, na.rm = TRUE)
          
          for(currMonth in 1:12) {
            fieldName <- paste(dataset, "_", method, "_corr_mean_", currMonth, sep = "")
            
            stations[[fieldName]][which(stations$id == i_st)] <- 
              cor(res$target.values.raw[month(res$target.dates.UTC) == currMonth], 
                  y = mea[month(res$target.dates.UTC) == currMonth], 
                  use = "complete.obs")
          }
          
        }
      }
      
      setTxtProgressBar(pb, value = getTxtProgressBar(pb)+1)
    }
  }
  
  close(pb)
  
  stations
}


#' Process the annual correlation coefficients
#'
#' Process the annual correlation coefficients of multiple methods / datasets
#'
#' @param directory Root directory of multiple runs.
#' @param predictandDB Path to the predictand DB.
#' @param datasets List of datasets (must be used as folder names - e.g. /JRA-55/)
#' @param methods List of methods (must be used as folder names - e.g. /4Z/)
#' @param startYear First year of the total period.
#' @param endYear Last year of the total period.
#'
#' @return Dataframe with annual correlation coefficients for every method / dataset.
#'
#' @examples
#' \dontrun{
#' datasets <- c('CFSR', 'ERA-20C', 'JRA-55')
#' methods <- c('2Z', '4Z', '4Z-2MI')
#' stations <- atmoswing::getCorrelYearly('path/to/runs', datasets, methods, 1981, 2010)
#' }
#' 
#' @export
#' 
getCorrelYearly <- function(directory, predictandDB, datasets, methods, startYear, endYear) {
  
  stations <- atmoswing::createStationsDataframe(predictandDB)
  
  # List run dirs
  dirs <- list.dirs(c(directory, ""), full.names = TRUE, recursive = TRUE)
  dirs <- dirs[ grep("results", dirs) ]
  dirs <- dirs[ grep("results/", dirs, invert = TRUE) ]
  
  # Parse results for all datasets
  pb <- txtProgressBar(max = length(datasets)*length(methods))
  for (dataset in datasets) {
    dirsSlctDat <- dirs[ grep(paste("/", dataset, "/", sep = ""), dirs) ]
    for (method in methods) {
      dirsSlct <- dirsSlctDat[ grep(paste("/", method, "/", sep = ""), dirsSlctDat) ]
      
      for(y in startYear:endYear) {
        fieldName <- paste(dataset, "_", method, "_corr_mean_", y, sep = "")
        stations[fieldName] <- NA
      }
      
      # Parse files
      for (dirPath in dirsSlct) {
        for (i_st in 1:nrow(stations)) {
          
          pattern <- paste('AnalogValues_id_', i_st, '_step_.{1}', sep='')
          resFileCP <- list.files(paste(dirPath, "calibration", sep = "/"), 
                                  pattern = pattern, full.names = TRUE, 
                                  recursive = TRUE)
          assertthat::assert_that(length(resFileCP) == 1)
          resCP <- atmoswing::parseValuesNcFile(resFileCP)
          
          pattern <- paste('AnalogValues_id_', i_st, '_step_.{1}', sep='')
          resFileVP <- list.files(paste(dirPath, "validation", sep = "/"), 
                                  pattern = pattern, full.names = TRUE, 
                                  recursive = TRUE)
          assertthat::assert_that(length(resFileVP) == 1)
          resVP <- atmoswing::parseValuesNcFile(resFileVP)
          
          for(y in startYear:endYear) {
            fieldName <- paste(dataset, "_", method, "_corr_mean_", y, sep = "")
            
            if (sum(year(resVP$target.dates.UTC) == y) > 0) {
              meaVP <- apply(resVP$analog.values.raw[year(resVP$target.dates.UTC) == y,], 1, mean, na.rm = TRUE)
              stations[[fieldName]][which(stations$id == i_st)] <- 
                cor(resVP$target.values.raw[year(resVP$target.dates.UTC) == y], 
                    y = meaVP, 
                    use = "complete.obs")
            } else {
              meaCP <- apply(resCP$analog.values.raw[year(resCP$target.dates.UTC) == y,], 1, mean, na.rm = TRUE)
              stations[[fieldName]][which(stations$id == i_st)] <- 
                cor(resCP$target.values.raw[year(resCP$target.dates.UTC) == y], 
                    y = meaCP, 
                    use = "complete.obs")
            }
            
          }
          
        }
      }
      
      setTxtProgressBar(pb, value = getTxtProgressBar(pb)+1)
    }
  }
  
  close(pb)
  
  stations
}


#' Process the inter-annual correlation coefficients
#'
#' Process the inter-annual correlation coefficients of multiple methods / datasets
#'
#' @param directory Root directory of multiple runs.
#' @param predictandDB Path to the predictand DB.
#' @param datasets List of datasets (must be used as folder names - e.g. /JRA-55/)
#' @param methods List of methods (must be used as folder names - e.g. /4Z/)
#' @param startYear First year of the total period.
#' @param endYear Last year of the total period.
#'
#' @return Dataframe with inter-annual correlation coefficients for every method / dataset.
#'
#' @examples
#' \dontrun{
#' datasets <- c('CFSR', 'ERA-20C', 'JRA-55')
#' methods <- c('2Z', '4Z', '4Z-2MI')
#' stations <- atmoswing::getInterannualCorrel('path/to/runs', datasets, methods, 1981, 2010)
#' }
#' 
#' @export
#' 
getInterannualCorrel <- function(directory, predictandDB, datasets, methods, startYear, endYear) {
  
  stations <- atmoswing::createStationsDataframe(predictandDB)
  
  # List run dirs
  dirs <- list.dirs(c(directory, ""), full.names = TRUE, recursive = TRUE)
  dirs <- dirs[ grep("results", dirs) ]
  dirs <- dirs[ grep("results/", dirs, invert = TRUE) ]
  
  # Parse results for all datasets
  pb <- txtProgressBar(max = length(datasets)*length(methods))
  for (dataset in datasets) {
    dirsSlctDat <- dirs[ grep(paste("/", dataset, "/", sep = ""), dirs) ]
    for (method in methods) {
      dirsSlct <- dirsSlctDat[ grep(paste("/", method, "/", sep = ""), dirsSlctDat) ]
      
      fieldName <- paste(dataset, "_", method, "_corr_1st", sep = "")
      stations[fieldName] <- NA
      fieldName <- paste(dataset, "_", method, "_corr_med", sep = "")
      stations[fieldName] <- NA
      fieldName <- paste(dataset, "_", method, "_corr_mean", sep = "")
      stations[fieldName] <- NA
      
      # Parse files
      for (dirPath in dirsSlct) {
        for (i_st in 1:nrow(stations)) {
          
          pattern <- paste('AnalogValues_id_', i_st, '_step_.{1}', sep='')
          resFileCP <- list.files(paste(dirPath, "calibration", sep = "/"), 
                                  pattern = pattern, full.names = TRUE, 
                                  recursive = TRUE)
          assertthat::assert_that(length(resFileCP) == 1)
          resCP <- atmoswing::parseValuesNcFile(resFileCP)
          
          pattern <- paste('AnalogValues_id_', i_st, '_step_.{1}', sep='')
          resFileVP <- list.files(paste(dirPath, "validation", sep = "/"), 
                                  pattern = pattern, full.names = TRUE, 
                                  recursive = TRUE)
          assertthat::assert_that(length(resFileVP) == 1)
          resVP <- atmoswing::parseValuesNcFile(resFileVP)
          
          years <- startYear:endYear
          vectorRef <- numeric(length(years))*NA
          vectorMed <- numeric(length(years))*NA
          vectorMean <- numeric(length(years))*NA
          vector1st <- numeric(length(years))*NA
          
          for(y in years) {
            ind <- y-years[1]+1
            
            vectorRef[ind] <- sum(resCP$target.values.raw[year(resCP$target.dates.UTC) == y], na.rm = TRUE)
            vectorMed[ind] <- sum(apply(resCP$analog.values.raw[year(resCP$target.dates.UTC) == y,], 1, stats::median, na.rm = TRUE), na.rm = TRUE)
            vectorMean[ind] <- sum(apply(resCP$analog.values.raw[year(resCP$target.dates.UTC) == y,], 1, mean, na.rm = TRUE), na.rm = TRUE)
            vector1st[ind] <- sum(resCP$analog.values.raw[year(resCP$target.dates.UTC) == y,1], na.rm = TRUE)
            if (vectorMed[ind] == 0) {
              vectorRef[ind] <- sum(resVP$target.values.raw[year(resVP$target.dates.UTC) == y], na.rm = TRUE)
              vectorMed[ind] <- sum(apply(resVP$analog.values.raw[year(resVP$target.dates.UTC) == y,], 1, stats::median, na.rm = TRUE), na.rm = TRUE)
              vectorMean[ind] <- sum(apply(resVP$analog.values.raw[year(resVP$target.dates.UTC) == y,], 1, mean, na.rm = TRUE), na.rm = TRUE)
              vector1st[ind] <- sum(resVP$analog.values.raw[year(resVP$target.dates.UTC) == y,1], na.rm = TRUE)
            }
          }
          
          fieldName <- paste(dataset, "_", method, "_corr_1st", sep = "")
          stations[[fieldName]][which(stations$id == i_st)] <- 
            cor(vectorRef, y = vector1st, use = "complete.obs")
          
          fieldName <- paste(dataset, "_", method, "_corr_med", sep = "")
          stations[[fieldName]][which(stations$id == i_st)] <- 
            cor(vectorRef, y = vectorMed, use = "complete.obs")
          
          fieldName <- paste(dataset, "_", method, "_corr_mean", sep = "")
          stations[[fieldName]][which(stations$id == i_st)] <- 
            cor(vectorRef, y = vectorMean, use = "complete.obs")
          
        }
      }
      
      setTxtProgressBar(pb, value = getTxtProgressBar(pb)+1)
    }
  }
  
  close(pb)
  
  stations
}


#' Process the precipitation annual totals
#'
#' Process the precipitation annual totals of multiple methods / datasets
#'
#' @param directory Root directory of multiple runs.
#' @param predictandDB Path to the predictand DB.
#' @param datasets List of datasets (must be used as folder names - e.g. /JRA-55/)
#' @param methods List of methods (must be used as folder names - e.g. /4Z/)
#' @param startYear First year of the total period.
#' @param endYear Last year of the total period.
#'
#' @return Dataframe with precip annual totals for every method / dataset.
#'
#' @examples
#' \dontrun{
#' datasets <- c('CFSR', 'ERA-20C', 'JRA-55')
#' methods <- c('2Z', '4Z', '4Z-2MI')
#' stations <- atmoswing::getAnnualTotals('path/to/runs', datasets, methods, 1981, 2010)
#' }
#' 
#' @export
#' 
getAnnualTotals <- function(directory, predictandDB, datasets, methods, startYear, endYear) {
  
  stations <- atmoswing::createStationsDataframe(predictandDB)
  
  totals <- list()
  
  # List run dirs
  dirs <- list.dirs(c(directory, ""), full.names = TRUE, recursive = TRUE)
  dirs <- dirs[ grep("results", dirs) ]
  dirs <- dirs[ grep("results/", dirs, invert = TRUE) ]
  
  # Parse results for all datasets
  pb <- txtProgressBar(max = length(datasets)*length(methods))
  for (dataset in datasets) {
    dirsSlctDat <- dirs[ grep(paste("/", dataset, "/", sep = ""), dirs) ]
    for (method in methods) {
      dirsSlct <- dirsSlctDat[ grep(paste("/", method, "/", sep = ""), dirsSlctDat) ]
      
      # Parse files
      for (dirPath in dirsSlct) {
        for (i_st in 1:nrow(stations)) {
          
          pattern <- paste('AnalogValues_id_', i_st, '_step_.{1}', sep='')
          resFileCP <- list.files(paste(dirPath, "calibration", sep = "/"), 
                                  pattern = pattern, full.names = TRUE, 
                                  recursive = TRUE)
          assertthat::assert_that(length(resFileCP) == 1)
          resCP <- atmoswing::parseValuesNcFile(resFileCP)
          
          pattern <- paste('AnalogValues_id_', i_st, '_step_.{1}', sep='')
          resFileVP <- list.files(paste(dirPath, "validation", sep = "/"), 
                                  pattern = pattern, full.names = TRUE, 
                                  recursive = TRUE)
          assertthat::assert_that(length(resFileVP) == 1)
          resVP <- atmoswing::parseValuesNcFile(resFileVP)
          
          years <- startYear:endYear
          vectorRef <- numeric(length(years))*NA
          vectorMed <- numeric(length(years))*NA
          vectorMean <- numeric(length(years))*NA
          vector1st <- numeric(length(years))*NA
          
          for(y in years) {
            ind <- y-years[1]+1
            
            vectorRef[ind] <- sum(resCP$target.values.raw[year(resCP$target.dates.UTC) == y], na.rm = TRUE)
            vectorMed[ind] <- sum(apply(resCP$analog.values.raw[year(resCP$target.dates.UTC) == y,], 1, stats::median, na.rm = TRUE), na.rm = TRUE)
            vectorMean[ind] <- sum(apply(resCP$analog.values.raw[year(resCP$target.dates.UTC) == y,], 1, mean, na.rm = TRUE), na.rm = TRUE)
            vector1st[ind] <- sum(resCP$analog.values.raw[year(resCP$target.dates.UTC) == y,1], na.rm = TRUE)
            if (vectorMed[ind] == 0) {
              vectorRef[ind] <- sum(resVP$target.values.raw[year(resVP$target.dates.UTC) == y], na.rm = TRUE)
              vectorMed[ind] <- sum(apply(resVP$analog.values.raw[year(resVP$target.dates.UTC) == y,], 1, stats::median, na.rm = TRUE), na.rm = TRUE)
              vectorMean[ind] <- sum(apply(resVP$analog.values.raw[year(resVP$target.dates.UTC) == y,], 1, mean, na.rm = TRUE), na.rm = TRUE)
              vector1st[ind] <- sum(resVP$analog.values.raw[year(resVP$target.dates.UTC) == y,1], na.rm = TRUE)
            }
          }
          
          fieldName <- paste(dataset, "_", method, "_", i_st, "_ref", sep = "")
          totals[[fieldName]] <- vectorRef
          
          fieldName <- paste(dataset, "_", method, "_", i_st, "_1st", sep = "")
          totals[[fieldName]] <- vector1st
          
          fieldName <- paste(dataset, "_", method, "_", i_st, "_med", sep = "")
          totals[[fieldName]] <- vectorMed
          
          fieldName <- paste(dataset, "_", method, "_", i_st, "_mean", sep = "")
          totals[[fieldName]] <- vectorMean
          
        }
      }
      
      setTxtProgressBar(pb, value = getTxtProgressBar(pb)+1)
    }
  }
  
  close(pb)
  
  totals
}


#' Process the precipitation monthly totals
#'
#' Process the precipitation monthly totals of multiple methods / datasets
#'
#' @param directory Root directory of multiple runs.
#' @param predictandDB Path to the predictand DB.
#' @param datasets List of datasets (must be used as folder names - e.g. /JRA-55/)
#' @param methods List of methods (must be used as folder names - e.g. /4Z/)
#' @param startYear First year of the total period.
#' @param endYear Last year of the total period.
#'
#' @return Dataframe with precip monthly totals for every method / dataset.
#'
#' @examples
#' \dontrun{
#' datasets <- c('CFSR', 'ERA-20C', 'JRA-55')
#' methods <- c('2Z', '4Z', '4Z-2MI')
#' stations <- atmoswing::getMonthlyTotals('path/to/runs', datasets, methods, 1981, 2010)
#' }
#' 
#' @export
#' 
getMonthlyTotals <- function(directory, predictandDB, datasets, methods, startYear, endYear) {
  
  stations <- atmoswing::createStationsDataframe(predictandDB)
  
  totals <- list()
  
  # List run dirs
  dirs <- list.dirs(c(directory, ""), full.names = TRUE, recursive = TRUE)
  dirs <- dirs[ grep("results", dirs) ]
  dirs <- dirs[ grep("results/", dirs, invert = TRUE) ]
  
  # Parse results for all datasets
  pb <- txtProgressBar(max = length(datasets)*length(methods))
  for (dataset in datasets) {
    dirsSlctDat <- dirs[ grep(paste("/", dataset, "/", sep = ""), dirs) ]
    for (method in methods) {
      dirsSlct <- dirsSlctDat[ grep(paste("/", method, "/", sep = ""), dirsSlctDat) ]
      
      # Parse files
      for (dirPath in dirsSlct) {
        for (i_st in 1:nrow(stations)) {
          
          pattern <- paste('AnalogValues_id_', i_st, '_step_.{1}', sep='')
          resFileCP <- list.files(paste(dirPath, "calibration", sep = "/"), 
                                  pattern = pattern, full.names = TRUE, 
                                  recursive = TRUE)
          assertthat::assert_that(length(resFileCP) == 1)
          resCP <- atmoswing::parseValuesNcFile(resFileCP)
          
          pattern <- paste('AnalogValues_id_', i_st, '_step_.{1}', sep='')
          resFileVP <- list.files(paste(dirPath, "validation", sep = "/"), 
                                  pattern = pattern, full.names = TRUE, 
                                  recursive = TRUE)
          assertthat::assert_that(length(resFileVP) == 1)
          resVP <- atmoswing::parseValuesNcFile(resFileVP)
          
          years <- startYear:endYear
          vectorRef <- numeric(length(12*years))*NA
          vectorMed <- numeric(length(12*years))*NA
          vectorMean <- numeric(length(12*years))*NA
          vector1st <- numeric(length(12*years))*NA
          
          for(y in years) {
            ind <- 12*(y-years[1])
            
            if (sum(resCP$target.values.raw[year(resCP$target.dates.UTC) == y], na.rm = TRUE) == 0) {
              for(m in 1:12) {
                ind2 <- ind+m
                vectorRef[ind2] <- sum(resVP$target.values.raw[year(resVP$target.dates.UTC) == y & month(resVP$target.dates.UTC) == m], na.rm = TRUE)
                vectorMed[ind2] <- sum(apply(resVP$analog.values.raw[year(resVP$target.dates.UTC) == y & month(resVP$target.dates.UTC) == m,], 1, stats::median, na.rm = TRUE), na.rm = TRUE)
                vectorMean[ind2] <- sum(apply(resVP$analog.values.raw[year(resVP$target.dates.UTC) == y & month(resVP$target.dates.UTC) == m,], 1, mean, na.rm = TRUE), na.rm = TRUE)
                vector1st[ind2] <- sum(resVP$analog.values.raw[year(resVP$target.dates.UTC) == y & month(resVP$target.dates.UTC) == m,1], na.rm = TRUE)
              }
            } else {
              for(m in 1:12) {
                ind2 <- ind+m
                vectorRef[ind2] <- sum(resCP$target.values.raw[year(resCP$target.dates.UTC) == y & month(resCP$target.dates.UTC) == m], na.rm = TRUE)
                vectorMed[ind2] <- sum(apply(resCP$analog.values.raw[year(resCP$target.dates.UTC) == y & month(resCP$target.dates.UTC) == m,], 1, stats::median, na.rm = TRUE), na.rm = TRUE)
                vectorMean[ind2] <- sum(apply(resCP$analog.values.raw[year(resCP$target.dates.UTC) == y & month(resCP$target.dates.UTC) == m,], 1, mean, na.rm = TRUE), na.rm = TRUE)
                vector1st[ind2] <- sum(resCP$analog.values.raw[year(resCP$target.dates.UTC) == y & month(resCP$target.dates.UTC) == m,1], na.rm = TRUE)
              }
            }
          }
          
          fieldName <- paste(dataset, "_", method, "_", i_st, "_ref", sep = "")
          totals[[fieldName]] <- vectorRef
          
          fieldName <- paste(dataset, "_", method, "_", i_st, "_1st", sep = "")
          totals[[fieldName]] <- vector1st
          
          fieldName <- paste(dataset, "_", method, "_", i_st, "_med", sep = "")
          totals[[fieldName]] <- vectorMed
          
          fieldName <- paste(dataset, "_", method, "_", i_st, "_mean", sep = "")
          totals[[fieldName]] <- vectorMean
          
        }
      }
      
      setTxtProgressBar(pb, value = getTxtProgressBar(pb)+1)
    }
  }
  
  close(pb)
  
  totals
}


#' Process the bias of multiple methods / datasets
#'
#' Process the bias of multiple methods / datasets
#'
#' @param directory Root directory of multiple runs.
#' @param predictandDB Path to the predictand DB.
#' @param datasets List of datasets (must be used as folder names - e.g. /JRA-55/)
#' @param methods List of methods (must be used as folder names - e.g. /4Z/)
#' @param startYear First year of the total period.
#' @param endYear Last year of the total period.
#'
#' @return Dataframe with the bias for every method / dataset.
#'
#' @examples
#' \dontrun{
#' datasets <- c('CFSR', 'ERA-20C', 'JRA-55')
#' methods <- c('2Z', '4Z', '4Z-2MI')
#' stations <- atmoswing::getBias('path/to/runs', datasets, methods, 1981, 2010)
#' }
#' 
#' @export
#' 
getBias <- function(directory, predictandDB, datasets, methods, startYear, endYear) {
  
  stations <- atmoswing::createStationsDataframe(predictandDB)
  
  # List run dirs
  dirs <- list.dirs(c(directory, ""), full.names = TRUE, recursive = TRUE)
  dirs <- dirs[ grep("results", dirs) ]
  dirs <- dirs[ grep("results/", dirs, invert = TRUE) ]
  
  # Parse results for all datasets
  pb <- txtProgressBar(max = length(datasets)*length(methods))
  for (dataset in datasets) {
    dirsSlctDat <- dirs[ grep(paste("/", dataset, "/", sep = ""), dirs) ]
    for (method in methods) {
      dirsSlct <- dirsSlctDat[ grep(paste("/", method, "/", sep = ""), dirsSlctDat) ]
      
      fieldName <- paste(dataset, "_", method, "_bias_1st", sep = "")
      stations[fieldName] <- NA
      fieldName <- paste(dataset, "_", method, "_bias_med", sep = "")
      stations[fieldName] <- NA
      fieldName <- paste(dataset, "_", method, "_bias_mean", sep = "")
      stations[fieldName] <- NA
      
      # Parse files
      for (dirPath in dirsSlct) {
        for (i_st in 1:nrow(stations)) {
          
          pattern <- paste('AnalogValues_id_', i_st, '_step_.{1}', sep='')
          resFileCP <- list.files(paste(dirPath, "calibration", sep = "/"), 
                                  pattern = pattern, full.names = TRUE, 
                                  recursive = TRUE)
          assertthat::assert_that(length(resFileCP) == 1)
          resCP <- atmoswing::parseValuesNcFile(resFileCP)
          
          pattern <- paste('AnalogValues_id_', i_st, '_step_.{1}', sep='')
          resFileVP <- list.files(paste(dirPath, "validation", sep = "/"), 
                                  pattern = pattern, full.names = TRUE, 
                                  recursive = TRUE)
          assertthat::assert_that(length(resFileVP) == 1)
          resVP <- atmoswing::parseValuesNcFile(resFileVP)
          
          years <- startYear:endYear
          vectorRef <- numeric(length(years))*NA
          vectorMed <- numeric(length(years))*NA
          vectorMean <- numeric(length(years))*NA
          vector1st <- numeric(length(years))*NA
          
          for(y in years) {
            ind <- y-years[1]+1
            
            vectorRef[ind] <- sum(resCP$target.values.raw[year(resCP$target.dates.UTC) == y], na.rm = TRUE)
            vectorMed[ind] <- sum(apply(resCP$analog.values.raw[year(resCP$target.dates.UTC) == y,], 1, stats::median, na.rm = TRUE), na.rm = TRUE)
            vectorMean[ind] <- sum(apply(resCP$analog.values.raw[year(resCP$target.dates.UTC) == y,], 1, mean, na.rm = TRUE), na.rm = TRUE)
            vector1st[ind] <- sum(resCP$analog.values.raw[year(resCP$target.dates.UTC) == y,1], na.rm = TRUE)
            if (vectorMed[ind] == 0) {
              vectorRef[ind] <- sum(resVP$target.values.raw[year(resVP$target.dates.UTC) == y], na.rm = TRUE)
              vectorMed[ind] <- sum(apply(resVP$analog.values.raw[year(resVP$target.dates.UTC) == y,], 1, stats::median, na.rm = TRUE), na.rm = TRUE)
              vectorMean[ind] <- sum(apply(resVP$analog.values.raw[year(resVP$target.dates.UTC) == y,], 1, mean, na.rm = TRUE), na.rm = TRUE)
              vector1st[ind] <- sum(resVP$analog.values.raw[year(resVP$target.dates.UTC) == y,1], na.rm = TRUE)
            }
          }
          
          fieldName <- paste(dataset, "_", method, "_bias_1st", sep = "")
          stations[[fieldName]][which(stations$id == i_st)] <- 
            (mean(vector1st) - mean(vectorRef))/mean(vectorRef)*100
          
          fieldName <- paste(dataset, "_", method, "_bias_med", sep = "")
          stations[[fieldName]][which(stations$id == i_st)] <- 
            (mean(vectorMed) - mean(vectorRef))/mean(vectorRef)*100
          
          fieldName <- paste(dataset, "_", method, "_bias_mean", sep = "")
          stations[[fieldName]][which(stations$id == i_st)] <- 
            (mean(vectorMean) - mean(vectorRef))/mean(vectorRef)*100
          
        }
      }
      
      setTxtProgressBar(pb, value = getTxtProgressBar(pb)+1)
    }
  }
  
  close(pb)
  
  stations
}


#' Get the percentage of similar analogue dates between reanalyses
#'
#' Get the percentage of similar analogue dates between reanalyses per method
#'
#' @param directory Root directory of multiple runs.
#' @param predictandDB Path to the predictand DB.
#' @param datasets List of datasets (must be used as folder names - e.g. /JRA-55/)
#' @param methods List of methods (must be used as folder names - e.g. /4Z/)
#' @param period Either "calibration" or "validation".
#' @param q.threshold Optional precipitation threshold (for target dates). It is expressed as a quantile of days with precipitation (put -1 to ignore)
#'
#' @return Dataframe with the percentage of similar analogue dates
#'
#' @examples
#' \dontrun{
#' datasets <- c('CFSR', 'ERA-20C', 'JRA-55')
#' methods <- c('2Z', '4Z', '4Z-2MI')
#' stations <- atmoswing::getBias('path/to/runs', 'path/to/db', datasets, methods, 'validation', 0.99)
#' }
#' 
#' @export
#' 
getPercentageSimilarDates <- function(directory, predictandDB, datasets, methods, period, q.threshold) {
  
  stations <- atmoswing::createStationsDataframe(predictandDB)
  
  # Initilize cluster
  cl <- makeCluster(detectCores())
  
  # List run dirs
  dirs <- list.dirs(c(dir.data, ""), full.names = TRUE, recursive = TRUE)
  dirs <- dirs[ grep("results", dirs) ]
  dirs <- dirs[ grep("results/", dirs, invert = TRUE) ]
  
  
  # Functions
  compareDates <- function(dataset2, stations, dirs, datasets, dataset1, methods, method, period, q.threshold) {
    
    comparePerStation <- function(i_st, dirsSlct1, dirsSlct2, period, method, q.threshold) {
      
      # Parse analog dates files
      pattern <- paste('AnalogDates_id_', i_st, '_step_.{1}', sep='')
      
      resFile1 <- list.files(paste(dirsSlct1[1], period, sep = "/"), 
                             pattern = pattern, full.names = TRUE, 
                             recursive = TRUE)
      assertthat::assert_that(length(resFile1) == 1)
      res1 <- atmoswing::parseDatesNcFile(resFile1)
      
      resFile2 <- list.files(paste(dirsSlct2[1], period, sep = "/"), 
                             pattern = pattern, full.names = TRUE, 
                             recursive = TRUE)
      assertthat::assert_that(length(resFile2) == 1)
      res2 <- atmoswing::parseDatesNcFile(resFile2)
      
      # Parse analog values files
      if (q.threshold >= 0) {
        pattern <- paste('AnalogValues_id_', i_st, '_step_.{1}', sep='')
        
        prcFile <- list.files(paste(dirsSlct1[1], period, sep = "/"), 
                              pattern = pattern, full.names = TRUE, 
                              recursive = TRUE)
        assertthat::assert_that(length(prcFile) == 1)
        prc <- atmoswing::parseValuesNcFile(prcFile)
      }
      
      if (q.threshold == 0) {
        threshold <- 0
      } else if (q.threshold < 0) {
        threshold <- -99
      } else {
        threshold <- quantile(prc$target.values.raw[prc$target.values.raw>0.1], q.threshold)[[1]]
      }
      
      counter <- 0
      tot <- 0
      for (i_target in seq(1,nrow(res1$analog.dates.MJD))) {
        if (threshold < 0 || prc$target.values.raw[i_target] > threshold) {
          analogs1 <- res1$analog.dates.MJD[i_target,]
          analogs2 <- res2$analog.dates.MJD[i_target,]
          tot <- tot + length(analogs1)
          for(analog in analogs1) {
            if(match(analog, analogs2, nomatch = 0) > 0) {
              counter <- counter+1
            }
          }
        }
      }
      
      counter / tot
    }
    
    dirMeth <- dirs[ grep(paste("/", method, '/', sep = ""), dirs) ]
    
    dirsSlct1 <- dirMeth[ grep(paste("/", dataset1, '/', sep = ""), dirMeth) ]
    dirsSlct2 <- dirMeth[ grep(paste("/", dataset2, '/', sep = ""), dirMeth) ]
    
    assertthat::assert_that(length(dirsSlct1) == 1)
    assertthat::assert_that(length(dirsSlct2) == 1)
    
    simiPerStation <- sapply(1:nrow(stations), comparePerStation, dirsSlct1, dirsSlct2, period, method, q.threshold)
    
    mean(simiPerStation)
  }
  
  results <- array(data = NA, dim=c(length(methods), length(datasets), length(datasets)), dimnames = list(methods, datasets, datasets))
  
  # Parse files
  i_m <- 0
  pb <- txtProgressBar(max = length(datasets)*length(methods))
  for (method in methods) {
    i_m <- i_m + 1
    i_col <- 0
    for (dataset1 in datasets) {
      i_col <- i_col + 1
      
      res <- parLapply(cl, datasets, compareDates, stations, dirs, datasets, dataset1, methods, method, period, q.threshold)
      results[i_m, , i_col] <- unlist(res)
      
      setTxtProgressBar(pb, value = getTxtProgressBar(pb)+1)
    }
  }
  
  stopCluster(cl)
  close(pb)
  
  results <- results*100
}