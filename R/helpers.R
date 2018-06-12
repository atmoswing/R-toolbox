#' Build a dataframe with the predefined color per reanalysis
#'
#' Build a dataframe with the predefined color per reanalysis based on "The 
#' S-RIP colour definitions for reanalysis datasets" (Updated on 11 July 2017, 
#' for the S-RIP Full Report 2018)
#' https://s-rip.ees.hokudai.ac.jp/mediawiki/index.php/Notes_for_Authors
#'
#' @return The colors dataframe.
#' 
#' @note 
#' The S-RIP colour definitions for reanalysis datasets (Updated on 11 July 2017, for the S-RIP Full Report 2018)
#' https://s-rip.ees.hokudai.ac.jp/mediawiki/index.php/Notes_for_Authors
#' Reanalyses	R, G, B	Hexadecimal	Notes
#' MERRA-2	226, 31, 38	#E21F26	
#' MERRA	246, 153, 153	#F69999	
#' ERA-Interim	41, 95, 138	#295F8A	
#' ERA5	95, 152, 198	#5F98C6	
#' ERA-40	175, 203, 227	#AFCBE3	
#' JRA-55	114, 59, 122	#723B7A	
#' JRA-55C, JRA-55AMIP	173, 113, 181	#AD71B5	
#' JRA-25	214, 184, 218	#D6B8DA	
#' NCEP-NCAR R1	245, 126, 32	#F57E20	
#' NCEP-DOE R2	253, 191, 110	#FDBF6E	
#' 20CR v2c	236, 0, 140	#EC008C	
#' 20CR v2	247, 153, 209	#F799D1	
#' CERA-20C	0, 174, 239	#00AEEF	
#' ERA-20C	96, 200, 232	#60C8E8	
#' CFSR	52, 160, 72	#34A048	
#' REM	179, 91, 40	#B35B28	reanalysis ensemble mean
#' Other	255, 215, 0	#FFD700	
#' Observations	0, 0, 0	#000000	observations - black
#' Other observations	119, 119, 119	#777777	observations - grey
#' 
#' @export
#' 
reanalysis.colors <- data.frame(
  id = c(   "NR-1",    "NR-2",    "ERA-INT", "CFSR",    "JRA-55",  "JRA-55C", "20CR-2c", "ERA-20C", "MERRA-2", "CERA-20C"),
  color = c("#F57E20", "#FDBF6E", "#295F8A", "#34A048", "#723B7A", "#AD71B5", "#EC008C", "#60C8E8", "#E21F26", "#00AEEF")
)


#' Get the reference color for a giver reanalysis.
#'
#' Get the reference color for a giver reanalysis.
#'
#'
#' @return The hex color code
#'
#' @examples
#' \dontrun{
#' col <- atmoswing::getReanalysisColor('CFSR')
#' }
#' 
#' @export
#' 
getReanalysisColor <- function(reanalysis) {
  
  if(reanalysis == "NR-1") {
    return("#F57E20")
  } else if (reanalysis == "NR-2") {
    return("#FDBF6E")
  } else if (reanalysis == "ERA-INT") {
    return("#295F8A")
  } else if (reanalysis == "CFSR") {
    return("#34A048")
  } else if (reanalysis == "JRA-55") {
    return("#723B7A")
  } else if (reanalysis == "JRA-55C") {
    return("#AD71B5")
  } else if (reanalysis == "20CR-2c") {
    return("#EC008C")
  } else if (reanalysis == "ERA-20C") {
    return("#60C8E8")
  } else if (reanalysis == "MERRA-2") {
    return("#E21F26")
  } else if (reanalysis == "CERA-20C") {
    return("#00AEEF")
  } else {
    message("Reanalysis not found")
  }
}


#' Get the column max value.
#'
#' Get the column max value.
#'
#'
#' @return The maximum value found in the column
#'
#' @examples
#' \dontrun{
#' max <- atmoswing::colMax(data)
#' }
#' 
#' @export
#' 
colMax <- function(data) sapply(data, max, na.rm = TRUE)


#' Get the column min value.
#'
#' Get the column min value.
#'
#'
#' @return The minimum value found in the column
#'
#' @examples
#' \dontrun{
#' min <- atmoswing::colMin(data)
#' }
#' 
#' @export
#' 
colMin <- function(data) sapply(data, min, na.rm = TRUE)