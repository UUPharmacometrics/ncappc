
#' Prepare individual PK data
#'
#' Extracts time-conc data for a given individual.
#' 
#' Extracts time-conc data for a given individual.
#' 
#' @param pkData PK concentration-time data.
#' @param dvLog If \code{TRUE} concentration is in logarithmic scale. Default is
#'   \strong{\code{FALSE}}
#' @param dataType Indicates if the data is observed ("obs") or simulated 
#'   ("sim"). Since the simulated data is assumed to be obtained from NONMEM 
#'   output, DATE and clock time (\code{dateColNm}, \code{dateFormat}, 
#'   \code{timeFormat}) are not used for time data. Default is \strong{"obs"}.
#' @param idNm Column name for ID in PK data. Default is \strong{"ID"}
#' @param timeNm Column name for time in PK data. Default is 
#'   \strong{"TIME"}
#' @param concNm Column name for concentration in PK data. Default is 
#'   \strong{"DV"}
#' @param adminType Route of administration. Allowed options are iv-bolus,
#'   iv-infusion or extravascular. Default is \strong{"extravascular"}
#' @param TI Infusion duration. If TI is a single numeric value, TI is the same
#'   for all individuals. If TI is the name of a column with numeric data
#'   present in the data set, TI is set to the unique value of the column for a
#'   given individual. Default is \strong{\code{NULL}}
#' @param dateColNm column name for date if used (e.g. "Date", "DATE"). Default
#'   is \strong{\code{NULL}}
#' @param dateFormat date format (D-M-Y, D/M/Y or any other combination of
#'   D,M,Y). Default is \strong{\code{NULL}}
#' @param timeFormat time format (number, H:M, H:M:S). Default is
#'   \strong{"number"}
#' 
#' @return A list of objects with time-conc data and individual infusion
#'   duration for iv-infusion data
# @export
#'

nca_ind_data <- function(pkData, dvLog = FALSE,
                         dataType = "obs",
                         idNm="ID", timeNm="TIME", concNm="DV",
                         adminType="extravascular", TI=NULL,
                         dateColNm=NULL, dateFormat=NULL, timeFormat="number"){
  
  pkData <- data.frame(pkData)
  ID <- pkData[1,idNm]
  if(adminType == "iv-infusion"){
    if(is.null(TI)){
      amt  <- pkData[pkData$AMT > 0,"AMT"][1]
      rate <- pkData[pkData$RATE > 0,"RATE"][1]
      if (is.na(amt) | is.na(rate) | rate==0){
        stop(paste0("Incorrect AMT and/or RATE value for ", ID))
      }else{
        iTI <- amt/rate
      }
    }else{
      if(grepl("^[-]?[0-9]*[.]?[0-9]*[eE]?[-]?[0-9]*[.]?[0-9]*$", TI)){
        iTI <- TI
      }else{
        TItmp <- unique(pkData[, TI])
        if(length(TItmp)!=1 || !grepl("^[-]?[0-9]*[.]?[0-9]*[eE]?[-]?[0-9]*[.]?[0-9]*$", TItmp)){
          stop(paste0("Incorrect TI information. TI column has either more than one unique value or non-numeric data for the individual with ID = ", ID, ".\n"))
        }else{
          iTI <- as.numeric(as.character(TItmp))
        }
      }
    }
  }else{
    iTI <- NaN
  }
  
  if(dataType=="obs" && timeFormat != "number"){
    time <- numeric(0)
    if (!is.null(dateColNm)){
      tm <- as.POSIXct(paste(pkData[,dateColNm],pkData[, timeNm]), format=paste(dateFormat,timeFormat,sep=" "))
    }else{
      tm <- pkData[, timeNm]
    }
    for(j in 1:length(tm)){
      time[j] <- ifelse ((is.null(dateColNm)), as.numeric(difftime(strptime(tm[j], format=timeFormat), strptime(tm[1], format=timeFormat), units='hours')), as.numeric(difftime(strptime(tm[j], format="%Y-%m-%d %H:%M:%S"), strptime(tm[1], format="%Y-%m-%d %H:%M:%S"), units='hours')))
    }
  }else{
    time <- suppressWarnings(as.numeric(as.character(pkData[, timeNm])))
  }
  
  conc <- as.character(pkData[, concNm])
  conc <- unname(sapply(conc, function(x) ifelse(x==".", 0, as.numeric(x))))
  
  if(length(which(is.na(time)))!=0){
    zidx <- which(is.na(time))
    time <- time[-zidx]
    conc <- conc[-zidx]
  }
  if(length(which(is.na(conc)))!=0){
    zidx <- which(is.na(conc))
    time <- time[-zidx]
    conc <- conc[-zidx]
  }
  
  if(dvLog) conc <- sapply(conc, function(x) ifelse(x==0, 0, exp(x)))
  
  tc <- data.frame(time,conc)
  if(nrow(tc)>0) tc <- tc[order(tc$time),]
  
  tc <- tibble::as_tibble(tc)
  tc$iTI <- iTI
  return(tc)
}

