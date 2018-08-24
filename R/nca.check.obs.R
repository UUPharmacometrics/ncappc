# Check observed data
# Chayan, 02/2017
#
# roxygen comments
#' Check observed data
#'
#' \pkg{nca.check.obs} Checks observed data for compatibility with ncappc and
#' processes the data with various filtering criteria.
#' 
#' \pkg{nca.check.obs} Checks observed data for compatibility with ncappc.
#' 
#' @param obsData Observed concentration-time data.
#' @param idNmObs Column name for ID in observed data. Default is \strong{"ID"}
#' @param timeNmObs Column name for time in observed data. Default is 
#'   \strong{"TIME"}
#' @param concNmObs Column name for concentration in observed data. Default is 
#'   \strong{"DV"}
#' @param doseType Steady-state (ss) or non-steady-state (ns) dose. Default is 
#'   \strong{"ns"}
#' @param doseTime Dose time prior to the first observation for steady-state
#'   data. Default is \strong{\code{NULL}}
#' @param Tau Dosing interval for steady-state data. Default is
#'   \strong{\code{NULL}}
#' @param filterNm Column name to filter data. Default is \strong{\code{NULL}}
#' @param filterExcl Row exclusion criteria based on the column defined by
#'   \code{filterNm}. This can be numeric value or logical condition (e.g. c(1,
#'   2, "<20", ">=100", "!=100")). Default is \strong{\code{NULL}}
#' @param str1Nm Column name for 1st level population stratifier. Default is
#'   \strong{\code{NULL}}
#' @param str1 Stratification ID of the members within 1st level stratification
#'   (e.g c(1,2)). Default is \strong{\code{NULL}}
#' @param str2Nm Column name for 2nd level population stratifier. Default is
#'   \strong{\code{NULL}}
#' @param str2 Stratification ID of the members within 2nd level stratification
#'   (e.g c(1,2)). Default is \strong{\code{NULL}}
#' @param str3Nm Column name for 3rd level population stratifier. Default is
#'   \strong{\code{NULL}}
#' @param str3 Stratification ID of the members within 3rd level stratification
#'   (e.g c(1,2)). Default is \strong{\code{NULL}}
#' @param AUCTimeRange User-defined window of time used to estimate AUC. Default
#'   is \strong{\code{NULL}}
#' @param LambdaTimeRange User-defined window of time to estimate elimination
#'   rate-constant. This argument lets the user to choose a specific window of
#'   time to be used to estimate the elimination rate constant (Lambda) in the
#'   elimination phase. The accepted format for the input to this argument is a
#'   numeric array of two elements; \code{c(14,24)} will estimate the Lambda
#'   using the data within the time units 14 to 24. Default is
#'   \strong{\code{NULL}}
#' @param adminType Route of administration. Allowed options are iv-bolus,
#'   iv-infusion or extravascular. Default is \strong{"extravascular"}
#' @param TI Infusion duration. If TI is a single numeric value, TI is the same
#'   for all individuals. If TI is the name of a column with numeric data
#'   present in the data set, TI is set to the unique value of the column for a
#'   given individual. Default is \strong{\code{NULL}}
#' @param doseAmtNm Column name to specify dose amount. Default is
#'   \strong{\code{NULL}}
#' @param dateColNm column name for date if used (e.g. "Date", "DATE"). Default
#'   is \strong{\code{NULL}}
#' @param dateFormat date format (D-M-Y, D/M/Y or any other combination of
#'   D,M,Y). Default is \strong{\code{NULL}}
#' @param timeFormat time format (number, H:M, H:M:S). Default is
#'   \strong{"number"}
#' @param concUnit Unit of concentration (e.g. "ng/mL"). Default is 
#'   \strong{\code{NULL}}
#' @param timeUnit Unit of time (e.g. "h"). Default is \strong{\code{NULL}}
#' @param doseUnit Unit of dose amount (e.g. "ng"). Default is
#'   \strong{\code{NULL}}
#' @param blqNm Name of BLQ column if used to exclude data. Default is
#'   \strong{\code{NULL}}
#' @param blqExcl Excluded BLQ value; either a numeric value or a logical
#'   condition (e.g. 1 or ">=1" or c(1,">3")). Used only if the \code{blqNm} is
#'   not \code{NULL}. Default is \strong{"1"}
#' @param evid If \code{TRUE} EVID is used to filter data. Default is
#'   \strong{\code{TRUE}}
#' @param evidIncl Included values in EVID. Default is \strong{"0"}
#' @param mdv If \code{TRUE} MDV is used to include data when MDV=0. Default is
#'   \strong{\code{FALSE}}
#' 
#' @return A list of objects
#' @export
#'


nca.check.obs <- function(obsData,
                          idNmObs="ID",timeNmObs="TIME",concNmObs="DV",
                          doseType="ns",doseTime=NULL,Tau=NULL,
                          filterNm=NULL,filterExcl=NULL,
                          str1Nm=NULL,str1=NULL,
                          str2Nm=NULL,str2=NULL,
                          str3Nm=NULL,str3=NULL,
                          AUCTimeRange=NULL,LambdaTimeRange=NULL,
                          adminType="extravascular",TI=NULL,doseAmtNm=NULL,
                          dateColNm=NULL,dateFormat=NULL,timeFormat="number",
                          concUnit=NULL,timeUnit=NULL,doseUnit=NULL,
                          blqNm=NULL,blqExcl=1,
                          evid=TRUE,evidIncl=0,mdv=FALSE){
  
  # Remove duplicated rows and columns
  obsData <- obsData[!duplicated(obsData),]         # Remove duplicated rows
  obsData <- obsData[,!duplicated(names(obsData))]  # Remove duplicated columns
  
  # Check for column names
  if ((!idNmObs%in%names(obsData)) | (!timeNmObs%in%names(obsData)) | (!concNmObs%in%names(obsData))){
    stop("Incorrect column names of ID, TIME and/or DV.")
  }
  
  # Check steady state dosing interval
  if (doseType == "ss"){
    if(is.null(Tau)){stop("Time for dosing interval is required for steady-state data.")}
    if(is.null(doseTime)){print("Note: Dose time prior to the first observation for steady-state data is missing. Steady state observation period will be estimated based on the first non-zero obverved concentration and Tau.")}
  }
  
  # exclude data based on specific values on filter column (optional)
  # Values in filterExcl will be excluded from the analysis
  if (!is.null(filterNm)){
    if(filterNm%in%names(obsData) & !is.null(filterExcl)){
      for (i in 1:length(filterExcl)){
        # Check if the filterExcl is numeric
        if (grepl("^[-]?[0-9]*[.]?[0-9]*[eE]?[-]?[0-9]*[.]?[0-9]*$", filterExcl[i])){
          obsData <- obsData[obsData[,filterNm] != filterExcl[i],]
        }else{
          # Check if the filterExcl is logical
          obsData <- eval(parse(text=paste0("subset(obsData, !", filterNm, filterExcl[i], ")")))
        }
      }
    }else{
      print("Note: Incorrect filterNm or filterExcl specification. filterNm will not be used to process the observed data.")
    }
  }
  
  # copy the input data before processing
  refdf <- obsData
  
  # 1st level population stratification
  if (!is.null(str1Nm)){
    if (!str1Nm%in%names(obsData)){stop("Incorrect name for the 1st level stratification column.")}
    if (is.null(str1)) str1 <- unique(sort(as.character(obsData[,str1Nm])))
  }
  
  # 2nd level population stratification
  if (!is.null(str2Nm)){
    if (!str2Nm%in%names(obsData)){stop("Incorrect name for the 2nd level stratification column.")}
    if (is.null(str2)) str2 <- unique(sort(as.character(obsData[,str2Nm])))
  }
  
  # 3rd level population stratification
  if (!is.null(str3Nm)){
    if (!str3Nm%in%names(obsData)){stop("Incorrect name for the 3rd level stratification column.")}
    if (is.null(str3)) str3 <- unique(sort(as.character(obsData[,str3Nm])))
  }
  
  # check time range, if any
  if ((!is.null(AUCTimeRange)) && (length(AUCTimeRange) != 2 | class(AUCTimeRange) != "numeric")){
    print("Note: Incorrect time range for AUC calculation. AUCTimeRange will not be used.")
  }
  
  if ((!is.null(LambdaTimeRange)) && (length(LambdaTimeRange) != 2 | class(LambdaTimeRange) != "numeric")){
    print("Note: Incorrect time range for Lambda calculation. LambdaTimeRange will not be used.")
    LambdaTimeRange <- NULL
  }
  
  # check requirements for infusion data
  if(adminType == "iv-infusion"){
    if(is.null(TI)){
      if((!"AMT"%in%names(obsData)) | (!"RATE"%in%names(obsData))){
         stop("Duration of the infusion time is needed if AMT and RATE are absent in the input data.")
      }
      TInum <- NULL
    }else{
      if(length(TI)!=1){
         stop("Incorrect TI information. TI has to be either a single numeric value or a column name.")
      }else{
        if(grepl("^[-]?[0-9]*[.]?[0-9]*[eE]?[-]?[0-9]*[.]?[0-9]*$", TI)){
          TI    <- as.numeric(TI)
          TInum <- TRUE
        }else{
          if(!TI%in%names(obsData)){
             stop("Incorrect TI information. TI has to be either a single numeric value or a column name.")
          }else{
            TInum <- FALSE
          }
        }
      }
    }
  }
  
  # Dose amount is extracted from doseAmtNm column
  if(!is.null(doseAmtNm)){
    if(doseAmtNm%in%names(obsData)){
      doseAmtNm <- doseAmtNm
    }else if("AMT"%in%names(obsData)){
      doseAmtNm <- "AMT"
    }else{
      doseAmtNm <- NULL
      print("Note: Dose amount column name provided in doseAmtNm or AMT column does not exist in the observed data file. Dose related NCA metrics will not be estimated for the observed data.")
    }
  }else{
    if("AMT"%in%names(obsData)) doseAmtNm <- "AMT"
  }
  
  
  # Appropriate date format
  if (!is.null(dateColNm)){
    if (!dateColNm%in%names(obsData)){stop("Incorrect name for the date column.")}
    if (dateFormat == "D-M-Y"){
      dateFormat <- "%d-%m-%Y"
    }else if (dateFormat == "D-Y-M"){
      dateFormat <- "%d-%Y-%m"
    }else if (dateFormat == "M-D-Y"){
      dateFormat <- "%m-%d-%Y"
    }else if (dateFormat == "M-Y-D"){
      dateFormat <- "%m-%Y-%d"
    }else if (dateFormat == "Y-M-D"){
      dateFormat <- "%Y-%m-%d"
    }else if (dateFormat == "Y-D-M"){
      dateFormat <- "%Y-%d-%m"
    }else if (dateFormat == "D/M/Y"){
      dateFormat <- "%d/%m/%Y"
    }else if (dateFormat == "D/Y/M"){
      dateFormat <- "%d/%Y/%m"
    }else if (dateFormat == "M/D/Y"){
      dateFormat <- "%m/%d/%Y"
    }else if (dateFormat == "M/Y/D"){
      dateFormat <- "%m/%Y/%d"
    }else if (dateFormat == "Y/M/D"){
      dateFormat <- "%Y/%m/%d"
    }else if (dateFormat == "Y/D/M"){
      dateFormat <- "%Y/%d/%m"
    }else{stop("Incorrect date format. Currently allowed date formats are \"D-M-Y\", \"D/M/Y\", or any combination of D, M and Y separated by either - or /.")}
  }
  
  # Appropriate time format
  if (timeFormat != "number"){
    if(timeFormat == "H:M"){
      timeFormat <- "%H:%M"
    }else if(timeFormat == "H:M:S"){
      timeFormat <- "%H:%M:%S"
    }else{stop("Incorrect time format. Currently allowed date formats are \"number\", \"H:M\", \"H:M:S\".")}
  }
  
  
  # Units for dose, time and conc
  if(is.null(doseUnit)){dunit <- NULL}else{dunit <- doseUnit}
  if(is.null(timeUnit)){tunit <- NULL}else{tunit <- timeUnit}
  if(is.null(concUnit)){cunit <- NULL}else{cunit <- concUnit}
  if(is.null(tunit) | is.null(cunit))   {aucunit  <- NULL}else{aucunit  <- paste0(tunit,"*",cunit)}
  if(is.null(tunit) | is.null(cunit))   {aumcunit <- NULL}else{aumcunit <- paste0("(",tunit,"^2)*",cunit)}
  if(is.null(dunit) | is.null(aucunit)) {clunit   <- NULL}else{clunit   <- paste0(dunit,"/(",aucunit,")")}
  if(is.null(dunit) | is.null(cunit))   {vlunit   <- NULL}else{vlunit   <- paste0(dunit,"/",cunit)}
  
  
  # ignore data with BLQ = 1 or user specified value (optional)
  if(!is.null(blqNm)){
    if(blqNm%in%names(obsData) & !is.null(blqExcl)){
      for (i in 1:length(blqExcl)){
        if (grepl("^[-]?[0-9]*[.]?[0-9]*[eE]?[-]?[0-9]*[.]?[0-9]*$", blqExcl[i])){
          # Check if the blqExcl is numeric
          obsData <- obsData[obsData[,blqNm] != blqExcl[i],]
        }else{
          # Check if the blqExcl is logical
          obsData <- eval(parse(text=paste0("subset(obsData, !", blqNm, blqExcl[i], ")")))
        }
      }
    }else{
      print("Note: Incorrect BLQ column name. BLQ will not be used to process the observed data.")
    }
  }
  
  
  # include data based on specific values on EVID column (optional) but keep rows with TIME == 0
  if (evid){
    if("EVID"%in%names(obsData)){
      # uevid == unique values in EVID column
      # evidIncl == EVID values to be included
      # ievid == EVID values to be ignored
      uevid <- unique(as.numeric(as.character(obsData$EVID))); ievid <- setdiff(uevid, as.numeric(evidIncl))
      if (length(ievid) != 0) obsData <- obsData[!obsData$EVID%in%ievid,]
    }else{
      print("Note: EVID column is not present. EVID will not be used to process the observed data.")
    }
  }
  
  # if MDV fiter is present, exclude data for MDV == 1 but keep rows with TIME == 0
  if(mdv){
    if("MDV"%in%names(obsData)){
      obsData <- obsData[obsData$MDV == 0,]
    }else{
      print("Note: MDV column is not present. MDV will not be used to process the observed data.")
    }
  }
  
  if(!exists("TInum")) TInum <- NA
  
  # Remove rows with NA or empty observation/time
  if(length(which(is.na(obsData[,concNmObs]) | obsData[,concNmObs]=="")) != 0){
    obsData <- obsData[-which(is.na(obsData[,concNmObs]) | obsData[,concNmObs]==""),]
  }
  if(length(which(is.na(obsData[,timeNmObs]) | obsData[,timeNmObs]=="")) != 0){
    obsData <- obsData[-which(is.na(obsData[,timeNmObs]) | obsData[,timeNmObs]==""),]
  }
  
  return(list(obsData=obsData,refdf=refdf,str1=str1,str2=str2,str3=str3,LambdaTimeRange=LambdaTimeRange,
              TI=TI,TInum=TInum,dateFormat=dateFormat,timeFormat=timeFormat,dunit=dunit,tunit=tunit,cunit=cunit,
              aucunit=aucunit,aumcunit=aumcunit,clunit=clunit,vlunit=vlunit,doseAmtNm=doseAmtNm))
  
}


