# Check simulated data
# Chayan, 02/2017
#
# roxygen comments
#' Check simulated data
#'
#' \pkg{nca.check.sim} Checks simulated data for compatibility with ncappc and 
#' processes the data with various filtering criteria.
#' 
#' \pkg{nca.check.sim} Checks simulated data for compatibility with ncappc.
#' 
#' @param obsData Simulated concentration-time data.
#' @param idNmSim Column name for ID in simulated data. Default is \strong{"ID"}
#' @param timeNmSim Column name for time in simulated data. Default is 
#'   \strong{"TIME"}
#' @param concNmSim Column name for concentration in simulated data. Default is 
#'   \strong{"DV"}
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
#' @param adminType Route of administration. Allowed options are iv-bolus,
#'   iv-infusion or extravascular. Default is \strong{"extravascular"}
#' @param TI Infusion duration. If TI is a single numeric value, TI is the same
#'   for all individuals. If TI is the name of a column with numeric data
#'   present in the data set, TI is set to the unique value of the column for a
#'   given individual. Default is \strong{\code{NULL}}
#' @param doseAmtNm Column name to specify dose amount. Default is
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

nca.check.sim <- function(simData,
                          idNmSim="ID",timeNmSim="TIME",concNmSim="DV",
                          filterNm=NULL,filterExcl=NULL,
                          str1Nm=NULL,str1=NULL,
                          str2Nm=NULL,str2=NULL,
                          str3Nm=NULL,str3=NULL,
                          adminType="extravascular",TI=NULL,doseAmtNm=NULL,
                          blqNm=NULL,blqExcl=1,
                          evid=TRUE,evidIncl=0,mdv=FALSE){
  
  # Remove duplicated rows and columns
  simData <- simData[!duplicated(simData),]         # Remove duplicated rows
  simData <- simData[,!duplicated(names(simData))]  # Remove duplicated columns
  
  if ((!idNmSim%in%colnames(simData)) | (!timeNmSim%in%colnames(simData)) | (!concNmSim%in%colnames(simData))){
    stop("Incorrect column names of ID, TIME and/or DV in simulation output\n")
  }
  
  # Values in filterExcl will be excluded from the analysis
  if (!is.null(filterNm)){
    if(filterNm%in%colnames(simData) & !is.null(filterExcl)){
      for (i in 1:length(filterExcl)){
        if (grepl("^[-]?[0-9]*[.]?[0-9]*[eE]?[-]?[0-9]*[.]?[0-9]*$", filterExcl[i])){
          # Check if the filterExcl is numeric
          simData <- simData[simData[,filterNm] != filterExcl[i],]
        }else{
          # Check if the filterExcl is logical
          simData <- eval(parse(text=paste0("subset(simData, !", filterNm, filterExcl[i], ")")))
        }
      }
    }else{
      print("Note: Incorrect filterNm or filterExcl specification. filterNm will not be used to process the simulated data.\n")
    }
  }
  
  srdf <- simData[simData$NSUB == 1,]  # copy 1st replicate of simulated data before processing
  
  
  # 1st level population stratification
  if (!is.null(str1Nm)){
    if (!str1Nm%in%colnames(srdf)){stop("Incorrect name for the 1st level stratification column in simulation output\n")}
    if (is.null(str1)) str1 <- unique(sort(srdf[,str1Nm]))
  }
  
  # 2nd level population stratification
  if (!is.null(str2Nm)){
    if (!str2Nm%in%colnames(srdf)){stop("Incorrect name for the 2nd level stratification column in simulation output\n")}
    if (is.null(str2)) str2 <- unique(sort(srdf[,str2Nm]))
  }
  
  # 3rd level population stratification
  if (!is.null(str3Nm)){
    if (!str3Nm%in%colnames(srdf)){stop("Incorrect name for the 3rd level stratification column in simulation output\n")}
    if (is.null(str3)) str3 <- unique(sort(srdf[,str3Nm]))
  }
  
  
  # check requirements for infusion data
  if(adminType == "iv-infusion"){
    if(is.null(TI)){
      if((!"AMT"%in%names(srdf)) | (!"RATE"%in%names(srdf))){
        stop("Duration of the infusion time is needed if AMT and RATE are absent in the simulated data\n")
      }
    }else{
      if(!TI%in%names(srdf)){
        stop("Incorrect TI information for simulated data. TI has to be either a single numeric value or a column name.\n")
      }
    }
  }
  
  
  # Dose amount is extracted from doseAmtNm column
  if(!is.null(doseAmtNm)){
    if(doseAmtNm%in%colnames(srdf)){
      doseAmtNm <- doseAmtNm
    }else if("AMT"%in%colnames(srdf)){
      doseAmtNm <- "AMT"
    }else{
      doseAmtNm <- NULL
      print("Note: Dose amount column name provided in doseAmtNm or AMT column does not exist in the simulated data file. Dose related NCA metrics will not be estimated for the simulated data.\n")
    }
  }else{
    if("AMT"%in%colnames(srdf)) doseAmtNm <- "AMT"
  }
  
  
  # ignore data with BLQ = 1 or user specified value (optional)
  if (!is.null(blqNm)){
    if(blqNm%in%colnames(srdf) & !is.null(blqExcl)){
      for (i in 1:length(blqExcl)){
        if(grepl("^[-]?[0-9]*[.]?[0-9]*[eE]?[-]?[0-9]*[.]?[0-9]*$", blqExcl[i])){
          # Check if the blqExcl is numeric
          simData <- simData[simData[,blqNm] != blqExcl[i],]
        }else{
          # Check if the blqExcl is logical
          simData <- eval(parse(text=paste0("subset(simData, !", blqNm, blqExcl[i], ")")))
        }
      }
    }else{
      print("Note: Incorrect BLQ column name. BLQ will not be used to process the simulated data.\n")
    }
  }
  
  
  # include data based on specific values on EVID column (optional) but keep rows with TIME == 0
  if (evid){
    if("EVID"%in%colnames(srdf)){
      # uevid == unique values in EVID column
      # evidIncl == EVID values to be included
      # ievid == EVID values to be ignored
      uevid <- unique(as.numeric(as.character(srdf$EVID))); ievid <- setdiff(uevid, as.numeric(evidIncl))
      if(length(ievid) != 0) simData <- simData[!simData$EVID%in%ievid,]
    }else{
      print("Note: EVID column is not present. EVID will not be used to process the simulated data.\n")
    }
  }
  
  # if MDV fiter is present, exclude data for MDV == 1 but keep rows with TIME == 0
  if(mdv){
    if("MDV"%in%colnames(simData)){
      simData <- simData[simData$MDV == 0,]
    }else{
      print("Note: MDV column is not present. MDV will not be used to process the simulated data.\n")
    }
  }
  
  return(list(simData=simData,simRefData=srdf,str1=str1,str2=str2,str3=str3,doseAmtNm=doseAmtNm))
  
}
