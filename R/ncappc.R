# This R-script performs diagnostic tests on PK/PD models using non-compartmental analysis parameters
# Chayan, 11/2014

#  roxygen comments
#' Performs NCA calculations and population PK model diagnosis.
#'
#' \pkg{ncappc} is a flexible tool, to
#' \enumerate{
#'  \item perform a traditional NCA
#'  \item perform simulation-based posterior predictive checks for a
#' population PK model using NCA metrics.
#' }
#'
#' Non-compartmental analysis (NCA) calculates pharmacokinetic (PK) metrics 
#' related to the systemic exposure to a drug following administration, e.g. 
#' area under the concentration-time curve and peak concentration. \pkg{ncappc} 
#' performs a traditional NCA using the observed plasma concentration-time data.
#' In the presence of simulated plasma concentration-time data, \pkg{ncappc} 
#' also performs simulation-based posterior predictive checks (ppc) using NCA 
#' metrics for the corresponding population PK (PopPK) model used to generate 
#' the simulated data. The diagnostic analysis is performed at the population as
#' well as the individual level. The distribution of the simulated population 
#' means of each NCA metric is compared with the corresponding observed 
#' population mean. The individual level comparison is performed based on the 
#' deviation of the mean of any NCA metric based on simulations for an 
#' individual from the corresponding NCA metric obtained from the observed data.
#' Additionaly, \pkg{ncappc} reports the normalized prediction distribution 
#' error (NPDE) of the simulated NCA metrics for each individual and their 
#' distribution within a population. \pkg{ncappc} produces two default outputs 
#' depending on the type of analysis performed, i.e., traditional NCA and PopPK 
#' diagnosis. The PopPK diagnosis feature of \pkg{ncappc} produces 7 sets of 
#' graphical outputs to assess the ability of a population model to simulate the
#' concentration-time profile of a drug and thereby identify model 
#' misspecification. In addition, tabular outputs are generated showing the 
#' values of the NCA metrics estimated from the observed and the simulated data,
#' along with the deviation, NPDE, regression parameters used to estimate the 
#' elimination rate constant and the related population statistics. The default 
#' values of the arguments used in \pkg{ncappc} are shown in the \strong{Useage}
#' section of this document and/or in \strong{bold} in the \strong{Arguments} 
#' section.
#'
#' @param obsFile Observed concentration-time data from an internal data frame
#'   or an external table with comma, tab or space as separators. 
#' @param simFile NONMEM simulation output with the simulated concentration-time
#'   data from an internal data frame or an external table. \code{NULL} produces
#'   just the NCA output, a filename or data frame prduces the NCA output as
#'   well as the PopPK diagnosis. If \code{new_data_method=TRUE} then this can
#'   be a compressed file as well. 
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
#' @param concUnit Unit of concentration (e.g. "ng/mL"). Default is 
#'   \strong{\code{NULL}}
#' @param timeUnit Unit of time (e.g. "h"). Default is \strong{\code{NULL}}
#' @param doseUnit Unit of dose amount (e.g. "ng"). Default is
#'   \strong{\code{NULL}}
#' @param obsLog If \code{TRUE} concentration in observed data is in logarithmic
#'   scale. Default is \strong{\code{FALSE}}
#' @param simLog If \code{TRUE} concentration in simulated data is in 
#'   logarithmic scale. Default is \strong{\code{FALSE}}
#' @param psnOut If \code{TRUE} observed data is an output from PsN or in NONMEM
#'   output format. Default is \strong{\code{TRUE}}
#' @param idNmObs Column name for ID in observed data. Default is \strong{"ID"}
#' @param timeNmObs Column name for time in observed data. Default is
#'   \strong{"TIME"}
#' @param concNmObs Column name for concentration in observed data. Default is
#'   \strong{"DV"}
#' @param idNmSim Column name for ID in simulated data. Default is \strong{"ID"}
#' @param timeNmSim Column name for time in simulated data. Default is
#'   \strong{"TIME"}
#' @param concNmSim Column name for concentration in simulated data. Default is
#'   \strong{"DV"}
#' @param onlyNCA If \code{TRUE} only NCA is performed and ppc part is ignored
#'   although simFile is not \code{NULL}. Default is \strong{\code{FALSE}}
#' @param AUCTimeRange User-defined window of time used to estimate AUC. Default
#'   is \strong{\code{NULL}}
#' @param backExtrp If \code{TRUE} back-extrapolation is performed while
#'   estimating AUC. Default is \strong{\code{FALSE}}
#' @param LambdaTimeRange User-defined window of time to estimate elimination
#'   rate-constant. This argument lets the user to choose a specific window of
#'   time to be used to estimate the elimination rate constant (Lambda) in the
#'   elimination phase. The accepted format for the input to this argument is a
#'   numeric array of two elements; \code{c(14,24)} will estimate the Lambda
#'   using the data within the time units 14 to 24. Default is
#'   \strong{\code{NULL}}
#' @param LambdaExclude User-defined excluded observation time points for
#'   estimation of Lambda. This can be numeric value or logical condition (e.g.
#'   c(1, 2, "<20", ">=100", "!=100")). Default is \strong{\code{NULL}}
#' @param doseAmtNm Column name to specify dose amount. Default is
#'   \strong{\code{NULL}}
#' @param adminType Route of administration. Allowed options are iv-bolus,
#'   iv-infusion or extravascular. Default is \strong{"extravascular"}
#' @param doseType Steady-state (ss) or nonsteady-state (ns) dose. Default is
#'   \strong{"ns"}
#' @param doseTime Dose time prior to the first observation for steady-state
#'   data. Default is \strong{\code{NULL}}
#' @param Tau Dosing interval for steady-state data. Default is
#'   \strong{\code{NULL}}
#' @param TI Infusion duration. If TI is a single numeric value, TI is the same
#'   for all individuals. If TI is the name of a column with numeric data
#'   present in the data set, TI is set to the unique value of the column for a
#'   given individual. Default is \strong{\code{NULL}}
#' @param method Method to estimate AUC. \code{linear} method applies the linear
#'   trapezoidal rule to estimate the area under the curve. \code{log} method 
#'   applies the logarithmic trapezoidal rule to estimate the area under the 
#'   curve. \code{linearup-logdown} method applies the linear trapezoidal rule
#'   to estimate the area under the curve for the ascending part of the curve
#'   and the logarithmic trapezoidal rule to estimate the area under the curve
#'   for the descending part of the curve. Default is
#'   \strong{"linearup-logdown"}
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
#' @param filterNm Column name to filter data. Default is \strong{\code{NULL}}
#' @param filterExcl Row exclusion criteria based on the column defined by
#'   \code{filterNm}. This can be numeric value or logical condition (e.g. c(1,
#'   2, "<20", ">=100", "!=100")). Default is \strong{\code{NULL}}
#' @param negConcExcl If \code{TRUE} negative concentrations are excluded.
#'   Default is \strong{\code{FALSE}}
#' @param param NCA parameters (AUClast, AUClower_upper, AUCINF_obs,
#'   AUCINF_pred, AUMClast, Cmax, Tmax, HL_Lambda_z). Default is
#'   \strong{(c"AUClast", "Cmax")}
#' @param timeFormat time format (number, H:M, H:M:S). Default is
#'   \strong{"number"}
#' @param dateColNm colunm name for date if used (e.g. "Date", "DATE"). Default
#'   is \strong{\code{NULL}}
#' @param dateFormat date format (D-M-Y, D/M/Y or any other combination of
#'   D,M,Y). Default is \strong{\code{NULL}}
#' @param spread Measure of the spread of simulated data (\code{"ppi"} (95\%
#'   parametric prediction interval) or \code{"npi"} (95\% nonparametric
#'   prediction interval)). Default is \strong{"npi"}
#' @param tabCol Output columns to be printed in the report in addition to ID,
#'   dose and population strata information (list of NCA metrics in a string
#'   array). Default is \strong{c("AUClast", "Cmax", "Tmax", "AUCINF_obs",
#'   "Vz_obs", "Cl_obs", "HL_Lambda_z")}
#' @param figFormat format of the produced figures (bmp, jpeg, tiff, png).
#'   Default is \strong{"tiff"}
#' @param noPlot If \code{TRUE} only NCA calculations are performed without any
#'   plot generation. Default is \strong{\code{FALSE}}
#' @param printOut If \code{TRUE} tabular and graphical outputs are saved on the
#'   disk. Default is \strong{\code{TRUE}}
#' @param studyName Name of the study to be added as a description in the
#'   report. Default is \strong{\code{NULL}}
#' @param new_data_method If \code{TRUE} a faster method of reading data is
#'   tested. Default is \strong{\code{TRUE}}
#' @param overwrite_SIMDATA If \code{TRUE} new information is created in the
#'   SIMDATA directory. If \code{FALSE} the information in the SIMDATA directory
#'   is used. If \code{NULL} a dialog will come up to ask the user what to do.
#'   Default is \strong{\code{NULL}}
#' @param outFileNm Additional tag to the name of the output html and pdf output
#'   file hyphenated to the standard ncappc report file name standard ncappc
#'   report file name. Default is the \strong{\code{NULL}}
#' 
#' @import ggplot2
#' @import grid
#' @import gridExtra
#' @import scales
#' @import gtable
#' @import knitr
#' @import Cairo
#' @import xtable
#' @import reshape2
#' 
#' @return NCA results and diagnostic test results
#' @export
#' @examples
#' ncappc(obsFile=system.file("extdata","pkdata.csv",package="ncappc"), 
#' psnOut=FALSE,noPlot=TRUE,printOut=FALSE)
#'

ncappc <- function(obsFile="nca_original.npctab.dta",
                   simFile="nca_simulation.1.npctab.dta.zip",
                   str1Nm=NULL,str1=NULL,
                   str2Nm=NULL,str2=NULL,
                   str3Nm=NULL,str3=NULL,
                   concUnit=NULL,timeUnit=NULL,doseUnit=NULL,
                   obsLog=FALSE,simLog=obsLog,
                   psnOut=TRUE,
                   idNmObs="ID",timeNmObs="TIME",concNmObs="DV",
                   idNmSim="ID",timeNmSim="TIME",concNmSim="DV",
                   onlyNCA=FALSE,
                   AUCTimeRange=NULL,backExtrp=FALSE,
                   LambdaTimeRange=NULL,LambdaExclude=NULL,doseAmtNm=NULL,
                   adminType="extravascular",doseType="ns",doseTime=NULL,Tau=NULL,
                   TI=NULL,method="linearup-logdown",blqNm=NULL,blqExcl=1,
                   evid=TRUE,evidIncl=0,mdv=FALSE,filterNm=NULL,filterExcl=NULL,
                   negConcExcl=FALSE,param=c("AUClast","Cmax"),timeFormat="number",
                   dateColNm=NULL,dateFormat=NULL,spread="npi",
                   tabCol=c("AUClast","Cmax","Tmax","AUCINF_obs","Vz_obs","Cl_obs","HL_Lambda_z"),
                   figFormat="tiff",noPlot=FALSE,printOut=TRUE,studyName=NULL,new_data_method=TRUE,
                   overwrite_SIMDATA=NULL,outFileNm=NULL,
                   out_format = "html"){
  
  "..density.." <- "meanObs" <- "sprlow" <- "sprhgh" <- "AUClast" <- "AUCINF_obs" <- "Cmax" <- "Tmax" <- "FCT" <- "ID" <- "STR1" <- "STR2" <- "STR3" <- "NPDE" <- "mcil" <- "mciu" <- "sdu" <- "sducil" <- "sduciu" <- "scale_linetype_manual" <- "scale_color_manual" <- "xlab" <- "ylab" <- "guides" <- "guide_legend" <- "theme" <- "element_text" <- "unit" <- "element_rect" <- "geom_histogram" <- "aes" <- "geom_vline" <- "grid.arrange" <- "unit.c" <- "grid.grab" <- "ggsave" <- "facet_wrap" <- "ggplot" <- "labs" <- "geom_point" <- "geom_errorbarh" <- "knit2html" <- "knit2pdf" <- "knit" <- "file_test" <- "tail" <- "read.csv" <- "read.table" <- "dev.off" <- "write.table" <- "head" <- "write.csv" <- "coef" <- "dist" <- "lm" <- "median" <- "na.omit" <- "percent" <- "qchisq" <- "qnorm" <- "qt" <- "quantile" <- "scale_y_continuous" <- "sd" <- "STRAT1" <- "STRAT2" <- "STRAT3" <- "sdcil" <- "sdciu" <- "str" <- NULL
  rm(list=c("..density..","meanObs","sprlow","sprhgh","AUClast","AUCINF_obs","Cmax","Tmax","FCT","ID","STR1","STR2","STR3","NPDE","mcil","mciu","sdu","sducil","sduciu","scale_linetype_manual","scale_color_manual","xlab","ylab","guides","guide_legend","theme","element_text","unit","element_rect","geom_histogram","aes","geom_vline","grid.arrange","unit.c","grid.grab","ggsave","facet_wrap","ggplot","labs","geom_point","geom_errorbarh","knit2html","knit2pdf","knit","file_test","tail","read.csv","read.table","dev.off","write.table","head","write.csv","coef","dist","lm","median","na.omit","percent","qchisq","qnorm","qt","quantile","scale_y_continuous","sd","STRAT1","STRAT2","STRAT3","sdcil","sdciu","str"))
  
  options(warning.length=5000)
  options(scipen=999)
  theme_set(theme_bw(base_size=22))
  usrdir <- getwd()
  alwprm <- c("AUClast","AUClower_upper","AUCINF_obs","AUCINF_pred","AUMClast","Cmax","Tmax","HL_Lambda_z")
  
  
  # Observed data
  if (is.null(obsFile)){stop("Name of the file with observed data is required.")}
  if (!is.data.frame(obsFile)){
    if (!file_test("-f", obsFile)){stop("File for the observed data does not exist.")}
    # read observed data file
    if (!psnOut){
      extn <- tail(unlist(strsplit(obsFile, ".", fixed=T)), n=1)
      if(extn=="csv"){indf <- read.csv(obsFile)}else{indf <- read.table(obsFile, header=T)}
    }else{
      message("Note: The observed data file is expected to be generated by PsN since psnOut is set to TRUE.")
      #indf <- read.table(obsFile, header=T, skip=1)
      indf <- read_nm_table(obsFile)
      indf <- as.data.frame(indf)
    }
  }else{
    indf <- obsFile
  }
  
  
  # Simulated data is not supplied
  if (is.null(simFile)){
    message("Note: Simulated data file is not supplied. Only NCA module will be executed.")
  }else{
    if (!is.data.frame(simFile) && !file_test("-f", simFile)){
      message(paste0("Note: Simulated data file, ",simFile,", is not found in the working directory. Only NCA module will be executed."))
      simFile <- NULL
    }else{
      backExtrp <- FALSE
    }
  }
  
  
  # Check observed data
  obsList <- nca.check.obs(obsData=indf,
                           idNmObs=idNmObs, timeNmObs=timeNmObs, concNmObs=concNmObs,
                           doseType=doseType, doseTime=doseTime, Tau=Tau,
                           filterNm=filterNm, filterExcl=filterExcl,
                           str1Nm=str1Nm, str1=str1,
                           str2Nm=str2Nm, str2=str2,
                           str3Nm=str3Nm, str3=str3,
                           AUCTimeRange=AUCTimeRange, LambdaTimeRange=LambdaTimeRange,
                           adminType=adminType, TI=TI, doseAmtNm=doseAmtNm,
                           dateColNm=dateColNm, dateFormat=dateFormat, timeFormat=timeFormat,
                           concUnit=concUnit, timeUnit=timeUnit, doseUnit=doseUnit,
                           blqNm=blqNm, blqExcl=blqExcl,evid=evid, evidIncl=evidIncl, mdv=mdv)
  
  indf            <- obsList$obsData
  refdf           <- obsList$refdf
  str1            <- obsList$str1
  str2            <- obsList$str2
  str3            <- obsList$str3
  LambdaTimeRange <- obsList$LambdaTimeRange
  TI              <- obsList$TI
  TInum           <- obsList$TInum
  dateFormat      <- obsList$dateFormat
  timeFormat      <- obsList$timeFormat
  dunit           <- obsList$dunit
  tunit           <- obsList$tunit
  cunit           <- obsList$cunit
  aucunit         <- obsList$aucunit
  aumcunit        <- obsList$aumcunit
  clunit          <- obsList$clunit
  vlunit          <- obsList$vlunit
  doseAmtNm       <- obsList$doseAmtNm
  rm(obsList)
  
  npr       <- length(param)
  obsFileNm <- ifelse(is.data.frame(obsFile), deparse(substitute(obsFile)), obsFile)
  outData   <- data.frame()  # Create empty data frame for output
  
  
  # cpopStrNm = Combined stratifyting column names
  # npopStr   = Number of stratification levels
  # popStrNm1 = 1st level stratifying column name
  # popStr1   = 1st level stratification ID names
  # npopStr1  = Number of 1st level stratification ID names
  
  if (is.null(str1Nm)  & is.null(str2Nm)  & is.null(str3Nm)) {case<-1; npopStr<-0} # No stratification
  
  if (!is.null(str1Nm) & is.null(str2Nm)  & is.null(str3Nm)) {case<-2; cpopStrNm<-str1Nm; npopStr<-1; popStrNm1<-str1Nm; popStr1<-str1; npopStr1<-length(str1)} # Str1
  if (is.null(str1Nm)  & !is.null(str2Nm) & is.null(str3Nm)) {case<-2; cpopStrNm<-str2Nm; npopStr<-1; popStrNm1<-str2Nm; popStr1<-str2; npopStr1<-length(str2)} # Str2
  if (is.null(str1Nm)  & is.null(str2Nm)  & !is.null(str3Nm)){case<-2; cpopStrNm<-str3Nm; npopStr<-1; popStrNm1<-str3Nm; popStr1<-str3; npopStr1<-length(str3)} # Str3
  
  # Str1 & Str2
  if (!is.null(str1Nm) & !is.null(str2Nm) & is.null(str3Nm)){
    case<-3; cpopStrNm<-paste(str1Nm,str2Nm,sep=", "); npopStr<-2; popStrNm1<-str1Nm; popStrNm2<-str2Nm; popStr1<-str1; popStr2<-str2; npopStr1<-length(str1); npopStr2<-length(str2)
  }
  # Str1 & Str3
  if (!is.null(str1Nm) & is.null(str2Nm) & !is.null(str3Nm)){
    case<-3; cpopStrNm<-paste(str1Nm,str3Nm,sep=", "); npopStr<-2; popStrNm1<-str1Nm; popStrNm2<-str3Nm; popStr1<-str1; popStr2<-str3; npopStr1<-length(str1); npopStr2<-length(str3)
  }
  # Str2 & Str3
  if (is.null(str1Nm) & !is.null(str2Nm) & !is.null(str3Nm)){
    case<-3; cpopStrNm<-paste(str2Nm,str3Nm,sep=", "); npopStr<-2; popStrNm1<-str2Nm; popStrNm2<-str3Nm; popStr1<-str2; popStr2<-str3; npopStr1<-length(str2); npopStr2<-length(str3)
  }
  
  # Str1 & Str2 & Str3
  if (!is.null(str1Nm) & !is.null(str2Nm) & !is.null(str3Nm)){
    case<-4; cpopStrNm<-paste(str1Nm,str2Nm,str3Nm,sep=", ")
    npopStr<-3
    popStrNm1<-str1Nm; popStrNm2<-str2Nm; popStrNm3<-str3Nm
    popStr1<-str1; popStr2<-str2; popStr3<-str3
    npopStr1<-length(str1); npopStr2<-length(str2); npopStr3<-length(str3)
  }
  
  
  if (!is.null(studyName)){
    txt <- paste0("Name of the study: \"",studyName,"\"")
    txt <- paste(txt,paste0("Name of the file with the observed data: \"",obsFileNm,"\""),sep="\n")
  }else{
    txt <- paste0("Name of the file with the observed data: \"",obsFileNm,"\"")
  }
  txt <- paste(txt,paste0("Route of administration: ",adminType),sep="\n")
  if (doseType == "ss"){
    txt <- paste(txt,paste0("Dose type: steady-state with dosing interval (Tau) of",Tau),sep="\n")
  }else{
    txt <- paste(txt,"Dose type: non-steady-state",sep="\n")
  }
  txt <- paste(txt,paste0("No. of population stratification level: ",npopStr),sep="\n")
  if (case==2){
    txt <- paste(txt,paste0("Population stratification column: ",cpopStrNm),sep="\n")
    txt <- paste(txt,paste0("Population stratification ID within ",popStrNm1,": ",paste(popStr1,collapse=", ")),sep="\n")
  }else if (case==3){
    txt <- paste(txt,paste0("Population stratification columns: ",cpopStrNm),sep="\n")
    txt <- paste(txt,paste0("1st level population stratification ID within ",popStrNm1,": ",paste(popStr1,collapse=", ")),sep="\n")
    txt <- paste(txt,paste0("2nd level population stratification ID within ",popStrNm2,": ",paste(popStr2,collapse=", ")),sep="\n")
  }else if (case==4){
    txt <- paste(txt,paste0("Population stratification columns: ",cpopStrNm),sep="\n")
    txt <- paste(txt,paste0("1st level population stratification ID within ",popStrNm1,": ",paste(popStr1,collapse=", ")),sep="\n")
    txt <- paste(txt,paste0("2nd level population stratification ID within ",popStrNm2,": ",paste(popStr2,collapse=", ")),sep="\n")
    txt <- paste(txt,paste0("3rd level population stratification ID within ",popStrNm3,": ",paste(popStr3,collapse=", ")),sep="\n")
  }
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # Initiate plot lists for pdf output
  concplot   <- list(); histobsplot <- list(); popplot      <- list(); devplot <- list(); outlierplot <- list()
  forestplot <- list(); npdeplot    <- list(); histnpdeplot <- list()
  
  # calculate the NCA parameters for the observed data
  dset <- "obs"
  pddf <- data.frame()   # Summary table
  
  if(is.null(cunit)){DVleg  <- "Concentration"}else{DVleg <- paste0("Concentration (",cunit,")")}    # DV legend
  if(is.null(tunit)){IDVleg <- "Time"}else{IDVleg <- paste0("Time (",tunit,")")}                     # IDV legend
  if(is.null(dunit)){Dcol   <- "Dose"}else{Dcol <- paste0("Dose (",dunit,")")}                       # Dose
  
  # Estimate NCA metrics for case = 1
  if (case == 1){
    nca_method=1
    if(nca_method==1){
      
      ifdf <- indf
      if (nrow(ifdf) == 0){next}
      
      idd  <- unique(as.character(ifdf[,idNmObs]))
      if (is.null(doseAmtNm)){
        doseAmount <- NA
      }else{
        doseData   <- as.numeric(as.character(refdf[,doseAmtNm]))
        doseData   <- doseData[complete.cases(doseData) & doseData>0]
        doseAmount <- ifelse(length(doseData)==0, NA, paste(unique(doseData), collapse=", "))
      }
      
      # Description
      pddf  <- rbind(pddf, data.frame(a=doseAmount, b=length(idd)))
      cdata <- data.frame()
      for (i in 1:length(idd)){
        if (is.null(doseAmtNm)){
          idzAmt <- NA
        }else{
          doseData <- as.numeric(as.character(refdf[refdf[,idNmObs]==idd[i], doseAmtNm]))
          doseData <- doseData[complete.cases(doseData) & doseData>0]
          idzAmt   <- ifelse(length(doseData)==0, NA, doseData[1])
        }
        tcTI <- nca.ind.data(pkData=ifdf, ID=idd[i], dvLog = obsLog, dataType=dset,
                             idNm=idNmObs, timeNm=timeNmObs, concNm=concNmObs,
                             adminType=adminType, TI=TI,
                             dateColNm=dateColNm, dateFormat=dateFormat, timeFormat=timeFormat)
        tc      <- tcTI$tc
        iTI     <- tcTI$iTI
        if (nrow(tc)==0) next
        time    <- as.numeric(tc$time)
        conc    <- as.numeric(tc$conc)
        cdata   <- rbind(cdata,cbind(Time=time,Conc=conc,ID=idd[i]))
        NCAprm  <- est.nca(time=time,conc=conc,backExtrp=backExtrp,negConcExcl=negConcExcl,doseType=doseType,adminType=adminType,
                           doseAmt=idzAmt,method=method,AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,
                           LambdaExclude=LambdaExclude,doseTime=doseTime,Tau=Tau,TI=iTI,simFile=simFile,dset=dset,onlyNCA=onlyNCA)
        outData <- rbind(outData, data.frame(ID=idd[i],Dose=idzAmt,t(NCAprm)))
      }
    }
    if(nca_method==2){
      obsData <- tibble::as_tibble(indf)
      fullData <- tibble::as_tibble(refdf)
      
      if (nrow(obsData) == 0){next}
      
      id_name <- quo(!!rlang::sym(idNmObs))
      
      if (is.null(doseAmtNm)){
        doseAmount <- NA
        tmp <- fullData %>%  
          dplyr::select(!!id_name) %>% 
          distinct() %>% 
          dplyr::summarise(n=n())
      }else{
        dose_amt_name <- quo(!!rlang::sym(doseAmtNm))
        
        tmp <- fullData %>%  
          dplyr::select(!!dose_amt_name,!!id_name) %>% 
          dplyr::filter((!!dose_amt_name) > 0) %>% 
          dplyr::group_by(!!dose_amt_name) %>% 
          distinct() %>% 
          dplyr::summarise(n=n())
      }
        
      
      if (is.null(doseAmtNm)){
        doseAmount <- NA
      }else{
        doseData   <- as.numeric(as.character(refdf[,doseAmtNm]))
        doseData   <- doseData[complete.cases(doseData) & doseData>0]
        doseAmount <- ifelse(length(doseData)==0, NA, paste(unique(doseData), collapse=", "))
      }
      
      # Description
      pddf  <- rbind(pddf, data.frame(a=doseAmount, b=length(idd)))
      
      pddf  <- rbind(pddf, data.frame(a=200, b=24))
      indf_t <- tibble::as_tibble(indf)
      #indf_t %>% dplyr::group_by_(idNmObs) %>% dplyr::summarise(mean=mean(DV))
      outData <- indf_t %>% dplyr::group_by(ID,DOSE) %>% dplyr::mutate(conc=exp(DV)) %>% 
        dplyr::do(as.data.frame(t(est.nca(time=.$TIME,conc=.$conc,backExtrp=backExtrp,negConcExcl=negConcExcl,doseType=doseType,adminType=adminType,
                                          doseAmt=.$DOSE[1],method=method,AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,
                                          LambdaExclude=LambdaExclude,doseTime=doseTime,Tau=Tau,TI=iTI,simFile=simFile,dset=dset,onlyNCA=onlyNCA))))
      outData <- as.data.frame(outData)
      
    }
    
    if(nca_method==3){
      
      ifdf <- indf
      if (nrow(ifdf) == 0){next}
      
      idd  <- unique(as.character(ifdf[,idNmObs]))
      if (is.null(doseAmtNm)){
        doseAmount <- NA
      }else{
        doseData   <- as.numeric(as.character(refdf[,doseAmtNm]))
        doseData   <- doseData[complete.cases(doseData) & doseData>0]
        doseAmount <- ifelse(length(doseData)==0, NA, paste(unique(doseData), collapse=", "))
      }
      
      # Description
      pddf  <- rbind(pddf, data.frame(a=doseAmount, b=length(idd)))
      cdata <- data.frame()
      nca_calc <- function(id,data) {
        if (is.null(doseAmtNm)){
          idzAmt <- NA
        }else{
          doseData <- as.numeric(as.character(data[data[,idNmObs]==id, doseAmtNm]))
          doseData <- doseData[complete.cases(doseData) & doseData>0]
          idzAmt   <- ifelse(length(doseData)==0, NA, doseData[1])
        }
        tcTI <- nca.ind.data(pkData=data, ID=id, dvLog = obsLog, dataType=dset,
                             idNm=idNmObs, timeNm=timeNmObs, concNm=concNmObs,
                             adminType=adminType, TI=TI,
                             dateColNm=dateColNm, dateFormat=dateFormat, timeFormat=timeFormat)
        tc      <- tcTI$tc
        iTI     <- tcTI$iTI
        if (nrow(tc)==0) next
        time    <- as.numeric(tc$time)
        conc    <- as.numeric(tc$conc)
        cdata   <- rbind(cdata,cbind(Time=time,Conc=conc,ID=id))
        NCAprm  <- est.nca(time=time,conc=conc,backExtrp=backExtrp,negConcExcl=negConcExcl,doseType=doseType,adminType=adminType,
                           doseAmt=idzAmt,method=method,AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,
                           LambdaExclude=LambdaExclude,doseTime=doseTime,Tau=Tau,TI=iTI,simFile=simFile,dset=dset,onlyNCA=onlyNCA)
        outData <- data.frame(ID=id,Dose=idzAmt,t(NCAprm))
        return(outData)
      }
      parallel=T
      if(parallel){
        parallel <- PopED::start_parallel(parallel) 
        on.exit(if(parallel && (attr(parallel,"type")=="snow")) parallel::stopCluster(attr(parallel,"cluster")))
      }  
      if(parallel && (attr(parallel,"type")=="multicore")){
        res <- parallel::mclapply(idd,nca_calc,mc.cores=attr(parallel, "cores"),data=ifdf)
      } else if(parallel && (attr(parallel,"type")=="snow")){
        res <- parallel::parLapply(attr(parallel, "cluster"),idd,nca_calc,data=ifdf)
      } else {
        res <- lapply(idd,nca_calc,data=ifdf)
      }
      outData <- do.call(rbind,res)
    }
    
    if (!noPlot){
      figlbl <- "All-data"
      gdr    <- dv.plot(df=cdata,xvar="Time",yvar="Conc",myXlab=IDVleg,myYlab=DVleg,color="ID")
      suppressMessages(suppressWarnings(grid.arrange(gdr)))
      ggr <- grid.grab()
      concplot[[length(concplot)+1]] <- ggr
      if (printOut){
        fl <- paste0(usrdir,"/TimeConc_",figlbl,".",figFormat)
        if (figFormat=="tiff"){
          eval(parse(text=paste0(figFormat,"(file=\"",fl,"\",height=10,width=15,units=\"cm\",compression=\"lzw\",res=120)")))
        }else{
          eval(parse(text=paste0(figFormat,"(file=\"",fl,"\",height=10,width=15,units=\"cm\",res=120)")))
        }
        suppressMessages(suppressWarnings(grid.arrange(gdr)))
        dev.off()
      }
      
      # Obs hist plot
      if (nrow(outData)>=5){
        plotData <- subset(outData, select=c(AUClast,AUCINF_obs,Cmax,Tmax))
        numPrm   <- sapply(plotData, FUN=function(x){x <- as.numeric(as.character(x)); length(x[complete.cases(x)])})
        if (length(numPrm[numPrm>=5]) == 0) next
        
        pltPrm      <- names(numPrm[numPrm>=5])
        figlbl      <- NULL
        histobsgrob <- histobs.plot(plotData=plotData,figlbl=figlbl,param=pltPrm,cunit=cunit,tunit=tunit,spread=spread)
        gdr         <- histobsgrob$gdr
        mylegend    <- histobsgrob$legend
        lheight     <- histobsgrob$lheight
        if (printOut){
          fl <- paste0(usrdir,"/HistObs")
          if (figFormat=="tiff"){
            eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=20,width=15,units=\"cm\",compression=\"lzw\",res=120)")))
          }else{
            eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=20,width=15,units=\"cm\",res=120)")))
          }
          suppressMessages(suppressWarnings(grid.arrange(gdr, mylegend, heights = unit.c(unit(1,"npc")-lheight, lheight))))
          dev.off()
        }
        suppressMessages(suppressWarnings(grid.arrange(gdr, mylegend, heights = unit.c(unit(1,"npc")-lheight, lheight))))
        ggr <- grid.grab()
        histobsplot[[length(histobsplot)+1]] <- ggr
      }
    }
    cnm         <- c(Dcol,"No. of individuals")
    names(pddf) <- cnm
  }
  
  # Estimate NCA metrics for case = 2
  if (case == 2){
    for (s1 in 1:npopStr1){
      cdata <- data.frame()
      ifdf  <- indf[indf[,popStrNm1]==popStr1[s1],]
      if (nrow(ifdf) == 0){next}
      idd <- unique(as.character(ifdf[,idNmObs]))
      if (!is.null(doseAmtNm)){
        doseData   <- as.numeric(as.character(refdf[refdf[,popStrNm1]==popStr1[s1], doseAmtNm]))
        doseData   <- doseData[complete.cases(doseData) & doseData>0]
        doseAmount <- ifelse(length(doseData)==0, NA, paste(unique(doseData), collapse=", "))
      }else{
        doseAmount <- NA
      }
      # Description
      pddf <- rbind(pddf, data.frame(a=popStr1[s1], b=doseAmount, c=length(idd)))
      for (i in 1:length(idd)){
        if (!is.null(doseAmtNm)){
          doseData <- as.numeric(as.character(refdf[refdf[,popStrNm1]==popStr1[s1] & refdf[,idNmObs]==idd[i], doseAmtNm]))
          doseData <- doseData[complete.cases(doseData) & doseData>0]
          idzAmt   <- ifelse(length(doseData)==0, NA, doseData[1])
        }else{
          idzAmt <- NA
        }
        tcTI <- nca.ind.data(pkData=ifdf, ID=idd[i], dvLog = obsLog, dataType=dset,
                             idNm=idNmObs, timeNm=timeNmObs, concNm=concNmObs,
                             adminType=adminType, TI=TI,
                             dateColNm=dateColNm, dateFormat=dateFormat, timeFormat=timeFormat)
        tc      <- tcTI$tc
        iTI     <- tcTI$iTI
        if (nrow(tc)==0) next
        time    <- as.numeric(tc$time)
        conc    <- as.numeric(tc$conc)
        cdata   <- rbind(cdata,cbind(Time=time,Conc=conc,ID=idd[i]))
        NCAprm  <- est.nca(time=time,conc=conc,backExtrp=backExtrp,negConcExcl=negConcExcl,
                           doseType=doseType,adminType=adminType,doseAmt=idzAmt,method=method,
                           AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,LambdaExclude=LambdaExclude,
                           doseTime=doseTime,Tau=Tau,TI=iTI,simFile=simFile,dset=dset,onlyNCA=onlyNCA)
        outData <- rbind(outData, data.frame(ID=idd[i],STRAT1=popStr1[s1],Dose=idzAmt,t(NCAprm)))
      }
      
      if (!noPlot){
        figlbl <- paste0(popStrNm1,"-",popStr1[s1])
        gdr    <- dv.plot(df=cdata,xvar="Time",yvar="Conc",myXlab=IDVleg,myYlab=DVleg,color="ID")
        suppressMessages(suppressWarnings(grid.arrange(gdr)))
        ggr <- grid.grab()
        concplot[[length(concplot)+1]] <- ggr
        if (printOut){
          fl <- paste0(usrdir,"/TimeConc_",figlbl,".",figFormat)
          if (figFormat=="tiff"){
            eval(parse(text=paste0(figFormat,"(file=\"",fl,"\",height=10,width=15,units=\"cm\",compression=\"lzw\",res=120)")))
          }else{
            eval(parse(text=paste0(figFormat,"(file=\"",fl,"\",height=10,width=15,units=\"cm\",res=120)")))
          }
          suppressMessages(suppressWarnings(grid.arrange(gdr)))
          dev.off()
        }
        
        # Obs hist plot
        plotData <- subset(outData, STRAT1==popStr1[s1], select=c(AUClast,AUCINF_obs,Cmax,Tmax))
        if (nrow(plotData)<5) next
        numPrm   <- sapply(plotData, FUN=function(x){x <- as.numeric(as.character(x)); length(x[complete.cases(x)])})
        if (length(numPrm[numPrm>=5]) == 0) next
        
        pltPrm      <- names(numPrm[numPrm>=5])
        figlbl      <- paste0(popStrNm1,"-",popStr1[s1])
        histobsgrob <- histobs.plot(plotData=plotData,figlbl=figlbl,param=pltPrm,cunit=cunit,tunit=tunit,spread=spread)
        gdr         <- histobsgrob$gdr
        mylegend    <- histobsgrob$legend
        lheight     <- histobsgrob$lheight
        if (printOut){
          fl <- paste0(usrdir,"/HistObs_",figlbl)
          if (figFormat=="tiff"){
            eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=20,width=15,units=\"cm\",compression=\"lzw\",res=120)")))
          }else{
            eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=20,width=15,units=\"cm\",res=120)")))
          }
          suppressMessages(suppressWarnings(grid.arrange(gdr, mylegend, heights = unit.c(unit(1,"npc")-lheight, lheight))))
          dev.off()
        }
        suppressMessages(suppressWarnings(grid.arrange(gdr, mylegend, heights = unit.c(unit(1,"npc")-lheight, lheight))))
        ggr <- grid.grab()
        histobsplot[[length(histobsplot)+1]] <- ggr
      }
    }
    cnm         <- c(popStrNm1,Dcol,"No. of individuals")
    names(pddf) <- cnm
  }
  
  # Estimate NCA metrics for case = 3
  if (case == 3){
    for (s1 in 1:npopStr1){
      for (s2 in 1:npopStr2){
        cdata <- data.frame()
        ifdf <- indf[indf[,popStrNm1]==popStr1[s1] & indf[,popStrNm2]==popStr2[s2],]
        if (nrow(ifdf) == 0){next}
        idd <- unique(as.character(ifdf[,idNmObs]))
        if (!is.null(doseAmtNm)){
          doseData   <- as.numeric(as.character(refdf[refdf[,popStrNm1]==popStr1[s1] & refdf[,popStrNm2]==popStr2[s2], doseAmtNm]))
          doseData   <- doseData[complete.cases(doseData) & doseData>0]
          doseAmount <- ifelse(length(doseData)==0, NA, paste(unique(doseData), collapse=", "))
        }else{
          doseAmount <- NA
        }
        # Description
        pddf   <- rbind(pddf, data.frame(a=popStr1[s1], b=popStr2[s2], c=doseAmount, d=length(idd)))
        for (i in 1:length(idd)){
          if (!is.null(doseAmtNm)){
            doseData <- as.numeric(as.character(refdf[refdf[,popStrNm1]==popStr1[s1] & refdf[,popStrNm2]==popStr2[s2] & refdf[,idNmObs]==idd[i], doseAmtNm]))
            doseData <- doseData[complete.cases(doseData) & doseData>0]
            idzAmt   <- ifelse(length(doseData)==0, NA, doseData[1])
          }else{
            idzAmt <- NA
          }
          tcTI <- nca.ind.data(pkData=ifdf, ID=idd[i], dvLog = obsLog, dataType=dset,
                               idNm=idNmObs, timeNm=timeNmObs, concNm=concNmObs,
                               adminType=adminType, TI=TI,
                               dateColNm=dateColNm, dateFormat=dateFormat, timeFormat=timeFormat)
          tc      <- tcTI$tc
          iTI     <- tcTI$iTI
          if (nrow(tc)==0) next
          time    <- as.numeric(tc$time)
          conc    <- as.numeric(tc$conc)
          cdata   <- rbind(cdata,cbind(Time=time,Conc=conc,ID=idd[i]))
          NCAprm  <- est.nca(time=time,conc=conc,backExtrp=backExtrp,negConcExcl=negConcExcl,
                             doseType=doseType,adminType=adminType,doseAmt=idzAmt,method=method,
                             AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,LambdaExclude=LambdaExclude,
                             doseTime=doseTime,Tau=Tau,TI=iTI,simFile=simFile,dset=dset,onlyNCA=onlyNCA)
          outData <- rbind(outData, data.frame(ID=idd[i],STRAT1=popStr1[s1],STRAT2=popStr2[s2],Dose=idzAmt,t(NCAprm)))
        }
        
        if (!noPlot){
          figlbl <- paste0(popStrNm1,"-",popStr1[s1],"_",popStrNm2,"-",popStr2[s2])
          gdr    <- dv.plot(df=cdata,xvar="Time",yvar="Conc",myXlab=IDVleg,myYlab=DVleg,color="ID")
          suppressMessages(suppressWarnings(grid.arrange(gdr)))
          ggr <- grid.grab()
          concplot[[length(concplot)+1]] <- ggr
          
          if (printOut){
            fl <- paste0(usrdir,"/TimeConc_",figlbl,".",figFormat)
            if (figFormat=="tiff"){
              eval(parse(text=paste0(figFormat,"(file=\"",fl,"\",height=10,width=15,units=\"cm\",compression=\"lzw\",res=120)")))
            }else{
              eval(parse(text=paste0(figFormat,"(file=\"",fl,"\",height=10,width=15,units=\"cm\",res=120)")))
            }
            suppressMessages(suppressWarnings(grid.arrange(gdr)))
            dev.off()
          }
          
          # Obs hist plot
          plotData <- subset(outData, STRAT1==popStr1[s1] & STRAT2==popStr2[s2], select=c(AUClast,AUCINF_obs,Cmax,Tmax))
          if (nrow(plotData)<5) next
          numPrm   <- sapply(plotData, FUN=function(x){x <- as.numeric(as.character(x)); length(x[complete.cases(x)])})
          if (length(numPrm[numPrm>=5]) == 0) next
          
          pltPrm      <- names(numPrm[numPrm>=5])
          figlbl      <- paste0(popStrNm1,"-",popStr1[s1],"_",popStrNm2,"-",popStr2[s2])
          histobsgrob <- histobs.plot(plotData=plotData,figlbl=figlbl,param=pltPrm,cunit=cunit,tunit=tunit,spread=spread)
          gdr         <- histobsgrob$gdr
          mylegend    <- histobsgrob$legend
          lheight     <- histobsgrob$lheight
          if (printOut){
            fl <- paste0(usrdir,"/HistObs_",figlbl)
            if (figFormat=="tiff"){
              eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=20,width=15,units=\"cm\",compression=\"lzw\",res=120)")))
            }else{
              eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=20,width=15,units=\"cm\",res=120)")))
            }
            suppressMessages(suppressWarnings(grid.arrange(gdr, mylegend, heights = unit.c(unit(1,"npc")-lheight, lheight))))
            dev.off()
          }
          suppressMessages(suppressWarnings(grid.arrange(gdr, mylegend, heights = unit.c(unit(1,"npc")-lheight, lheight))))
          ggr <- grid.grab()
          histobsplot[[length(histobsplot)+1]] <- ggr
        }
      }
    }
    cnm         <- c(popStrNm1,popStrNm2,Dcol,"No. of individuals")
    names(pddf) <- cnm
  }
  
  # Estimate NCA metrics for case = 4
  if (case == 4){
    for (s1 in 1:npopStr1){
      for (s2 in 1:npopStr2){
        for (s3 in 1:npopStr3){
          cdata <- data.frame()
          ifdf <- indf[(indf[,popStrNm1]==popStr1[s1] & indf[,popStrNm2]==popStr2[s2] & indf[,popStrNm3]==popStr3[s3]),]
          if (nrow(ifdf) == 0){next}
          idd <- unique(as.character(ifdf[,idNmObs]))
          if (!is.null(doseAmtNm)){
            doseData   <- as.numeric(as.character(refdf[refdf[,popStrNm1]==popStr1[s1] & refdf[,popStrNm2]==popStr2[s2] & refdf[,popStrNm3]==popStr3[s3], doseAmtNm]))
            doseData   <- doseData[complete.cases(doseData) & doseData>0]
            doseAmount <- ifelse(length(doseData)==0, NA, paste(unique(doseData), collapse=", "))
          }else{
            doseAmount <- NA
          }
          # Description
          pddf  <- rbind(pddf, data.frame(a=popStr1[s1], b=popStr2[s2], c=popStr3[s3], d=doseAmount, e=length(idd)))
          for (i in 1:length(idd)){
            if (!is.null(doseAmtNm)){
              doseData <- as.numeric(as.character(refdf[refdf[,popStrNm1]==popStr1[s1] & refdf[,popStrNm2]==popStr2[s2] & refdf[,popStrNm3]==popStr3[s3] & refdf[,idNmObs]==idd[i], doseAmtNm]))
              doseData <- doseData[complete.cases(doseData) & doseData>0]
              idzAmt   <- ifelse(length(doseData)==0, NA, doseData[1])
            }else{
              idzAmt <- NA
            }
            tcTI <- nca.ind.data(pkData=ifdf, ID=idd[i], dvLog = obsLog, dataType=dset,
                                 idNm=idNmObs, timeNm=timeNmObs, concNm=concNmObs,
                                 adminType=adminType, TI=TI,
                                 dateColNm=dateColNm, dateFormat=dateFormat, timeFormat=timeFormat)
            tc      <- tcTI$tc
            iTI     <- tcTI$iTI
            if (nrow(tc)==0) next
            time    <- as.numeric(tc$time)
            conc    <- as.numeric(tc$conc)
            cdata   <- rbind(cdata,cbind(Time=time,Conc=conc,ID=idd[i]))
            NCAprm  <- est.nca(time=time,conc=conc,backExtrp=backExtrp,negConcExcl=negConcExcl,
                               doseType=doseType,adminType=adminType,doseAmt=idzAmt,method=method,
                               AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,LambdaExclude=LambdaExclude,
                               doseTime=doseTime,Tau=Tau,TI=iTI,simFile=simFile,dset=dset,onlyNCA=onlyNCA)
            outData <- rbind(outData, data.frame(ID=idd[i],STRAT1=popStr1[s1],STRAT2=popStr2[s2],STRAT3=popStr3[s3],Dose=idzAmt,t(NCAprm)))
          }
          
          if (!noPlot){
            figlbl <- paste0(popStrNm1,"-",popStr1[s1],"_",popStrNm2,"-",popStr2[s2],"_",popStrNm3,"-",popStr3[s3])
            gdr    <- dv.plot(df=cdata,xvar="Time",yvar="Conc",myXlab=IDVleg,myYlab=DVleg,color="ID")
            suppressMessages(suppressWarnings(grid.arrange(gdr)))
            ggr <- grid.grab()
            concplot[[length(concplot)+1]] <- ggr
            if (printOut){
              fl <- paste0(usrdir,"/TimeConc_",figlbl,".",figFormat)
              if (figFormat=="tiff"){
                eval(parse(text=paste0(figFormat,"(file=\"",fl,"\",height=10,width=15,units=\"cm\",compression=\"lzw\",res=120)")))
              }else{
                eval(parse(text=paste0(figFormat,"(file=\"",fl,"\",height=10,width=15,units=\"cm\",res=120)")))
              }
              suppressMessages(suppressWarnings(grid.arrange(gdr)))
              dev.off()
            }
            
            # Obs hist plot
            plotData <- subset(outData, STRAT1==popStr1[s1] & STRAT2==popStr2[s2] & STRAT3==popStr3[s3], select=c(AUClast,AUCINF_obs,Cmax,Tmax))
            if (nrow(plotData)<5) next
            numPrm   <- sapply(plotData, FUN=function(x){x <- as.numeric(as.character(x)); length(x[complete.cases(x)])})
            if (length(numPrm[numPrm>=5]) == 0) next
            
            pltPrm      <- names(numPrm[numPrm>=5])
            figlbl      <- paste0(popStrNm1,"-",popStr1[s1],"_",popStrNm2,"-",popStr2[s2],"_",popStrNm3,"-",popStr3[s3])
            histobsgrob <- histobs.plot(plotData=plotData,figlbl=figlbl,param=pltPrm,cunit=cunit,tunit=tunit,spread=spread)
            gdr         <- histobsgrob$gdr
            mylegend    <- histobsgrob$legend
            lheight     <- histobsgrob$lheight
            if (printOut){
              fl <- paste0(usrdir,"/HistObs_",figlbl)
              if (figFormat=="tiff"){
                eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=20,width=15,units=\"cm\",compression=\"lzw\",res=120)")))
              }else{
                eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=20,width=15,units=\"cm\",res=120)")))
              }
              suppressMessages(suppressWarnings(grid.arrange(gdr, mylegend, heights = unit.c(unit(1,"npc")-lheight, lheight))))
              dev.off()
            }
            suppressMessages(suppressWarnings(grid.arrange(gdr, mylegend, heights = unit.c(unit(1,"npc")-lheight, lheight))))
            ggr <- grid.grab()
            histobsplot[[length(histobsplot)+1]] <- ggr
          }
        }
      }
    }
    cnm         <- c(popStrNm1,popStrNm2,popStrNm3,Dcol,"No. of individuals")
    names(pddf) <- cnm
  }
  
  # Statistical analysis for each patient group
  statData <- data.frame()
  statPrm  <- c("N","Nunique","Min","Max","Mean","Median","SD","SE","CVp","CI95l","CI95u","geoMean","gCVp")
  ncaPrm   <- c("C0","Tmax","Cmax","Cmax_D","Tlast","Clast","AUClast","AUMClast","MRTlast","No_points_Lambda_z",
                "AUC_pBack_Ext_obs","AUC_pBack_Ext_pred","AUClower_upper","Rsq","Rsq_adjusted","Corr_XY","Lambda_z",
                "HL_Lambda_z","AUCINF_obs","AUCINF_D_obs","AUC_pExtrap_obs","Vz_obs","Cl_obs","AUCINF_pred","AUCINF_D_pred",
                "AUC_pExtrap_pred","Vz_pred","Cl_pred","AUMCINF_obs","AUMC_pExtrap_obs","AUMCINF_pred","AUMC_pExtrap_pred",
                "MRTINF_obs","MRTINF_pred","Vss_obs","Vss_pred","Tau","Tmin","Cmin","Cavg","p_Fluctuation","Accumulation_Index","Clss")
  
  if (case == 1){
    tmpDF    <- outData[,ncaPrm]
    statData <- sapply(tmpDF, function(x){x<-as.numeric(as.character(x)); if(length(x[complete.cases(x)])<2){rep(NA,13)}else{out.digits(calc.stat(x[complete.cases(x)]),dig=4)}})
    statData <- data.frame(statData)
    if (nrow(statData) == 0) next
    statData <- cbind(data.frame(Stat=statPrm), statData)
  }
  if (case == 2){
    for (s1 in 1:npopStr1){
      tmpDF    <- outData[outData$STRAT1==popStr1[s1],ncaPrm]
      tmpStat  <- sapply(tmpDF, function(x){x<-as.numeric(as.character(x)); if(length(x[complete.cases(x)])<2){rep(NA,13)}else{out.digits(calc.stat(x[complete.cases(x)]),dig=4)}})
      tmpStat  <- data.frame(tmpStat)
      if (nrow(tmpStat) == 0) next
      tmpStat  <- cbind(data.frame(STRAT1=popStr1[s1],Stat=statPrm), tmpStat)
      statData <- rbind(statData, tmpStat)
    }
    names(statData)[names(statData)%in%"STRAT1"] <- popStrNm1
  }
  if (case == 3){
    for (s1 in 1:npopStr1){
      for (s2 in 1:npopStr2){
        tmpDF    <- outData[outData$STRAT1==popStr1[s1] & outData$STRAT2==popStr2[s2],ncaPrm]
        tmpStat  <- sapply(tmpDF, function(x){x<-as.numeric(as.character(x)); if(length(x[complete.cases(x)])<2){rep(NA,13)}else{out.digits(calc.stat(x[complete.cases(x)]),dig=4)}})
        tmpStat  <- data.frame(tmpStat)
        if (nrow(tmpStat) == 0) next
        tmpStat  <- cbind(data.frame(STRAT1=popStr1[s1],STRAT2=popStr2[s2],Stat=statPrm), tmpStat)
        statData <- rbind(statData, tmpStat)
      }
    }
    names(statData)[names(statData)%in%c("STRAT1","STRAT2")] <- c(popStrNm1,popStrNm2)
  }
  if (case == 4){
    for (s1 in 1:npopStr1){
      for (s2 in 1:npopStr2){
        for (s3 in 1:npopStr3){
          tmpDF    <- outData[outData$STRAT1==popStr1[s1] & outData$STRAT2==popStr2[s2] & outData$STRAT3==popStr3[s3],ncaPrm]
          tmpStat  <- sapply(tmpDF, function(x){x<-as.numeric(as.character(x)); if(length(x[complete.cases(x)])<2){rep(NA,13)}else{out.digits(calc.stat(x[complete.cases(x)]),dig=4)}})
          tmpStat  <- data.frame(tmpStat)
          if (nrow(tmpStat) == 0) next
          tmpStat  <- cbind(data.frame(STRAT1=popStr1[s1],STRAT2=popStr2[s2],STRAT3=popStr3[s3],Stat=statPrm), tmpStat)
          statData <- rbind(statData, tmpStat)
        }
      }
    }
    names(statData)[names(statData)%in%c("STRAT1","STRAT2","STRAT3")] <- c(popStrNm1,popStrNm2,popStrNm3)
  }
  
  if (printOut){
    if (is.null(outFileNm) || outFileNm==""){
      write.table(statData, file=paste0(usrdir,"/ObsStat.tsv"), sep="\t", col.names=T, row.names=F, quote=F)
    }else{
      write.table(statData, file=paste0(usrdir,"/ObsStat-",outFileNm,".tsv"), sep="\t", col.names=T, row.names=F, quote=F)
    }
  }
  
  
  # Print output files if simulated data does not exist
  if (is.null(simFile) || onlyNCA){
    # Raname ID and stratifier columns and format output table sigfig
    if (case == 1){
      names(outData)[names(outData)%in%c("ID")] <- c(idNmObs)
      outData[,c(2:ncol(outData))] <- as.data.frame(lapply(outData[,c(2:ncol(outData))], FUN=function(x) out.digits(x,dig=4)))
    }
    if (case == 2){
      names(outData)[names(outData)%in%c("ID","STRAT1")] <- c(idNmObs,popStrNm1)
      outData[,c(3:ncol(outData))] <- as.data.frame(lapply(outData[,c(3:ncol(outData))], FUN=function(x) out.digits(x,dig=4)))
    }
    if (case == 3){
      names(outData)[names(outData)%in%c("ID","STRAT1","STRAT2")] <- c(idNmObs,popStrNm1,popStrNm2)
      outData[,c(4:ncol(outData))] <- as.data.frame(lapply(outData[,c(4:ncol(outData))], FUN=function(x) out.digits(x,dig=4)))
    }
    if (case == 4){
      names(outData)[names(outData)%in%c("ID","STRAT1","STRAT2","STRAT3")] <- c(idNmObs,popStrNm1,popStrNm2,popStrNm3)
      outData[,c(5:ncol(outData))] <- as.data.frame(lapply(outData[,c(5:ncol(outData))], FUN=function(x) out.digits(x,dig=4)))
    }
    
    # Subset table to print in the report
    if (case == 1) {prnTab <- head(cbind(outData[,1:2], subset(outData, select = tabCol)),50); names(prnTab)[1:2] <- names(outData)[1:2]}
    if (case == 2) {prnTab <- head(cbind(outData[,1:3], subset(outData, select = tabCol)),50); names(prnTab)[1:3] <- names(outData)[1:3]}
    if (case == 3) {prnTab <- head(cbind(outData[,1:4], subset(outData, select = tabCol)),50); names(prnTab)[1:4] <- names(outData)[1:4]}
    if (case == 4) {prnTab <- head(cbind(outData[,1:5], subset(outData, select = tabCol)),50); names(prnTab)[1:5] <- names(outData)[1:5]}
    
    # Add unit to report table header
    tabUnit1 <- data.frame(NAME=c("Dose"))
    
    tabUnit2 <- data.frame(NAME=c("C0","Tmax","Cmax","Cmax_D","Tlast","Clast","AUClast","AUMClast","MRTlast","AUClower_upper","Lambda_z","Lambda_z_lower","Lambda_z_upper","HL_Lambda_z","AUCINF_obs","AUCINF_D_obs","Vz_obs","Cl_obs","AUCINF_pred","AUCINF_D_pred","Vz_pred","Cl_pred","AUMCINF_obs","AUMCINF_pred","MRTINF_obs","MRTINF_pred","Tau","Tmin","Cmin","Cavg","AUCtau","AUMCtau","Clss","Vss_obs","Vss_pred"))
    tabUnit <- rbind(tabUnit1,tabUnit2)
    
    streamsEnv <- parent.frame()
    if (exists("outData"))     assign("ncaOutput", outData, envir=streamsEnv)
    if (exists("statData"))    assign("ObsStat", statData, envir=streamsEnv)
    if (exists("concplot"))    assign("concPlot", concplot, envir=streamsEnv)
    if (exists("histobsplot")) assign("histobsPlot", histobsplot, envir=streamsEnv)
    
    if (printOut){
      if (is.null(outFileNm) || outFileNm==""){
        write.table(outData, file=paste0(usrdir,"/ncaOutput.tsv"), sep="\t", row.names=F, col.names=T, quote=F)
      }else{
        write.table(outData, file=paste0(usrdir,"/ncaOutput-",outFileNm,".tsv"), sep="\t", row.names=F, col.names=T, quote=F)
      }
      fnOut <- list(arglist=match.call(), case=case, TXT=txt, pddf=pddf, prnTab=prnTab, spread=spread, conc=concplot, histobs=histobsplot)   # Function output list
    }    
  }else{
    ########################################################################
    ############## Analyze the simulated data if exists ####################
    ########################################################################
    simDir <- paste0(usrdir,"/SIMDATA")
    
    # Check status of SIMDATA
    if (file.exists(simDir)){
      if (is.null(overwrite_SIMDATA)){
        if(interactive()){
          cat("\nDirectory \"SIMDATA\" already exists.\n") 
          cat("Please choose one of the following options:\n")
          dirTest <- menu(c("Overwrite it",
                            "Rename \"SIMDATA\" and create a new one",
                            "Use \"SIMDATA\" as it is."))
        } else {
          dirTest <- "1"
        }
      }else if (overwrite_SIMDATA){
        dirTest <- "1"
      }else if (!overwrite_SIMDATA) {
        dirTest <- "3"
      }
      
      if (dirTest == "1"){
        unlink(simDir, recursive=TRUE)
        dir.create(simDir)
      }else if (dirTest == "2"){
        print("\nRenaming \"SIMDATA\" to \"SIMDATA_PREVIOUS\"\n")
        file.rename(from="SIMDATA", to="SIMDATA_PREVIOUS")
        dir.create(simDir)
      }else if (dirTest == "3"){
        file_list <- list.files(path="./SIMDATA/", pattern="sim_[0-9]*.csv", full.names=T)
      }else{
        setwd(usrdir);stop("Don't know what to do with \"SIMDATA\" \n")
        #setwd(usrdir);stop("Bad choice!!! Please choose either 1 or 2 or 3\n")
      }
    }else{
      dir.create(simDir)
      dirTest <- "0"
    }
    
    # If SIMDATA is newly created, read NM output and perform NCA
    if (dirTest != 3){
      if (is.data.frame(simFile)){
        nmdf <- simFile
      }else{
        if (new_data_method){
          nmdf <- read_nm_table(simFile,sim_num = T,sim_name="NSUB")
          nmdf <- data.frame(nmdf)
        } else {
          nmdf <- nca.read.sim(simFile=simFile, MDV.rm=F)
        }
      }
      if (printOut){
        if (is.null(outFileNm) || outFileNm==""){
          #write.table(nmdf, file=paste0(usrdir,"/ncaSimData.tsv"), row.names=F, quote=F, sep="\t")
          readr::write_delim(nmdf, file.path(usrdir,"ncaSimData.tsv.gz"), delim = "\t")
        }else{
          #write.table(nmdf, file=paste0(usrdir,"/ncaSimData-",outFileNm,".tsv"), row.names=F, quote=F, sep="\t")
          readr::write_delim(nmdf, file.path(usrdir,paste0("ncaSimData-",outFileNm,".tsv.gz")), delim = "\t")
        }
      }
      
      simID <- unique(sort(nmdf$NSUB))
      nsim  <- length(simID)
      dset  <- "sim"
      
      # Perform checks on simulated NM output file (nmdf)
      simList <- nca.check.sim(simData=nmdf,
                               idNmSim=idNmSim,timeNmSim=timeNmSim,concNmSim=concNmSim,
                               filterNm=filterNm,filterExcl=filterExcl,
                               str1Nm=str1Nm,str1=str1,
                               str2Nm=str2Nm,str2=str2,
                               str3Nm=str3Nm,str3=str3,
                               adminType=adminType,TI=TI,doseAmtNm=doseAmtNm,
                               blqNm=blqNm,blqExcl=blqExcl,
                               evid=evid,evidIncl=evidIncl,mdv=mdv)
      
      nmdf      <- simList$simData
      srdf      <- simList$srdf
      str1      <- simList$str1
      str2      <- simList$str2
      str3      <- simList$str3
      doseAmtNm <- simList$doseAmtNm
      
      
      # cpopStrNm = Combined stratifyting column names
      # npopStr   = Number of stratification levels
      # popStrNm1 = 1st level stratifying column name
      # popStr1   = 1st level stratification ID names
      # npopStr1  = Number of 1st level stratification ID names
      if (!is.null(str1Nm) & is.null(str2Nm)  & is.null(str3Nm)) {case<-2; cpopStrNm<-str1Nm; npopStr<-1; popStrNm1<-str1Nm; popStr1<-str1; npopStr1<-length(str1)} # Str1
      if (is.null(str1Nm)  & !is.null(str2Nm) & is.null(str3Nm)) {case<-2; cpopStrNm<-str2Nm; npopStr<-1; popStrNm1<-str2Nm; popStr1<-str2; npopStr1<-length(str2)} # Str2
      if (is.null(str1Nm)  & is.null(str2Nm)  & !is.null(str3Nm)){case<-2; cpopStrNm<-str3Nm; npopStr<-1; popStrNm1<-str3Nm; popStr1<-str3; npopStr1<-length(str3)} # Str3
      
      # Str1 & Str2
      if (!is.null(str1Nm) & !is.null(str2Nm) & is.null(str3Nm)){
        case<-3; cpopStrNm<-paste(str1Nm,str2Nm,sep=", "); npopStr<-2; popStrNm1<-str1Nm; popStrNm2<-str2Nm; popStr1<-str1; popStr2<-str2; npopStr1<-length(str1); npopStr2<-length(str2)
      }
      # Str1 & Str3
      if (!is.null(str1Nm) & is.null(str2Nm) & !is.null(str3Nm)){
        case<-3; cpopStrNm<-paste(str1Nm,str3Nm,sep=", "); npopStr<-2; popStrNm1<-str1Nm; popStrNm2<-str3Nm; popStr1<-str1; popStr2<-str3; npopStr1<-length(str1); npopStr2<-length(str3)
      }
      # Str2 & Str3
      if (is.null(str1Nm) & !is.null(str2Nm) & !is.null(str3Nm)){
        case<-3; cpopStrNm<-paste(str2Nm,str3Nm,sep=", "); npopStr<-2; popStrNm1<-str2Nm; popStrNm2<-str3Nm; popStr1<-str2; popStr2<-str3; npopStr1<-length(str2); npopStr2<-length(str3)
      }
      # Str1 & Str2 & Str3
      if (!is.null(str1Nm) & !is.null(str2Nm) & !is.null(str3Nm)){
        case<-4; cpopStrNm<-paste(str1Nm,str2Nm,str3Nm,sep=", ")
        npopStr<-3
        popStrNm1<-str1Nm; popStrNm2<-str2Nm; popStrNm3<-str3Nm
        popStr1<-str1; popStr2<-str2; popStr3<-str3
        npopStr1<-length(str1); npopStr2<-length(str2); npopStr3<-length(str3)
      }
      
      
      # Calculate AUC parameters for the simulation output
      for (s in 1:nsim){
        simData <- data.frame()
        smdf    <- nmdf[nmdf$NSUB == simID[s],]
        # Calculate NCA parameters
        if (case == 1){
          ifdf <- smdf
          if (nrow(ifdf) == 0){next}
          idd  <- unique(as.character(ifdf[,idNmSim]))
          for (i in 1:length(idd)){
            if (!is.null(doseAmtNm)){
              doseData <- as.numeric(as.character(srdf[srdf[,idNmSim]==idd[i], doseAmtNm]))
              doseData <- doseData[complete.cases(doseData) & doseData>0]
              idzAmt   <- ifelse(length(doseData)==0, NA, doseData[1])
            }else{
              idzAmt <- NA
            }
            stcTI <- nca.ind.data(pkData=ifdf, ID=idd[i], dvLog = simLog, dataType=dset,
                                  idNm=idNmSim, timeNm=timeNmSim, concNm=concNmSim,
                                  adminType=adminType, TI=TI,
                                  dateColNm=dateColNm, dateFormat=dateFormat, timeFormat=timeFormat)
            stc     <- stcTI$tc
            iTI     <- stcTI$iTI
            if (nrow(stc)==0) next
            time    <- as.numeric(stc$time)
            conc    <- as.numeric(stc$conc)
            NCAprm  <- est.nca(time=time,conc=conc,backExtrp=backExtrp,negConcExcl=negConcExcl,
                               doseType=doseType,adminType=adminType,doseAmt=idzAmt,method=method,
                               AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,LambdaExclude=LambdaExclude,
                               doseTime=doseTime,Tau=Tau,TI=iTI,simFile=simFile,dset=dset,onlyNCA=onlyNCA)
            simData <- rbind(simData, data.frame(ID=idd[i],Dose=idzAmt,t(NCAprm),NSIM=s))
          }
        }
        if (case == 2){
          for (s1 in 1:npopStr1){
            ifdf <- smdf[smdf[,popStrNm1]==popStr1[s1],]
            if (nrow(ifdf) == 0){next}
            idd <- unique(as.character(ifdf[,idNmSim]))
            if (is.data.frame(idd)) idd <- idd[[1]]
            for (i in 1:length(idd)){
              if (!is.null(doseAmtNm)){
                doseData <- as.numeric(as.character(srdf[srdf[,popStrNm1]==popStr1[s1] & srdf[,idNmSim]==idd[i], doseAmtNm]))
                doseData <- doseData[complete.cases(doseData) & doseData>0]
                idzAmt   <- ifelse(length(doseData)==0, NA, doseData[1])
              }else{
                idzAmt <- NA
              }
              stcTI <- nca.ind.data(pkData=ifdf, ID=idd[i], dvLog = simLog, dataType=dset,
                                    idNm=idNmSim, timeNm=timeNmSim, concNm=concNmSim,
                                    adminType=adminType, TI=TI,
                                    dateColNm=dateColNm, dateFormat=dateFormat, timeFormat=timeFormat)
              stc     <- stcTI$tc
              iTI     <- stcTI$iTI
              if (nrow(stc)==0) next
              time    <- as.numeric(stc$time)
              conc    <- as.numeric(stc$conc)
              NCAprm  <- est.nca(time=time,conc=conc,backExtrp=backExtrp,negConcExcl=negConcExcl,
                                 doseType=doseType,adminType=adminType,doseAmt=idzAmt,method=method,
                                 AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,LambdaExclude=LambdaExclude,
                                 doseTime=doseTime,Tau=Tau,TI=iTI,simFile=simFile,dset=dset,onlyNCA=onlyNCA)
              simData <- rbind(simData, data.frame(ID=idd[i],STRAT1=popStr1[s1],Dose=idzAmt,t(NCAprm),NSIM=s))
            }
          }
        }
        if (case == 3){
          for (s1 in 1:npopStr1){
            for (s2 in 1:npopStr2){
              ifdf <- smdf[smdf[,popStrNm1]==popStr1[s1] & smdf[,popStrNm2]==popStr2[s2],]
              if (nrow(ifdf) == 0){next}
              idd <- unique(as.character(ifdf[,idNmSim]))
              for (i in 1:length(idd)){
                if (!is.null(doseAmtNm)){
                  doseData <- as.numeric(as.character(srdf[srdf[,popStrNm1]==popStr1[s1] & srdf[,popStrNm2]==popStr2[s2] & srdf[,idNmSim]==idd[i], doseAmtNm]))
                  doseData <- doseData[complete.cases(doseData) & doseData>0]
                  idzAmt   <- ifelse(length(doseData)==0, NA, doseData[1])
                }else{
                  idzAmt <- NA
                }
                stcTI <- nca.ind.data(pkData=ifdf, ID=idd[i], dvLog = simLog, dataType=dset,
                                      idNm=idNmSim, timeNm=timeNmSim, concNm=concNmSim,
                                      adminType=adminType, TI=TI,
                                      dateColNm=dateColNm, dateFormat=dateFormat, timeFormat=timeFormat)
                stc     <- stcTI$tc
                iTI     <- stcTI$iTI
                if (nrow(stc)==0) next
                time    <- as.numeric(stc$time)
                conc    <- as.numeric(stc$conc)
                NCAprm  <- est.nca(time=time,conc=conc,backExtrp=backExtrp,negConcExcl=negConcExcl,
                                   doseType=doseType,adminType=adminType,doseAmt=idzAmt,method=method,
                                   AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,LambdaExclude=LambdaExclude,
                                   doseTime=doseTime,Tau=Tau,TI=iTI,simFile=simFile,dset=dset,onlyNCA=onlyNCA)
                simData <- rbind(simData, data.frame(ID=idd[i],STRAT1=popStr1[s1],STRAT2=popStr2[s2],Dose=idzAmt,t(NCAprm),NSIM=s))
              }
            }
          }
        }
        if (case == 4){
          for (s1 in 1:npopStr1){
            for (s2 in 1:npopStr2){
              for (s3 in 1:npopStr3){
                ifdf <- smdf[smdf[,popStrNm1]==popStr1[s1] & smdf[,popStrNm2]==popStr2[s2] & smdf[,popStrNm3]==popStr3[s3],]
                if (nrow(ifdf) == 0){next}
                idd <- unique(as.character(ifdf[,idNmSim]))
                for (i in 1:length(idd)){
                  if (!is.null(doseAmtNm)){
                    doseData <- as.numeric(as.character(srdf[srdf[,popStrNm1]==popStr1[s1] & srdf[,popStrNm2]==popStr2[s2] & srdf[,popStrNm3]==popStr3[s3] & srdf[,idNmSim]==idd[i], doseAmtNm]))
                    doseData <- doseData[complete.cases(doseData) & doseData>0]
                    idzAmt   <- ifelse(length(doseData)==0, NA, doseData[1])
                  }else{
                    idzAmt <- NA
                  }
                  stcTI <- nca.ind.data(pkData=ifdf, ID=idd[i], dvLog = simLog, dataType=dset,
                                        idNm=idNmSim, timeNm=timeNmSim, concNm=concNmSim,
                                        adminType=adminType, TI=TI,
                                        dateColNm=dateColNm, dateFormat=dateFormat, timeFormat=timeFormat)
                  stc     <- stcTI$tc
                  iTI     <- stcTI$iTI
                  if (nrow(stc)==0) next
                  time    <- as.numeric(stc$time)
                  conc    <- as.numeric(stc$conc)
                  NCAprm  <- est.nca(time=time,conc=conc,backExtrp=backExtrp,negConcExcl=negConcExcl,doseType=doseType,adminType=adminType,doseAmt=idzAmt,method=method,AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,LambdaExclude=LambdaExclude,doseTime=doseTime,Tau=Tau,TI=iTI,simFile=simFile,dset=dset,onlyNCA=onlyNCA)
                  simData <- rbind(simData, data.frame(ID=idd[i],STRAT1=popStr1[s1],STRAT2=popStr2[s2],STRAT3=popStr3[s3],Dose=idzAmt,t(NCAprm),NSIM=s))
                }
              }
            }
          }
        }
        
        if (case == 1) names(simData)[names(simData)%in%c("ID")] <- c(idNmSim)
        if (case == 2) names(simData)[names(simData)%in%c("ID","STRAT1")] <- c(idNmSim,popStrNm1)
        if (case == 3) names(simData)[names(simData)%in%c("ID","STRAT1","STRAT2")] <- c(idNmSim,popStrNm1,popStrNm2)
        if (case == 4) names(simData)[names(simData)%in%c("ID","STRAT1","STRAT2","STRAT3")] <- c(idNmSim,popStrNm1,popStrNm2,popStrNm3)
        write.csv(simData, file=paste(simDir,"/sim_",s,".csv",sep=""), row.names=F, quote=F)
      }
    }
    
    
    # read all simulated NCA parameters to a list
    lasdf <- lapply(list.files(path = simDir, pattern="sim_[0-9]*.csv",full.names=T),function(i){read.csv(i, header=T)})
    nsim  <- length(lasdf)
    dasdf <- do.call(rbind, lapply(lasdf, as.data.frame))
    
    # Rename the ID and stratification column names to ID, STRAT1, STRAT2 and/STRAT3
    if (case==1) lasdf <- lapply(seq(lasdf), function(i){x <- data.frame(lasdf[[i]]); names(x)[match(c(idNmSim),names(x))] <- c("ID"); return(x)})
    if (case==2) lasdf <- lapply(seq(lasdf), function(i){x <- data.frame(lasdf[[i]]); names(x)[match(c(idNmSim,popStrNm1),names(x))] <- c("ID","STRAT1"); return(x)})
    if (case==3) lasdf <- lapply(seq(lasdf), function(i){x <- data.frame(lasdf[[i]]); names(x)[match(c(idNmSim,popStrNm1,popStrNm2),names(x))] <- c("ID","STRAT1","STRAT2"); return(x)})
    if (case==4) lasdf <- lapply(seq(lasdf), function(i){x <- data.frame(lasdf[[i]]); names(x)[match(c(idNmSim,popStrNm1,popStrNm2,popStrNm3),names(x))] <- c("ID","STRAT1","STRAT2","STRAT3"); return(x)})
    
    
    if (printOut){
      if (is.null(outFileNm) || outFileNm==""){
        write.table(dasdf, file=paste0(usrdir,"/ncaSimEst.tsv"), row.names=F, quote=F, sep="\t")
      }else{
        write.table(dasdf, file=paste0(usrdir,"/ncaSimEst-",outFileNm,".tsv"), row.names=F, quote=F, sep="\t")
      }
    }
    
    if (case==1) names(dasdf)[match(c(idNmSim),names(dasdf))] <- c("ID")
    if (case==2) names(dasdf)[match(c(idNmSim,popStrNm1),names(dasdf))] <- c("ID","STRAT1")
    if (case==3) names(dasdf)[match(c(idNmSim,popStrNm1,popStrNm2),names(dasdf))] <- c("ID","STRAT1","STRAT2")
    if (case==4) names(dasdf)[match(c(idNmSim,popStrNm1,popStrNm2,popStrNm3),names(dasdf))] <- c("ID","STRAT1","STRAT2","STRAT3")
    
    
    # Set plot output dimensions
    if (npr<=2){
      hth<-14; wth<-15; phth<-7; pwth<-6
    }else if (npr>2 & npr<=4){
      hth<-20; wth<-16; phth<-8; pwth<-7
    }else if (npr>4 & npr<=6){
      hth<-26; wth<-16; phth<-9; pwth<-8
    }else if (npr>6){
      hth<-26; wth<-22; phth<-9; pwth<-10
    }
    
    
    # Plot population histogram
    pop_hist_list <- list()
    if (!noPlot){
      # Population histogram (case=1)
      if (case == 1){
        smeanData   <- data.frame()
        smedianData <- data.frame()
        svarData <- data.frame()
        for (i in 1:length(lasdf)){
          tmdf        <- subset(data.frame(lasdf[[i]]), select=param)
          meanPrm     <- as.data.frame(lapply(tmdf, FUN=function(x) mean(as.numeric(x[!is.na(x)]))))
          smeanData   <- rbind(smeanData, meanPrm)
          medianPrm   <- as.data.frame(lapply(tmdf, FUN=function(x) median(as.numeric(x[!is.na(x)]))))
          smedianData <- rbind(smedianData, medianPrm)
          varPrm   <- as.data.frame(lapply(tmdf, FUN=function(x) var(as.numeric(x[!is.na(x)]))))
          svarData <- rbind(svarData, varPrm)
        }
        obsdata     <- subset(outData, select=param, ID!="")
        figlbl      <- NULL
        histpopgrob <- histpop.plot(obsdata=obsdata,simdata=smedianData,figlbl=figlbl,param=param,cunit=cunit,tunit=tunit,spread=spread)
        gdr         <- histpopgrob$gdr
        mylegend    <- histpopgrob$legend
        lheight     <- histpopgrob$lheight
        pop_hist_tmp <- pop_hist(obsdata,
                                 smeanData, 
                                 svarData,
                                 title=figlbl)
        if (printOut){
          fl <- paste0(usrdir,"/PopMean")
          if (figFormat=="tiff"){
            eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",compression=\"lzw\",res=120)")))
          }else{
            eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",res=120)")))
          }
          suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
          dev.off()
        }
        suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
        ggr <- grid.grab()
        popplot[[length(popplot)+1]] <- ggr
        pop_hist_list[[length(pop_hist_list)+1]] <- pop_hist_tmp
        
      }
      
      # Population histogram (case=2)
      if (case == 2){
        for (s1 in 1:npopStr1){
          if (nrow(dasdf[dasdf$STRAT1==popStr1[s1],]) == 0) next
          smeanData   <- data.frame()
          smedianData <- data.frame()
          svarData <- data.frame()
          
          for (i in 1:length(lasdf)){
            tmdf        <- subset(data.frame(lasdf[[i]]), select=param, STRAT1==popStr1[s1])
            meanPrm     <- as.data.frame(lapply(tmdf, FUN=function(x) mean(as.numeric(x[!is.na(x)]))))
            smeanData   <- rbind(smeanData, meanPrm)
            medianPrm   <- as.data.frame(lapply(tmdf, FUN=function(x) median(as.numeric(x[!is.na(x)]))))
            smedianData <- rbind(smedianData, medianPrm)
            varPrm   <- as.data.frame(lapply(tmdf, FUN=function(x) var(as.numeric(x[!is.na(x)]))))
            svarData <- rbind(svarData, varPrm)
          }
          obsdata     <- subset(outData, select=param, ID!="" & STRAT1==popStr1[s1])
          figlbl      <- paste0(popStrNm1,"-",popStr1[s1])
          histpopgrob <- histpop.plot(obsdata=obsdata,simdata=smedianData,figlbl=figlbl,param=param,cunit=cunit,tunit=tunit,spread=spread)
          gdr         <- histpopgrob$gdr
          mylegend    <- histpopgrob$legend
          lheight     <- histpopgrob$lheight
          pop_hist_tmp <- pop_hist(obsdata,
                                   smeanData, 
                                   svarData,
                                   title=figlbl)
          if (printOut){
            fl <- paste0(usrdir,"/PopMean_",figlbl)
            if (figFormat=="tiff"){
              eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",compression=\"lzw\",res=120)")))
            }else{
              eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",res=120)")))
            }
            suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
            dev.off()
          }
          suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
          ggr <- grid.grab()
          popplot[[length(popplot)+1]] <- ggr
          pop_hist_list[[length(pop_hist_list)+1]] <- pop_hist_tmp
        }
      }
      
      # Population histogram (case=3)
      if (case == 3){
        for (s1 in 1:npopStr1){
          for (s2 in 1:npopStr2){
            if (nrow(dasdf[dasdf$STRAT1==popStr1[s1] & dasdf$STRAT2==popStr2[s2],]) == 0) next
            smeanData   <- data.frame()
            smedianData <- data.frame()
            svarData <- data.frame()
            
            for (i in 1:length(lasdf)){
              tmdf        <- subset(data.frame(lasdf[[i]]), select=param, STRAT1==popStr1[s1] & STRAT2==popStr2[s2])
              meanPrm     <- as.data.frame(lapply(tmdf, FUN=function(x) mean(as.numeric(x[!is.na(x)]))))
              smeanData   <- rbind(smeanData, meanPrm)
              medianPrm   <- as.data.frame(lapply(tmdf, FUN=function(x) median(as.numeric(x[!is.na(x)]))))
              smedianData <- rbind(smedianData, medianPrm)
              varPrm   <- as.data.frame(lapply(tmdf, FUN=function(x) var(as.numeric(x[!is.na(x)]))))
              svarData <- rbind(svarData, varPrm)
            }
            obsdata     <- subset(outData, select=param, ID!="" & STRAT1==popStr1[s1] & STRAT2==popStr2[s2])
            figlbl      <- paste0(popStrNm1,"-",popStr1[s1],"_",popStrNm2,"-",popStr2[s2])
            histpopgrob <- histpop.plot(obsdata=obsdata,simdata=smedianData,figlbl=figlbl,param=param,cunit=cunit,tunit=tunit,spread=spread)
            gdr         <- histpopgrob$gdr
            mylegend    <- histpopgrob$legend
            lheight     <- histpopgrob$lheight
            pop_hist_tmp <- pop_hist(obsdata,
                                     smeanData, 
                                     svarData,
                                     title=figlbl)
            if (printOut){
              fl <- paste0(usrdir,"/PopMean_",figlbl)
              if (figFormat=="tiff"){
                eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",compression=\"lzw\",res=120)")))
              }else{
                eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",res=120)")))
              }
              suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
              dev.off()
            }
            suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
            ggr <- grid.grab()
            popplot[[length(popplot)+1]] <- ggr
            pop_hist_list[[length(pop_hist_list)+1]] <- pop_hist_tmp
          }
        }
      }
      
      # Population histogram (case=4)
      if (case == 4){
        for (s1 in 1:npopStr1){
          for (s2 in 1:npopStr2){
            for (s3 in 1:npopStr3){
              if (nrow(dasdf[dasdf$STRAT1==popStr1[s1] & dasdf$STRAT2==popStr2[s2] & dasdf$STRAT3==popStr3[s3],]) == 0) next
              smeanData   <- data.frame()
              smedianData <- data.frame()
              svarData <- data.frame()
              
              for (i in 1:length(lasdf)){
                tmdf        <- subset(data.frame(lasdf[[i]]), select=param, STRAT1==popStr1[s1] & STRAT2==popStr2[s2] & STRAT3==popStr3[s3])
                meanPrm     <- as.data.frame(lapply(tmdf, FUN=function(x) mean(as.numeric(x[!is.na(x)]))))
                smeanData   <- rbind(smeanData, meanPrm)
                medianPrm   <- as.data.frame(lapply(tmdf, FUN=function(x) median(as.numeric(x[!is.na(x)]))))
                smedianData <- rbind(smedianData, medianPrm)
                varPrm   <- as.data.frame(lapply(tmdf, FUN=function(x) var(as.numeric(x[!is.na(x)]))))
                svarData <- rbind(svarData, varPrm)
              }
              obsdata     <- subset(outData, select=param, ID!="" & STRAT1==popStr1[s1] & STRAT2==popStr2[s2] & STRAT3==popStr3[s3])
              figlbl      <- paste0(popStrNm1,"-",popStr1[s1],"_",popStrNm2,"-",popStr2[s2],"_",popStrNm3,"-",popStr3[s3])
              histpopgrob <- histpop.plot(obsdata=obsdata,simdata=smedianData,figlbl=figlbl,param=param,cunit=cunit,tunit=tunit,spread=spread)
              gdr         <- histpopgrob$gdr
              mylegend    <- histpopgrob$legend
              lheight     <- histpopgrob$lheight
              pop_hist_tmp <- pop_hist(obsdata,
                                       smeanData, 
                                       svarData,
                                       title=figlbl)
              if (printOut){
                fl <- paste0(usrdir,"/PopMean_",figlbl)
                if (figFormat=="tiff"){
                  eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",compression=\"lzw\",res=120)")))
                }else{
                  eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",res=120)")))
                }
                suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
                dev.off()
              }
              suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
              ggr <- grid.grab()
              popplot[[length(popplot)+1]] <- ggr
              pop_hist_list[[length(pop_hist_list)+1]] <- pop_hist_tmp
            }
          }
        }
      }
    }
    
    
    # Perform ppc; sensitive to onlyNCA
    if (!onlyNCA){
      devcol  <- paste0("d",param)
      npdecol <- paste0("npde",param)
      
      # ggplot options for the forest plot
      ggOpt_forest <- list(xlab("\nNPDE"),ylab(""),
                           labs(title = "Forest plot of NPDE\nErrorbar = 95% confidence interval\n"),
                           scale_color_manual(name="",values=c("mean"="red","SD"="darkgreen")),
                           theme(axis.text.x = element_text(vjust=1,hjust=1,size=10),
                                 axis.text.y = element_text(hjust=0,size=10),
                                 strip.text.x = element_text(size=10),
                                 legend.text = element_text(size=12),
                                 title = element_text(size=14,face="bold"),
                                 legend.position = "bottom", legend.direction = "horizontal",
                                 legend.background = element_rect(),
                                 legend.key.height = unit(1,"cm")),
                           facet_wrap(~type, scales="fixed", ncol=2))
      
      OTL   <- data.frame(No_of_outliers=numeric(0),ID_metric=character(0))
      npde  <- data.frame()
      fpval <- data.frame(type=character(0),mean=numeric(0),mcil=numeric(0),mciu=numeric(0),sdu=numeric(0),sducil=numeric(0),sduciu=numeric(0),str=character(0))
    }
    
    if (case == 1){
      tdasdf <- dasdf
      id     <- unique(tdasdf$ID)
      pde    <- data.frame()
      metric <- ""
      nout   <- 0
      for (i in 1:length(id)){
        obsdata <- subset(outData, ID==id[i])
        simdata <- subset(tdasdf, ID==id[i])
        figlbl  <- NULL
        pdeout  <- nca.pde.deviation.outlier(obsdata=obsdata,simdata=simdata,idNm="ID",id=id[i],spread=spread,figlbl=figlbl,calcparam=alwprm,diagparam=param,cunit=cunit,tunit=tunit,noPlot=noPlot,onlyNCA=onlyNCA)
        pde     <- rbind(pde, cbind(data.frame(ID=id[i]), pdeout$pde))
        outData[outData$ID==id[i],] <- pdeout$obsdata
        if (!onlyNCA && pdeout$metric != ""){
          nout   <- nout + 1
          metric <- paste(metric,pdeout$metric,sep=", ")
          if (!noPlot){
            gdr      <- pdeout$grob
            mylegend <- pdeout$legend
            lheight  <- pdeout$lheight
            if (printOut){
              fl <- paste0(usrdir,"/Outlier_ID-",id[i])
              if (figFormat=="tiff"){
                eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",compression=\"lzw\",res=120)")))
              }else{
                eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",res=120)")))
              }
              suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
              dev.off()
            }
            suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
            ggr <- grid.grab()
            outlierplot[[length(outlierplot)+1]] <- ggr
          }
        }
      }
      if (!onlyNCA){
        if (metric != "") metric <- gsub("^, ", "", metric)
        OTL  <- rbind(OTL, data.frame(No_of_outlier=nout,ID_metric=metric))
        npde <- rbind(npde,pde)
        
        npde   <- nca.npde(pdedata=npde,pdecol=alwprm)
        npdeNm <- paste0("npde",alwprm)
        for (r in 1:nrow(outData)){
          if (nrow(npde[npde$ID==outData$ID[r],])!=1){
            outData[r,npdeNm] <- NaN
          }else{
            outData[r,npdeNm] <- npde[npde$ID==outData$ID[r],npdeNm]
          }
        }
        plotdata <- outData
        if (nrow(plotdata) == 0) next
        figlbl <- "All data"
        # Deviation plot
        if (!noPlot){
          ggdev <- nca.deviation.plot(plotdata=plotdata,xvar="ID",devcol=devcol,figlbl=figlbl,spread=spread,cunit=cunit,tunit=tunit)
          if (!is.null(ggdev)){
            suppressMessages(suppressWarnings(print(ggdev)))
            if (printOut) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/Deviation_",figlbl,".",figFormat),height=hth,width=wth,units="cm",dpi=120)))
            devplot[[length(devplot)+1]] <- ggdev
          }
        }
        # NPDE plot
        npdeout <- nca.npde.plot(plotdata=plotdata,xvar="ID",npdecol=npdecol,figlbl=figlbl,cunit=cunit,tunit=tunit)
        if (is.null(npdeout$forestdata)) next
        forestdata <- npdeout$forestdata
        forestdata$str <- figlbl
        fpval <- rbind(fpval, forestdata)
        if (!noPlot){
          npdeplot[[length(npdeplot)+1]] <- npdeout$ggnpde
          suppressMessages(suppressWarnings(print(npdeout$ggnpde)))
          if (printOut) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/NPDE_",figlbl,".",figFormat),height=hth,width=wth,units="cm",dpi=120)))
          histnpdeplot[[length(histnpdeplot)+1]] <- npdeout$gghnpde
          suppressMessages(suppressWarnings(print(npdeout$gghnpde)))
          if (printOut) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/HistNPDE_",figlbl,".",figFormat),height=hth,width=wth,units="cm",dpi=120)))
        }
        if (!noPlot){
          # Forest plot for NPDE
          fpval$FCT <- paste0("mean=",signif(fpval$mean,2),"+/-CI=",signif(fpval$mcil,2),",",signif(fpval$mciu,2),", SD=",signif(fpval$sd,2),"+/-CI=",signif(fpval$sdcil,2),",",signif(fpval$sdciu,2))
          xlim1 <- floor(min(fpval$mcil)); xlim2 <- ceiling(fpval$sdciu)
          ggplt <- ggplot(fpval) + ggOpt_forest +
            geom_point(aes(mean,str,color="mean"), show.legend=T, size=2) +
            geom_errorbarh(aes(x=mean,y=str,xmin=mcil,xmax=mciu),size=0.4, color="red",height=0.1) +
            geom_point(aes(sd,str,color="SD"), size=2) +
            geom_errorbarh(aes(x=sd,y=str,xmin=sdcil,xmax=sdciu), size=0.4, color="darkgreen", height=0.1) +
            geom_text(aes(label=out.digits(mean,dig=2),x=mean,y=str,color="mean",vjust=-1),size=4,show.legend=F) +
            geom_text(aes(label=out.digits(mcil,dig=2),x=mcil,y=str,color="mean",vjust=-2.5),size=4,show.legend=F) +
            geom_text(aes(label=out.digits(mciu,dig=2),x=mciu,y=str,color="mean",vjust=-2.5),size=4,show.legend=F) +
            geom_text(aes(label=out.digits(sd,dig=2),x=sd,y=str,color="SD",vjust=1.5),size=4,show.legend=F) +
            geom_text(aes(label=out.digits(sdcil,dig=2),x=sdcil,y=str,color="SD",vjust=2.5),size=4,show.legend=F) +
            geom_text(aes(label=out.digits(sdciu,dig=2),x=sdciu,y=str,color="SD",vjust=2.5),size=4,show.legend=F)
          suppressMessages(suppressWarnings(print(ggplt)))
          forestplot[[length(forestplot)+1]] <- ggplt
          if (printOut) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/ForestNPDE.",figFormat),height=hth,width=wth,units="cm",dpi=120)))
        }
      }
    }else if (case==2){
      for (s1 in 1:npopStr1){
        if (nrow(dasdf[dasdf$STRAT1==popStr1[s1],]) == 0) next
        tdasdf <- subset(dasdf, STRAT1==popStr1[s1])
        id     <- unique(tdasdf$ID)
        pde    <- data.frame()
        metric <- ""
        nout   <- 0
        figlbl <- paste0(popStrNm1,"-",popStr1[s1])
        for (i in 1:length(id)){
          obsdata <- subset(outData, ID==id[i] & STRAT1==popStr1[s1])
          simdata <- subset(tdasdf, ID==id[i])
          pdeout  <- nca.pde.deviation.outlier(obsdata=obsdata,simdata=simdata,idNm="ID",id=id[i],spread=spread,figlbl=figlbl,calcparam=alwprm,diagparam=param,cunit=cunit,tunit=tunit,noPlot=noPlot,onlyNCA=onlyNCA)
          pde     <- rbind(pde, cbind(data.frame(ID=id[i],STRAT1=popStr1[s1]), pdeout$pde))
          outData[(outData$ID==id[i] & outData$STRAT1==popStr1[s1]),] <- pdeout$obsdata
          if (!onlyNCA && pdeout$metric != ""){
            nout     <- nout + 1
            metric   <- paste(metric,pdeout$metric,sep=", ")
            if (!noPlot){
              gdr      <- pdeout$grob
              mylegend <- pdeout$legend
              lheight  <- pdeout$lheight
              if (printOut){
                fl <- paste0(usrdir,"/Outlier_ID-",id[i],"_",figlbl)
                if (figFormat=="tiff"){
                  eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",compression=\"lzw\",res=120)")))
                }else{
                  eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",res=120)")))
                }
                
                suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
                dev.off()
              }
              suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
              ggr <- grid.grab()
              outlierplot[[length(outlierplot)+1]] <- ggr
            }
          }
        }
        if (!onlyNCA){
          if (metric != "") metric <- gsub("^, ", "", metric)
          OTL  <- rbind(OTL, data.frame(No_of_outlier=nout,ID_metric=metric))
          npde <- rbind(npde,pde)
        }
      }
      if (!onlyNCA){
        npde   <- nca.npde(pdedata=npde,pdecol=alwprm)
        npdeNm <- paste0("npde",alwprm)
        for (r in 1:nrow(outData)){
          if (nrow(npde[(npde$ID==outData$ID[r] & npde$STRAT1==outData$STRAT1[r]),])!=1){
            outData[r,npdeNm] <- NaN
          }else{
            outData[r,paste0("npde",alwprm)] <- npde[(npde$ID==outData$ID[r] & npde$STRAT1==outData$STRAT1[r]),paste0("npde",alwprm)]
          }
        }
        for (s1 in 1:npopStr1){
          plotdata <- subset(outData, STRAT1==popStr1[s1])
          if (nrow(plotdata) == 0) next
          figlbl <- paste0(popStrNm1,"-",popStr1[s1])
          # Deviation plot
          if (!noPlot){
            ggdev <- nca.deviation.plot(plotdata=plotdata,xvar="ID",devcol=devcol,figlbl=figlbl,spread=spread,cunit=cunit,tunit=tunit)
            if (!is.null(ggdev)){
              suppressMessages(suppressWarnings(print(ggdev)))
              if (printOut) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/Deviation_",figlbl,".",figFormat),height=hth,width=wth,units="cm",dpi=120)))
              devplot[[length(devplot)+1]] <- ggdev
            }
          }
          # NPDE plot
          npdeout <- nca.npde.plot(plotdata=plotdata,xvar="ID",npdecol=npdecol,figlbl=figlbl,cunit=cunit,tunit=tunit)
          if (is.null(npdeout$forestdata)) next
          forestdata <- npdeout$forestdata
          forestdata$str <- figlbl
          fpval <- rbind(fpval, forestdata)
          if (!noPlot){
            npdeplot[[length(npdeplot)+1]] <- npdeout$ggnpde
            suppressMessages(suppressWarnings(print(npdeout$ggnpde)))
            if (printOut) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/NPDE_",figlbl,".",figFormat),height=hth,width=wth,units="cm",dpi=120)))
            histnpdeplot[[length(histnpdeplot)+1]] <- npdeout$gghnpde
            suppressMessages(suppressWarnings(print(npdeout$gghnpde)))
            if (printOut) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/HistNPDE_",figlbl,".",figFormat),height=hth,width=wth,units="cm",dpi=120)))
          }
        }
        if (!noPlot){
          # Forest plot for NPDE
          fpval$FCT <- paste0("mean=",signif(fpval$mean,2),"+/-CI=",signif(fpval$mcil,2),",",signif(fpval$mciu,2),", SD=",signif(fpval$sd,2),"+/-CI=",signif(fpval$sdcil,2),",",signif(fpval$sdciu,2))
          xlim1 <- floor(min(fpval$mcil)); xlim2 <- ceiling(fpval$sdciu)
          ggplt <- ggplot(fpval) + ggOpt_forest +
            geom_point(aes(mean,str,color="mean"), show.legend=T, size=2) +
            geom_errorbarh(aes(x=mean,y=str,xmin=mcil,xmax=mciu),size=0.4, color="red",height=0.1) +
            geom_point(aes(sd,str,color="SD"), size=2) +
            geom_errorbarh(aes(x=sd,y=str,xmin=sdcil,xmax=sdciu), size=0.4, color="darkgreen", height=0.1) +
            geom_text(aes(label=out.digits(mean,dig=2),x=mean,y=str,color="mean",vjust=-1),size=4,show.legend=F) +
            geom_text(aes(label=out.digits(mcil,dig=2),x=mcil,y=str,color="mean",vjust=-2.5),size=4,show.legend=F) +
            geom_text(aes(label=out.digits(mciu,dig=2),x=mciu,y=str,color="mean",vjust=-2.5),size=4,show.legend=F) +
            geom_text(aes(label=out.digits(sd,dig=2),x=sd,y=str,color="SD",vjust=1.5),size=4,show.legend=F) +
            geom_text(aes(label=out.digits(sdcil,dig=2),x=sdcil,y=str,color="SD",vjust=2.5),size=4,show.legend=F) +
            geom_text(aes(label=out.digits(sdciu,dig=2),x=sdciu,y=str,color="SD",vjust=2.5),size=4,show.legend=F)
          suppressMessages(suppressWarnings(print(ggplt)))
          forestplot[[length(forestplot)+1]] <- ggplt
          if (printOut) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/ForestNPDE.",figFormat),height=hth,width=wth,units="cm",dpi=120)))
        }
      }
    }else if (case == 3){
      for (s1 in 1:npopStr1){
        for (s2 in 1:npopStr2){
          if (nrow(dasdf[dasdf$STRAT1==popStr1[s1] & dasdf$STRAT2==popStr2[s2],]) == 0) next
          tdasdf <- subset(dasdf, STRAT1==popStr1[s1] & STRAT2==popStr2[s2])
          id     <- unique(tdasdf$ID)
          pde    <- data.frame()
          metric <- ""
          nout   <- 0
          for (i in 1:length(id)){
            obsdata <- subset(outData, ID==id[i] & STRAT1==popStr1[s1] & STRAT2==popStr2[s2])
            simdata <- subset(tdasdf, ID==id[i])
            figlbl  <- paste0(popStrNm1,"-",popStr1[s1],"_",popStrNm2,"-",popStr2[s2])
            pdeout  <- nca.pde.deviation.outlier(obsdata=obsdata,simdata=simdata,idNm="ID",id=id[i],spread=spread,figlbl=figlbl,calcparam=alwprm,diagparam=param,cunit=cunit,tunit=tunit,noPlot=noPlot,onlyNCA=onlyNCA)
            pde     <- rbind(pde, cbind(data.frame(ID=id[i],STRAT1=popStr1[s1],STRAT2=popStr2[s2]), pdeout$pde))
            outData[(outData$ID==id[i] & outData$STRAT1==popStr1[s1] & outData$STRAT2==popStr2[s2]),] <- pdeout$obsdata
            if (!onlyNCA && pdeout$metric != ""){
              nout     <- nout + 1
              metric   <- paste(metric,pdeout$metric,sep=", ")
              if (!noPlot){
                gdr      <- pdeout$grob
                mylegend <- pdeout$legend
                lheight  <- pdeout$lheight
                if (printOut){
                  fl <- paste0(usrdir,"/Outlier_ID-",id[i],"_",figlbl)
                  if (figFormat=="tiff"){
                    eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",compression=\"lzw\",res=120)")))
                  }else{
                    eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",res=120)")))
                  }
                  suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
                  dev.off()
                }
                suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
                ggr <- grid.grab()
                outlierplot[[length(outlierplot)+1]] <- ggr
              }
            }
          }
          if (!onlyNCA){
            if (metric != "") metric <- gsub("^, ", "", metric)
            OTL  <- rbind(OTL, data.frame(No_of_outlier=nout,ID_metric=metric))
            npde <- rbind(npde,pde)
          }
        }
      }
      if (!onlyNCA){
        npde   <- nca.npde(pdedata=npde,pdecol=alwprm)
        npdeNm <- paste0("npde",alwprm)
        for (r in 1:nrow(outData)){
          if (nrow(npde[(npde$ID==outData$ID[r] & npde$STRAT1==outData$STRAT1[r] & npde$STRAT2==outData$STRAT2[r]),])!=1){
            outData[r,npdeNm] <- NaN
          }else{
            outData[r,paste0("npde",alwprm)] <- npde[(npde$ID==outData$ID[r] & npde$STRAT1==outData$STRAT1[r] & npde$STRAT2==outData$STRAT2[r]),paste0("npde",alwprm)]
          }
        }
        for (s1 in 1:npopStr1){  
          for (s2 in 1:npopStr2){
            plotdata <- subset(outData, STRAT1==popStr1[s1] & STRAT2==popStr2[s2])
            if (nrow(plotdata) == 0) next
            figlbl <- paste0(popStrNm1,"-",popStr1[s1],"_",popStrNm2,"-",popStr2[s2])
            # Deviation plot
            if (!noPlot){
              ggdev <- nca.deviation.plot(plotdata=plotdata,xvar="ID",devcol=devcol,figlbl=figlbl,spread=spread,cunit=cunit,tunit=tunit)
              if (!is.null(ggdev)){
                suppressMessages(suppressWarnings(print(ggdev)))
                if (printOut) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/Deviation_",figlbl,".",figFormat),height=hth,width=wth,units="cm",dpi=120)))
                devplot[[length(devplot)+1]] <- ggdev
              }
            }
            # NPDE plot
            npdeout <- nca.npde.plot(plotdata=plotdata,xvar="ID",npdecol=npdecol,figlbl=figlbl,cunit=cunit,tunit=tunit)
            if (is.null(npdeout$forestdata)) next
            forestdata <- npdeout$forestdata
            forestdata$str <- figlbl
            fpval <- rbind(fpval, forestdata)
            if (!noPlot){
              npdeplot[[length(npdeplot)+1]] <- npdeout$ggnpde
              suppressMessages(suppressWarnings(print(npdeout$ggnpde)))
              if (printOut) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/NPDE_",figlbl,".",figFormat),height=hth,width=wth,units="cm",dpi=120)))
              histnpdeplot[[length(histnpdeplot)+1]] <- npdeout$gghnpde
              suppressMessages(suppressWarnings(print(npdeout$gghnpde)))
              if (printOut) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/HistNPDE_",figlbl,".",figFormat),height=hth,width=wth,units="cm",dpi=120)))
            }
          }
        }
        if (!noPlot){
          # Forest plot for NPDE
          fpval$FCT <- paste0("mean=",signif(fpval$mean,2),"+/-CI=",signif(fpval$mcil,2),",",signif(fpval$mciu,2),", SD=",signif(fpval$sd,2),"+/-CI=",signif(fpval$sdcil,2),",",signif(fpval$sdciu,2))
          xlim1 <- floor(min(fpval$mcil)); xlim2 <- ceiling(fpval$sdciu)
          ggplt <- ggplot(fpval) + ggOpt_forest +
            geom_point(aes(mean,str,color="mean"), show.legend=T, size=2) +
            geom_errorbarh(aes(x=mean,y=str,xmin=mcil,xmax=mciu),size=0.4, color="red",height=0.1) +
            geom_point(aes(sd,str,color="SD"), size=2) +
            geom_errorbarh(aes(x=sd,y=str,xmin=sdcil,xmax=sdciu), size=0.4, color="darkgreen", height=0.1) +
            geom_text(aes(label=out.digits(mean,dig=2),x=mean,y=str,color="mean",vjust=-1),size=4,show.legend=F) +
            geom_text(aes(label=out.digits(mcil,dig=2),x=mcil,y=str,color="mean",vjust=-2.5),size=4,show.legend=F) +
            geom_text(aes(label=out.digits(mciu,dig=2),x=mciu,y=str,color="mean",vjust=-2.5),size=4,show.legend=F) +
            geom_text(aes(label=out.digits(sd,dig=2),x=sd,y=str,color="SD",vjust=1.5),size=4,show.legend=F) +
            geom_text(aes(label=out.digits(sdcil,dig=2),x=sdcil,y=str,color="SD",vjust=2.5),size=4,show.legend=F) +
            geom_text(aes(label=out.digits(sdciu,dig=2),x=sdciu,y=str,color="SD",vjust=2.5),size=4,show.legend=F)
          suppressMessages(suppressWarnings(print(ggplt)))
          forestplot[[length(forestplot)+1]] <- ggplt
          if (printOut) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/ForestNPDE.",figFormat),height=hth,width=wth,units="cm",dpi=120)))
        }
      }
    }else if (case == 4){
      for (s1 in 1:npopStr1){
        for (s2 in 1:npopStr2){
          for (s3 in 1:npopStr3){
            if (nrow(dasdf[dasdf$STRAT1==popStr1[s1] & dasdf$STRAT2==popStr2[s2] & dasdf$STRAT3==popStr3[s3],]) == 0) next
            tdasdf <- subset(dasdf, STRAT1==popStr1[s1] & STRAT2==popStr2[s2] & STRAT3==popStr3[s3])
            id     <- unique(tdasdf$ID)
            pde    <- data.frame()
            metric <- ""
            nout   <- 0
            for (i in 1:length(id)){
              obsdata <- subset(outData, ID==id[i] & STRAT1==popStr1[s1] & STRAT2==popStr2[s2] & STRAT3==popStr3[s3])
              simdata <- subset(tdasdf, ID==id[i])
              if (nrow(obsdata)==0 | nrow(simdata)==0) next
              figlbl  <- paste0(popStrNm1,"-",popStr1[s1],"_",popStrNm2,"-",popStr2[s2],"_",popStrNm3,"-",popStr3[s3])
              pdeout  <- nca.pde.deviation.outlier(obsdata=obsdata,simdata=simdata,idNm="ID",id=id[i],spread=spread,figlbl=figlbl,calcparam=alwprm,diagparam=param,cunit=cunit,tunit=tunit,noPlot=noPlot,onlyNCA=onlyNCA)
              pde     <- rbind(pde, cbind(data.frame(ID=id[i],STRAT1=popStr1[s1],STRAT2=popStr2[s2],STRAT3=popStr3[s3]), pdeout$pde))
              outData[(outData$ID==id[i] & outData$STRAT1==popStr1[s1] & outData$STRAT2==popStr2[s2] & outData$STRAT3==popStr3[s3]),] <- pdeout$obsdata
              if (!onlyNCA && pdeout$metric != ""){
                nout     <- nout + 1
                metric   <- paste(metric,pdeout$metric,sep=", ")
                if (!noPlot){
                  gdr      <- pdeout$grob
                  mylegend <- pdeout$legend
                  lheight  <- pdeout$lheight
                  if (printOut){
                    fl <- paste0(usrdir,"/Outlier_ID-",id[i],"_",figlbl)
                    if (figFormat=="tiff"){
                      eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",compression=\"lzw\",res=120)")))
                    }else{
                      eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",res=120)")))
                    }
                    suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
                    dev.off()
                  }
                  suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
                  ggr <- grid.grab()
                  outlierplot[[length(outlierplot)+1]] <- ggr
                }
              }
            }
            if (!onlyNCA){
              if (metric != "") metric <- gsub("^, ", "", metric)
              OTL  <- rbind(OTL, data.frame(No_of_outlier=nout,ID_metric=metric))
              npde <- rbind(npde,pde)
            }
          }
        }
      }
      if (!onlyNCA){
        npde   <- nca.npde(pdedata=npde,pdecol=alwprm)
        npdeNm <- paste0("npde",alwprm)
        for (r in 1:nrow(outData)){
          if (nrow(npde[(npde$ID==outData$ID[r] & npde$STRAT1==outData$STRAT1[r] & npde$STRAT2==outData$STRAT2[r] & npde$STRAT3==outData$STRAT3[r]),])!=1){
            outData[r,npdeNm] <- NaN
          }else{
            outData[r,paste0("npde",alwprm)] <- npde[(npde$ID==outData$ID[r] & npde$STRAT1==outData$STRAT1[r] & npde$STRAT2==outData$STRAT2[r] & npde$STRAT3==outData$STRAT3[r]),paste0("npde",alwprm)]
          }
        }
        for (s1 in 1:npopStr1){
          for (s2 in 1:npopStr2){
            for (s3 in 1:npopStr3){
              plotdata <- subset(outData, STRAT1==popStr1[s1] & STRAT2==popStr2[s2] & STRAT3==popStr3[s3])
              if (nrow(plotdata) == 0) next
              figlbl  <- paste0(popStrNm1,"-",popStr1[s1],"_",popStrNm2,"-",popStr2[s2],"_",popStrNm3,"-",popStr3[s3])
              # Deviation plot
              if (!noPlot){
                ggdev <- nca.deviation.plot(plotdata=plotdata,xvar="ID",devcol=devcol,figlbl=figlbl,spread=spread,cunit=cunit,tunit=tunit)
                if (!is.null(ggdev)){
                  suppressMessages(suppressWarnings(print(ggdev)))
                  if (printOut) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/Deviation_",figlbl,".",figFormat),height=hth,width=wth,units="cm",dpi=120)))
                  devplot[[length(devplot)+1]] <- ggdev
                }
              }
              # NPDE plot
              npdeout <- nca.npde.plot(plotdata=plotdata,xvar="ID",npdecol=npdecol,figlbl=figlbl,cunit=cunit,tunit=tunit)
              if (is.null(npdeout$forestdata)) next
              forestdata <- npdeout$forestdata
              forestdata$str <- figlbl
              fpval <- rbind(fpval, forestdata)
              if (!noPlot){
                npdeplot[[length(npdeplot)+1]] <- npdeout$ggnpde
                suppressMessages(suppressWarnings(print(npdeout$ggnpde)))
                if (printOut) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/NPDE_",figlbl,".",figFormat),height=hth,width=wth,units="cm",dpi=120)))
                histnpdeplot[[length(histnpdeplot)+1]] <- npdeout$gghnpde
                suppressMessages(suppressWarnings(print(npdeout$gghnpde)))
                if (printOut) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/HistNPDE_",figlbl,".",figFormat),height=hth,width=wth,units="cm",dpi=120)))
              }
            }
          }
        }
        if (!noPlot){
          # Forest plot for NPDE
          fpval$FCT <- paste0("mean=",signif(fpval$mean,2),"+/-CI=",signif(fpval$mcil,2),",",signif(fpval$mciu,2),", SD=",signif(fpval$sd,2),"+/-CI=",signif(fpval$sdcil,2),",",signif(fpval$sdciu,2))
          xlim1 <- floor(min(fpval$mcil)); xlim2 <- ceiling(fpval$sdciu)
          ggplt <- ggplot(fpval) + ggOpt_forest +
            geom_point(aes(mean,str,color="mean"), show.legend=T, size=2) +
            geom_errorbarh(aes(x=mean,y=str,xmin=mcil,xmax=mciu),size=0.4, color="red",height=0.1) +
            geom_point(aes(sd,str,color="SD"), size=2) +
            geom_errorbarh(aes(x=sd,y=str,xmin=sdcil,xmax=sdciu), size=0.4, color="darkgreen", height=0.1) +
            geom_text(aes(label=out.digits(mean,dig=2),x=mean,y=str,color="mean",vjust=-1),size=4,show.legend=F) +
            geom_text(aes(label=out.digits(mcil,dig=2),x=mcil,y=str,color="mean",vjust=-2.5),size=4,show.legend=F) +
            geom_text(aes(label=out.digits(mciu,dig=2),x=mciu,y=str,color="mean",vjust=-2.5),size=4,show.legend=F) +
            geom_text(aes(label=out.digits(sd,dig=2),x=sd,y=str,color="SD",vjust=1.5),size=4,show.legend=F) +
            geom_text(aes(label=out.digits(sdcil,dig=2),x=sdcil,y=str,color="SD",vjust=2.5),size=4,show.legend=F) +
            geom_text(aes(label=out.digits(sdciu,dig=2),x=sdciu,y=str,color="SD",vjust=2.5),size=4,show.legend=F)
          suppressMessages(suppressWarnings(print(ggplt)))
          forestplot[[length(forestplot)+1]] <- ggplt
          if (printOut) suppressMessages(suppressWarnings(ggsave(filename=paste0(usrdir,"/ForestNPDE.",figFormat),height=hth,width=wth,units="cm",dpi=120)))
        }
      }
    }
    
    
    # Statistical analysis for each patient group
    statDataSim <- data.frame()
    statPrm     <- c("N","Nunique","Min","Max","Mean","Median","SD","SE","CVp","CI95l","CI95u","geoMean","gCVp")
    ncaPrm <- c("Tmax","Cmax","AUClast","AUClower_upper","AUCINF_obs","AUC_pExtrap_obs","AUCINF_pred","AUC_pExtrap_pred","AUMClast","AUMCINF_obs","AUMC_pExtrap_obs","AUMCINF_pred","AUMC_pExtrap_pred","HL_Lambda_z","Rsq","Rsq_adjusted","No_points_Lambda_z")
    
    if (case == 1){
      tmpDF       <- dasdf[,ncaPrm]
      statDataSim <- sapply(tmpDF, function(x){x<-as.numeric(as.character(x)); if(length(x[complete.cases(x)])<2){rep(NA,13)}else{out.digits(calc.stat(x[complete.cases(x)]),dig=4)}})
      statDataSim <- data.frame(statDataSim)
      if (nrow(statDataSim) == 0) next
      statDataSim <- cbind(data.frame(Stat=statPrm), statDataSim)
    }
    if (case == 2){
      for (s1 in 1:npopStr1){
        tmpDF    <- dasdf[dasdf$STRAT1==popStr1[s1],ncaPrm]
        tmpStat  <- sapply(tmpDF, function(x){x<-as.numeric(as.character(x)); if(length(x[complete.cases(x)])<2){rep(NA,13)}else{out.digits(calc.stat(x[complete.cases(x)]),dig=4)}})
        tmpStat  <- data.frame(tmpStat)
        if (nrow(tmpStat) == 0) next
        tmpStat     <- cbind(data.frame(STRAT1=popStr1[s1],Stat=statPrm), tmpStat)
        statDataSim <- rbind(statDataSim, tmpStat)
      }
      names(statDataSim)[names(statDataSim)%in%"STRAT1"] <- popStrNm1
    }
    if (case == 3){
      for (s1 in 1:npopStr1){
        for (s2 in 1:npopStr2){
          tmpDF    <- dasdf[dasdf$STRAT1==popStr1[s1] & dasdf$STRAT2==popStr2[s2],ncaPrm]
          tmpStat  <- sapply(tmpDF, function(x){x<-as.numeric(as.character(x)); if(length(x[complete.cases(x)])<2){rep(NA,13)}else{out.digits(calc.stat(x[complete.cases(x)]),dig=4)}})
          tmpStat  <- data.frame(tmpStat)
          if (nrow(tmpStat) == 0) next
          tmpStat     <- cbind(data.frame(STRAT1=popStr1[s1],STRAT2=popStr2[s2],Stat=statPrm), tmpStat)
          statDataSim <- rbind(statDataSim, tmpStat)
        }
      }
      names(statDataSim)[names(statDataSim)%in%c("STRAT1","STRAT2")] <- c(popStrNm1,popStrNm2)
    }
    if (case == 4){
      for (s1 in 1:npopStr1){
        for (s2 in 1:npopStr2){
          for (s3 in 1:npopStr3){
            tmpDF    <- dasdf[dasdf$STRAT1==popStr1[s1] & dasdf$STRAT2==popStr2[s2] & dasdf$STRAT3==popStr3[s3],ncaPrm]
            tmpStat  <- sapply(tmpDF, function(x){x<-as.numeric(as.character(x)); if(length(x[complete.cases(x)])<2){rep(NA,13)}else{out.digits(calc.stat(x[complete.cases(x)]),dig=4)}})
            tmpStat  <- data.frame(tmpStat)
            if (nrow(tmpStat) == 0) next
            tmpStat     <- cbind(data.frame(STRAT1=popStr1[s1],STRAT2=popStr2[s2],STRAT3=popStr3[s3],Stat=statPrm), tmpStat)
            statDataSim <- rbind(statDataSim, tmpStat)
          }
        }
      }
      names(statDataSim)[names(statDataSim)%in%c("STRAT1","STRAT2","STRAT3")] <- c(popStrNm1,popStrNm2,popStrNm3)
    }
    
    
    if (printOut){
      if (is.null(outFileNm) || outFileNm==""){
        write.table(statDataSim, file=paste0(usrdir,"/SimStat.tsv"), sep="\t", col.names=T, row.names=F, quote=F)
      }else{
        write.table(statDataSim, file=paste0(usrdir,"/SimStat-",outFileNm,".tsv"), sep="\t", col.names=T, row.names=F, quote=F)
      }
    }
    
    # Print output table
    # Raname ID and stratifier columns and format output table sigfig
    if (case == 1){
      names(outData)[names(outData)%in%c("ID")] <- c(idNmObs)
      outData[,c(2:ncol(outData))] <- as.data.frame(lapply(outData[,c(2:ncol(outData))], FUN=function(x) out.digits(x,dig=4)))
    }
    if (case == 2){
      names(outData)[names(outData)%in%c("ID","STRAT1")] <- c(idNmObs,popStrNm1)
      outData[,c(3:ncol(outData))] <- as.data.frame(lapply(outData[,c(3:ncol(outData))], FUN=function(x) out.digits(x,dig=4)))
    }
    if (case == 3){
      names(outData)[names(outData)%in%c("ID","STRAT1","STRAT2")] <- c(idNmObs,popStrNm1,popStrNm2)
      outData[,c(4:ncol(outData))] <- as.data.frame(lapply(outData[,c(4:ncol(outData))], FUN=function(x) out.digits(x,dig=4)))
    }
    if (case == 4){
      names(outData)[names(outData)%in%c("ID","STRAT1","STRAT2","STRAT3")] <- c(idNmObs,popStrNm1,popStrNm2,popStrNm3)
      outData[,c(5:ncol(outData))] <- as.data.frame(lapply(outData[,c(5:ncol(outData))], FUN=function(x) out.digits(x,dig=4)))
    }
    
    
    # Subset table to print in the report
    if (case == 1) prnTab <- head(cbind(outData[,1:2], subset(outData, select = tabCol)),50)
    if (case == 2) prnTab <- head(cbind(outData[,1:3], subset(outData, select = tabCol)),50)
    if (case == 3) prnTab <- head(cbind(outData[,1:4], subset(outData, select = tabCol)),50)
    if (case == 4) prnTab <- head(cbind(outData[,1:5], subset(outData, select = tabCol)),50)
    
    streamsEnv <- parent.frame()
    if (exists("outData"))      assign("ncaOutput", outData, envir=streamsEnv)
    if (exists("statData"))     assign("ObsStat", statData, envir=streamsEnv)
    if (exists("statDataSim"))  assign("SimStat", statDataSim, envir=streamsEnv)
    if (exists("nmdf"))         assign("ncaSimData", nmdf, envir=streamsEnv)
    if (exists("dasdf"))        assign("ncaSimEst", dasdf, envir=streamsEnv)
    if (exists("concplot"))     assign("concPlot", concplot, envir=streamsEnv)
    if (exists("histobsplot"))  assign("histobsPlot", histobsplot, envir=streamsEnv)
    if (exists("popplot"))      assign("popPlot", popplot, envir=streamsEnv)
    
    if (!onlyNCA){
      if (exists("devplot"))      assign("devPlot", devplot, envir=streamsEnv)
      if (exists("outlierplot"))  assign("outlierPlot", outlierplot, envir=streamsEnv)
      if (exists("forestplot"))   assign("forestPlot", forestplot, envir=streamsEnv)
      if (exists("npdeplot"))     assign("npdePlot", npdeplot, envir=streamsEnv)
      if (exists("histnpdeplot")) assign("histnpdePlot", histnpdeplot, envir=streamsEnv)
    }
    
    if (printOut){
      # Create HTML output
      simFileNm <- ifelse(is.data.frame(simFile), deparse(substitute(simFile)), simFile)
      txt <- paste(txt,paste0("Name of the NONMEM simulation output file: \"",simFileNm,"\""),sep="\n")
      txt <- paste(txt,paste0("Number of simulations performed: ",nsim),sep="\n")
      
      if (is.null(outFileNm) || outFileNm==""){
        write.table(outData, file=paste0(usrdir,"/ncaOutput.tsv"), sep="\t", row.names=F, col.names=T, quote=F)
      }else{
        write.table(outData, file=paste0(usrdir,"/ncaOutput-",outFileNm,".tsv"), sep="\t", row.names=F, col.names=T, quote=F)
      }
      
      if (!onlyNCA){
        pddf <- cbind(pddf,OTL); names(pddf)[c((ncol(pddf)-1),ncol(pddf))] <- c("No. of outliers","Outlier IDs and outlying metrics")
        fnOut <- list(arglist=match.call(), case=case, TXT=txt, pddf=pddf, prnTab=prnTab, NSIM=nsim, spread=spread, conc=concplot, histobs=histobsplot, pop=popplot, dev=devplot, outlier=outlierplot, forest=forestplot, npde=npdeplot, histnpde=histnpdeplot, phth=phth, pwth=pwth)
      }else{
        fnOut <- list(arglist=match.call(), case=case, TXT=txt, pddf=pddf, prnTab=prnTab, NSIM=nsim, spread=spread, conc=concplot, histobs=histobsplot, pop=popplot, phth=phth, pwth=pwth)
      }
    }
  }
  #setwd(usrdir)
  
  if (printOut){
    misc <- system.file("misc", package = "ncappc")
    
    if (is.null(simFile) || onlyNCA){
      mdFile <- paste(misc,"ncappcReport-NCA.Rmd",sep="/")
      nwFile <- paste(misc,"ncappcReport-NCA.Rnw",sep="/")
      if (is.null(outFileNm) || outFileNm==""){
        outNm <- paste0("ncappcReport-NCA.tex")
      }else{
        outNm <- paste0("ncappcReport-NCA-",outFileNm,".tex")
      }
    }else{
      #if (!onlyNCA){
      mdFile <- paste(misc,"ncappcReport-NCA-PPC.Rmd",sep="/")
      nwFile <- paste(misc,"ncappcReport-NCA-PPC.Rnw",sep="/")
      if (is.null(outFileNm) || outFileNm==""){
        outNm <- paste0("ncappcReport-NCA-PPC.tex")
      }else{
        outNm <- paste0("ncappcReport-NCA-PPC-",outFileNm,".tex")
      }
      # }else{
      #   mdFile <- paste(misc,"ncappcReport-NCA-NOPPC.Rmd",sep="/")
      #   nwFile <- paste(misc,"ncappcReport-NCA-NOPPC.Rnw",sep="/")
      #   if (is.null(outFileNm) || outFileNm==""){
      #     outNm <- paste0("ncappcReport-NCA-NOPPC.tex")
      #   }else{
      #     outNm <- paste0("ncappcReport-NCA-NOPPC-",outFileNm,".tex")
      #   }
      # }
    }
    
    file.copy(paste(misc,"ncappc_report.Rmd",sep="/"),"ncappc_report.Rmd",overwrite = TRUE)
    file.copy(paste(misc,"styles.tex",sep="/"),"styles.tex",overwrite = TRUE)
    file.copy(paste(misc,"styles2.tex",sep="/"),"styles2.tex",overwrite = TRUE)
    file.copy(paste(misc,"custom.css",sep="/"),"custom.css",overwrite = TRUE)
    
    
    if(out_format == "html"){
      rmarkdown::render("ncappc_report.Rmd",output_format = "bookdown::html_document2")
    }  

    if(out_format == "pdf"){
      rmarkdown::render("ncappc_report.Rmd",output_format = "bookdown::pdf_document2")
    }  
    
    if(out_format == "all"){
      rmarkdown::render("ncappc_report.Rmd",output_format = "all")
    }
    
    if(out_format == "first"){
      rmarkdown::render("ncappc_report.Rmd")
    }
    
    # knit2html(input=mdFile, output=outNm, style=paste(misc,"custom.css",sep="/"))#, force_v1 = TRUE)
    
    # knit(input=nwFile, output=outNm)#, force_v1 = TRUE)
    # if (.Platform$OS.type == "unix"){
    #   texcomp <- system('which texi2pdf')
    #   if (texcomp == 0){
    #     knit2pdf(input=nwFile, output=outNm)#, force_v1 = TRUE)
    #   }else{
    #     print("Please install \"texi2pdf\" to compile the produced tex file into a PDF report")
    #   }
    # }else if (.Platform$OS.type == "windows"){
    #   texcomp <- system('kpsewhich pdftex --version')
    #   if (texcomp == 0){
    #     knit2pdf(input=nwFile, output=outNm)#, force_v1 = TRUE)
    #   }else{
    #     print("Please install \"pdftex\" to compile the produced tex file into a PDF report")
    #   }
    # }
  }
  
  unlink(list.files(pattern = "^.*ncappcReport-.*aux$|^.*ncappcReport-.*log$|^.*ncappcReport-.*out$|^.*ncappcReport-.*toc$|^.*ncappcReport-.*md$"))
  unlink(list.files(pattern = "sum.tex"))
  unlink(list.files(pattern = "tab.tex"))
  unlink(list.files(pattern = "knitr.sty"))
}
