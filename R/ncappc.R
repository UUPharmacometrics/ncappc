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
#' values of the arguments used in \pkg{ncappc} are shown in \strong{bold}.
#'
#' @param obsFile Observed concentration-time data from an internal data frame
#'   or an external table with comma, tab or space as separator
#' @param simFile NONMEM simulation output with the simulated concentration-time
#'   data from an internal data frame or an external table (\strong{"NULL"})
#' @param grNm Column name for population stratifier (\strong{"NULL"})
#' @param grp Stratification ID (e.g c(1,2)) (\strong{"NULL"})
#' @param flNm Column name for popualtion stratifier (\strong{"NULL"})
#' @param flag Stratification ID (e.g. c(1,2)) (\strong{"NULL"})
#' @param doseNm Column name to specify dose identifiers (\strong{"NULL"})
#' @param dose Dose identifiers to be used (c(1,2)) (\strong{"NULL"})
#' @param concUnit Unit of concentration ("ng/mL") (\strong{"[M].[L]^-3"})
#' @param timeUnit Unit of time ("h") (\strong{"[T]})
#' @param doseUnit Unit of dose amount ("ng") (\strong{"[M]"})
#' @param doseNormUnit Normalization factor of dose amount if used (kg)
#'   (\strong{"NULL"})
#' @param obsLog Concentration in observed data in logarithmic form (TRUE,
#'   FALSE) (\strong{"FALSE"})
#' @param simLog Concentration in simulated data in logarithmic form (TRUE,
#'   FALSE) (\strong{"FALSE"})
#' @param psnOut observed data is an output from PsN or in NONMEM output format 
#'   (TRUE, FALSE) (\strong{"FALSE"})
#' @param idNmObs Column name for ID in observed data (\strong{"ID"})
#' @param timeNmObs Column name for time in observed data (\strong{"TIME"})
#' @param concNmObs Column name for concentration in observed data
#'   (\strong{"DV"})
#' @param idNmSim Column name for ID in simulated data (\strong{"ID"})
#' @param timeNmSim Column name for time in simulated data (\strong{"TIME"})
#' @param concNmSim Column name for concentration in simulated data
#'   (\strong{"DV"})
#' @param AUCTimeRange User-defined window of time used to estimate AUC
#'   (\strong{"NULL"})
#' @param backExtrp If back-extrapolation is needed for AUC (TRUE or FALSE)
#'   (\strong{"FALSE"})
#' @param LambdaTimeRange User-defined window of time to estimate elimination 
#'   rate-constant (\strong{"NULL"})
#' @param LambdaExclude User-defined excluded observation time points for
#'   estimation of elimination rate-constant (\strong{"NULL"})
#' @param doseAmtNm Column name to specify dose amount (\strong{"NULL"})
#' @param adminType Route of administration
#'   (iv-bolus,iv-infusion,extravascular) (\strong{"extravascular"})
#' @param doseType Steady-state (ss) or nonsteady-state (ns) dose
#'   (\strong{"ns"})
#' @param Tau Dosing interval for steady-state data (\strong{"NULL"})
#' @param TI Infusion duration (\strong{"NULL"})
#' @param method linear, loglinear or mixed (\strong{"mixed"})
#' @param blqNm Name of BLQ column if used (\strong{"NULL"})
#' @param blqExcl Excluded BLQ value or logical condition (e.g. 1 or ">=1" or 
#'   c(1,">3")) (\strong{"1"})
#' @param evid Use EVID (TRUE, FALSE) (\strong{"FALSE"})
#' @param evidIncl Included EVID (\strong{"0"})
#' @param mdv Use MDV (TRUE(includes data for MDV==0), FALSE) (\strong{"FALSE"})
#' @param filterNm Column name for filter (\strong{"NULL"})
#' @param filterExcl Filter identifier or logical condition used for row
#'   exclusion (e.g. c(1, 2, "<20", ">=100", "!=100")) (\strong{"NULL"})
#' @param negConcExcl Exclude -ve conc (\strong{"FALSE"})
#' @param param NCA parameters (AUClast, AUClower_upper, AUCINF_obs, 
#'   AUCINF_pred, AUMClast, Cmax, Tmax, HL_Lambda_z) (c(\strong{"AUClast",
#'   "Cmax"}))
#' @param timeFormat time format (number, H:M, H:M:S) (\strong{"number"})
#' @param dateColNm colunm name for date if used (Date, DATE) (\strong{"NULL"})
#' @param dateFormat date format (D-M-Y, D/M/Y or any other combination of
#'   D,M,Y) (\strong{"NULL"})
#' @param spread Measure of the spread of simulated data (ppi (95\% parametric
#'   prediction interval) or npi (95\% nonparametric prediction interval))
#'   (\strong{"npi"})
#' @param tabCol Output columns to be printed in the report in addition to ID, 
#'   dose and population strata information (list of NCA metrics in a string 
#'   array) (\strong{c("AUClast", "Cmax", "Tmax", "AUCINF_obs", "Vz_obs",
#'   "Cl_obs", "HL_Lambda_z")})
#' @param figFormat format of the produced figures (bmp, jpeg, tiff, png)
#'   (\strong{"tiff"})
#' @param noPlot Perform only NCA calculations without any plot generation
#'   (TRUE, FALSE) (\strong{"FALSE"})
#' @param printOut Write/print output on the disk. No plot will be saved if
#'   noPlot is set to TRUE (TRUE, FALSE) (\strong{"TRUE"})
#' @param studyName Name of the study to be added as a description in the report
#'   (\strong{"NULL"})
#'
#' @import ggplot2
#' @import gridExtra
#' @import scales
#' @import gtable
#' @import knitr
#' @import xtable
#' @import reshape2
#' @import testthat
#' 
#' @return NCA results and diagnostic test results
#' @export
#'

ncappc <- function(obsFile=NULL,simFile=NULL,grNm=NULL,grp=NULL,flNm=NULL,flag=NULL,doseNm=NULL,dose=NULL,concUnit=NULL,timeUnit=NULL,doseUnit=NULL,doseNormUnit=NULL,obsLog="FALSE",simLog="FALSE",psnOut="FALSE",idNmObs="ID",timeNmObs="TIME",concNmObs="DV",idNmSim="ID",timeNmSim="TIME",concNmSim="DV",AUCTimeRange=NULL,backExtrp="FALSE",LambdaTimeRange=NULL,LambdaExclude=NULL,doseAmtNm=NULL,adminType="extravascular",doseType="ns",Tau=NULL,TI=NULL,method="mixed",blqNm=NULL,blqExcl=1,evid="FALSE",evidIncl=0,mdv="FALSE",filterNm=NULL,filterExcl=NULL,negConcExcl="FALSE",param=c("AUClast","Cmax"),timeFormat="number",dateColNm=NULL,dateFormat=NULL,spread="npi",tabCol=c("AUClast","Cmax","Tmax","AUCINF_obs","Vz_obs","Cl_obs","HL_Lambda_z"),figFormat="tiff",noPlot="FALSE",printOut="TRUE",studyName=NULL){
  
  "..density.." <- "meanObs" <- "sprlow" <- "sprhgh" <- "AUClast" <- "AUCINF_obs" <- "Cmax" <- "Tmax" <- "FCT" <- "ID" <- "GROUP" <- "FLAG" <- "NPDE" <- "mcil" <- "mciu" <- "sdu" <- "sducil" <- "sduciu" <- "scale_linetype_manual" <- "scale_color_manual" <- "xlab" <- "ylab" <- "guides" <- "guide_legend" <- "theme" <- "element_text" <- "unit" <- "element_rect" <- "geom_histogram" <- "aes" <- "geom_vline" <- "grid.arrange" <- "unit.c" <- "grid.grab" <- "ggsave" <- "facet_wrap" <- "ggplot" <- "labs" <- "geom_point" <- "geom_errorbarh" <- "knit2html" <- "knit2pdf" <- "knit" <- "file_test" <- "tail" <- "read.csv" <- "read.table" <- "dev.off" <- "write.table" <- "head" <- "write.csv" <- "coef" <- "dist" <- "lm" <- "median" <- "na.omit" <- "percent" <- "qchisq" <- "qnorm" <- "qt" <- "quantile" <- "scale_y_continuous" <- "sd" <- NULL
  rm(list=c("..density..","meanObs","sprlow","sprhgh","AUClast","AUCINF_obs","Cmax","Tmax","FCT","ID","GROUP","FLAG","NPDE","mcil","mciu","sdu","sducil","sduciu","scale_linetype_manual","scale_color_manual","xlab","ylab","guides","guide_legend","theme","element_text","unit","element_rect","geom_histogram","aes","geom_vline","grid.arrange","unit.c","grid.grab","ggsave","facet_wrap","ggplot","labs","geom_point","geom_errorbarh","knit2html","knit2pdf","knit","file_test","tail","read.csv","read.table","dev.off","write.table","head","write.csv","coef","dist","lm","median","na.omit","percent","qchisq","qnorm","qt","quantile","scale_y_continuous","sd"))
  
  options(warning.length=5000)
  options(scipen=999)
  usrdir <- getwd()
  
  # Observed data
  if (is.null(obsFile)){stop("Name of the file with observed data is required\n")}
  if (!is.data.frame(obsFile)){
    if (!file_test("-f", obsFile)){stop("File for the observed data does not exist\n")}
    # read observed data file
    if (psnOut == "FALSE"){
      extn <- tail(unlist(strsplit(obsFile, ".", fixed=T)), n=1)
      if(extn=="csv"){indf <- read.csv(obsFile)}else{indf <- read.table(obsFile, header=T)}
    }else if (psnOut == "TRUE"){
      indf <- read.table(obsFile, header=T, skip=1)
    }
  }else{
    indf <- obsFile
  }
  refdf <- indf # copy the input data before processing
  
  # Simulated data
  if ((!is.null(simFile)) && (!is.data.frame(simFile)) && (!file_test("-f", simFile))){stop("File for the simulated data does not exist\n")}
  
  # Check steady state dosing interval
  if (doseType == "ss" && (is.null(Tau))){setwd(usrdir);stop("Dosing interval is required for steady-state data\n")}
  
  # Check for column names
  if (idNmObs%in%colnames(indf)==F | timeNmObs%in%colnames(indf)==F | concNmObs%in%colnames(indf)==F){
    setwd(usrdir);stop("Incorrect column names of ID, TIME and/or DV\n")
  }else{
    idCol   <- which(colnames(indf) == idNmObs); id <- unique(indf[,idCol])
    timeCol <- which(colnames(indf) == timeNmObs)
    concCol <- which(colnames(indf) == concNmObs)
  }
  
  # exclude data based on specific values on filter column (optional)
  if (!is.null(filterNm)){
    if (filterNm%in%colnames(indf)==T & !is.null(filterExcl)){
      # filterExcl  == values to be excluded
      filterCol <- which(colnames(indf) == filterNm)
      for (i in 1:length(filterExcl)){
        if (grepl("^[-]?[0-9]*[.]?[0-9]*[eE]?[-]?[0-9]*[.]?[0-9]*$", filterExcl[i])){
          indf <- indf[indf[,filterCol] != filterExcl[i],]
        }else{
          indf <- eval(parse(text=paste("subset(indf, !",filterNm,"%in% indf[indf[,",filterCol,"]",filterExcl[i],",filterCol])",sep="")))
        }
      }
    }else{setwd(usrdir);stop("Incorrect filterNm or filterExcl specification\n")}
  }
  
  if (!is.null(grNm)){
    if (grNm%in%colnames(indf)==F){setwd(usrdir);stop("Incorrect name for the group column\n")}else{grCol <- which(colnames(indf) == grNm)}
    if (is.null(grp)){grp <- unique(sort(indf[,grCol]))}
    ngrp <- length(grp)
    for (i in 1:ngrp){
      if (nrow(indf[indf[,grCol]==grp[i],]) == 0){setwd(usrdir);stop("Group identifier does not match the group column\n")}
    }
  }
  
  if (!is.null(flNm)){
    if (flNm%in%colnames(indf)==F){setwd(usrdir);stop("Incorrect name for the flag column\n")}else{flCol <- which(colnames(indf) == flNm)}
    if (is.null(flag)){flag <- unique(indf[,flCol])}
    nflag <- length(flag)
    for (i in 1:nflag){
      if (nrow(indf[indf[,flCol]==flag[i],]) == 0){setwd(usrdir);stop("Flag identifier does not match the flag column\n")}
    }
  }
  
  # check time range, if any
  if ((!is.null(AUCTimeRange)) && (length(AUCTimeRange) != 2 | class(AUCTimeRange) != "numeric")){
    setwd(usrdir);stop("Incorrect time range for AUC calculation\n")
  }
  if ((!is.null(LambdaTimeRange)) && (length(LambdaTimeRange) != 2 | class(LambdaTimeRange) != "numeric")){
    setwd(usrdir);stop("Incorrect time range for Lambda calculation\n")
  }
  
  # check requirements for infusion data
  if (adminType == "iv-infusion" & is.null(TI) & ("AMT"%in%colnames(indf)==F | "RATE"%in%colnames(indf)==F)){setwd(usrdir);stop("Duration of the infusion time is needed if AMT and RATE are absent in the input data\n")}
  
  # Set backExtrp to FALSE in the presence of simulated data
  if (!is.null(simFile)){backExtrp <- "FALSE"}
  
  # Dose identifiers. For missing doseNm argument, data is assumed to have single dose.
  # Dose amount is extracted from doseAmtNm column
  if (!is.null(doseNm)){
    if (doseNm%in%colnames(indf)==F){setwd(usrdir);stop("Incorrect name for the dose column\n")}else{doseCol <- which(colnames(indf) == doseNm)}
    oidNm <- doseNm
    if (is.null(dose)){dose <- unique(sort(indf[,doseCol]))}
    ndose <- length(dose)
    for (i in 1:ndose){
      if (nrow(indf[indf[,doseCol]==dose[i],]) == 0){setwd(usrdir);stop("Dose identifier does not match the dose column\n")}
    }
  }else{
    oidNm <- "OID"; dose <- 1; ndose <- 1
  }
  
  if (is.null(doseAmtNm)){
    #if ("AMT"%in%colnames(indf) == F){setwd(usrdir);stop("Dose amount column is required as AMT column is absent\n")}
    if ("AMT"%in%colnames(indf)) doseAmtNm <- "AMT"
  }
  
  # Dose unit
  if (is.null(doseUnit)) doseUnit <- "[M]"
  
  # Preliminary description
  # Units for dose, time and conc
  dunit <- ifelse(is.null(doseNormUnit), doseUnit, paste(doseUnit,"/",doseNormUnit))
  tunit <- ifelse(is.null(timeUnit), "T", timeUnit)
  cunit <- ifelse(is.null(concUnit), "M.L^-3", concUnit)
  
  # ignore data with BLQ = 1 or user specified value (optional)
  if (!is.null(blqNm)){
    if (blqNm%in%colnames(indf) == T){
      blqCol <- which(colnames(indf) == blqNm)
      for (i in 1:length(blqExcl)){
        if (grepl("^[-]?[0-9]*[.]?[0-9]*[eE]?[-]?[0-9]*[.]?[0-9]*$", blqExcl[i]) == T) {indf <- indf[indf[,blqNm] != blqExcl[i],]}
        if (grepl("^[-]?[0-9]*[.]?[0-9]*[eE]?[-]?[0-9]*[.]?[0-9]*$", blqExcl[i]) == F) {indf <- eval(parse(text=paste("subset(indf, !",blqNm,"%in% indf[as.numeric(as.character(indf[,",blqCol,"]))",blqExcl[i],",blqCol])",sep="")))}
      }
    }else{setwd(usrdir);stop("Incorrect BLQ column name\n")}
  }
  
  # include data based on specific values on EVID column (optional) but keep rows with TIME == 0
  if ((evid == "TRUE") & ("EVID"%in%colnames(indf) == T)){
    # uevid == unique values in EVID column
    # evidIncl == EVID values to be included
    # ievid == EVID values to be ignored
    uevid <- unique(as.numeric(as.character(indf$EVID))); ievid <- setdiff(uevid, as.numeric(evidIncl))
    if (length(ievid) != 0){for (i in 1:length(ievid)){indf <- indf[indf$EVID != ievid[i],]}}
  }
  
  # if MDV fiter is present, exclude data for MDV == 1 but keep rows with TIME == 0
  if (mdv == "TRUE"){
    if ("MDV"%in%colnames(indf) == T){indf <- indf[indf$MDV == 0,]}
  }
  
  # Appropriate date format
  if (!is.null(dateColNm)){
    if (dateColNm%in%colnames(indf)==F){setwd(usrdir);stop("Incorrect name for the date column\n")}
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
    }else{setwd(usrdir);stop("Incorrect date format. Currently allowed date formats are \"D-M-Y\", \"D/M/Y\", or any combination of D, M and Y separated by either - or /\n")}
  }
  
  # Appropriate time format
  if (timeFormat != "number"){
    if(timeFormat == "H:M"){
      timeFormat <- "%H:%M"
    }else if(timeFormat == "H:M:S"){
      timeFormat <- "%H:%M:%S"
    }else{setwd(usrdir);stop("Incorrect time format. Currently allowed date formats are \"number\", \"H:M\", \"H:M:S\"\n")}
  }
  
  # create 4 possible cases based on input file columns
  # case 1: "GROUP" = FALSE  and "FLAG" = FALSE
  # case 2: "GROUP" = TRUE and "FLAG" = FALSE
  # case 3: "GROUP" = FALSE  and "FLAG" = TRUE
  # case 4: "GROUP" = TRUE and "FLAG" = TRUE
  
  outData <- data.frame()
  
  if(is.null(grNm)  & is.null(flNm))  case <- 1
  if(!is.null(grNm) & is.null(flNm))  case <- 2
  if(is.null(grNm)  & !is.null(flNm)) case <- 3
  if(!is.null(grNm) & !is.null(flNm)) case <- 4
  
  
  # Allowed NCA parameters
  alwprm <- c("AUClast","AUClower_upper","AUCINF_obs","AUCINF_pred","AUMClast","Cmax","Tmax","HL_Lambda_z")
  npr    <- length(param)
  fctNm  <- data.frame()
  for (p in 1:npr){
    if (param[p]%in%alwprm == F){setwd(usrdir);stop("Incorrect NCA parameters. Please select NCA parameters from \"AUClast\", \"AUClower_upper\", \"AUCINF_obs\", \"AUCINF_pred\", \"AUMClast\", \"Cmax\", \"Tmax\", \"HL_Lambda_z\"\n")}
    if (param[p] == "AUClast" | param[p] == "AUClower_upper" | param[p] == "AUCINF_obs" | param[p] == "AUCINF_pred"){
      fctNm <- rbind(fctNm, data.frame(prmNm=param[p],prmUnit=paste(param[p]," (",cunit,"*",tunit,")",sep="")))
    }else if (param[p] == "AUMClast"){
      fctNm <- rbind(fctNm, data.frame(prmNm=param[p],prmUnit=paste(param[p]," (",cunit,"*",tunit,"^2)",sep="")))
    }else if (param[p] == "Cmax"){
      fctNm <- rbind(fctNm, data.frame(prmNm=param[p],prmUnit=paste(param[p]," (",cunit,")",sep="")))
    }else if (param[p] == "Tmax"){
      fctNm <- rbind(fctNm, data.frame(prmNm=param[p],prmUnit=paste(param[p]," (",tunit,")",sep="")))
    }else if (param[p] == "HL_Lambda_z"){
      fctNm <- rbind(fctNm, data.frame(prmNm=param[p],prmUnit=paste(param[p]," (",tunit,")",sep="")))
    }
  }
  
  obsFileNm <- ifelse(is.data.frame(obsFile), deparse(substitute(obsFile)), obsFile)
  
  if(!is.null(studyName)){
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
  txt <- paste(txt,paste0("No. of dosing occasions: ",ndose),sep="\n")
  if (ndose == 1){
    txt <- paste(txt,"Occasion flag or OID: OID",sep="\n")
    txt <- paste(txt,"Occasion ID: 1",sep="\n")
  }else if (ndose > 1){
    txt <- paste(txt,paste0("Occasion flag: ",oidNm),sep="\n")
    txt <- paste(txt,paste0("Occasion ID (OID): ",paste(dose,collapse=", ")),sep="\n")
  }
  txt <- paste(txt,paste0("Total number of subjects in the input file = ",length(unique(indf[,idNmObs]))),sep="\n")
  if (case == 2){
    txt <- paste(txt,paste0("Population group flag or GID: ",grNm),sep="\n")
    txt <- paste(txt,paste0("Population group ID: ",paste(grp,collapse=", ")),sep="\n")
  }
  if (case == 3){
    txt <- paste(txt,paste0("Population group flag or GID: ",flNm),sep="\n")
    txt <- paste(txt,paste0("Population group ID: ",paste(flag,collapse=", ")),sep="\n")
  }
  if (case == 4){
    txt <- paste(txt,paste0("Population group flag or GID: ",grNm),sep="\n")
    txt <- paste(txt,paste0("Population group ID: ",paste(grp,collapse=", ")),sep="\n")
    txt <- paste(txt,paste0("Population sub-group flag or SGID: ",flNm),sep="\n")
    txt <- paste(txt,paste0("Population sub-group ID: ",paste(flag,collapse=", ")),sep="\n")
  }
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # Set plot output dimensions
  if (npr<=2){ncol<-2;hth<-12;wth<-16;phth<-6;pwth<-7}else if(npr>2 & npr<=4){ncol<-2;hth<-18;wth<-18;phth<-9;pwth<-7}else if(npr>4 & npr<=6){ncol<-3;hth<-18;wth<-25;phth<-9;pwth<-9}else if(npr>6){ncol<-3;hth<-26;wth<-25;phth<-10;pwth<-9}
  
  # Initiate plot lists for pdf output
  concplot <- list(); histobsplot=list(); popplot <- list(); devplot <- list(); outlierplot <- list(); forestplot <- list(); npdeplot <- list(); histnpdeplot <- list()
  
  # calculate the NCA parameters for the observed data
  dset = "obs"
  
  # Function to generate time and conc data to calculate NCA metrics
  ncaId <- function(ifdf,ID){
    if(adminType == "iv-infusion" & is.null(TI)){
      amt  <- ifdf[ifdf[,idCol]==ID & ifdf$AMT > 0,"AMT"][1]
      rate <- ifdf[ifdf[,idCol]==ID & ifdf$RATE > 0,"RATE"][1]
      if (is.na(amt) | is.na(rate) | rate==0){setwd(usrdir);stop(paste("Incorrect AMT and/or RATE value for ",ID,sep=""))}else{TI <- amt/rate}
    }else{
      TI <- "NaN"
    }
    if (timeFormat != "number"){
      time <- numeric(0)
      if (!is.null(dateColNm)){
        tm <- as.POSIXct(paste(ifdf[ifdf[,idCol]==ID,dateColNm],ifdf[ifdf[,idCol]==ID,timeCol]),format=paste(dateFormat,timeFormat,sep=" "))
      }else{
        tm <- ifdf[ifdf[,idCol]==ID,timeCol]
      }
      for (j in 1:length(tm)){
        time[j] <- ifelse ((is.null(dateColNm)), as.numeric(difftime(strptime(tm[j], format=timeFormat), strptime(tm[1], format=timeFormat), units='hours')), as.numeric(difftime(strptime(tm[j], format="%Y-%m-%d %H:%M:%S"), strptime(tm[1], format="%Y-%m-%d %H:%M:%S"), units='hours')))
      }
    }else{
      time <- suppressWarnings(as.numeric(as.character(ifdf[ifdf[,idCol]==ID,timeCol])))
    }
    tconc <- ifdf[ifdf[,idCol]==ID, concCol]; conc <- numeric(0)
    for (c in 1:length(tconc)){conc <- c(conc, ifelse ((tconc[c]=="."), 0, as.numeric(as.character(tconc[c]))))}; rm(tconc)
    if (length(which(is.na(time)))!=0){
      zidx <- which(is.na(time))
      time <- time[-zidx]
      conc <- conc[-zidx]
    }
    if (length(which(is.na(conc)))!=0){
      zidx <- which(is.na(conc))
      time <- time[-zidx]
      conc <- conc[-zidx]
    }
    if (obsLog == "TRUE"){
      for (c in 1:length(conc)){conc[c] <- ifelse ((conc[c] != 0), exp(conc[c]), 0)}
    }
    
    tc <- data.frame(time,conc)
    if(nrow(tc)>0){
      tc <- tc[order(tc$time),]
      tc$time <- tc$time + abs(min(tc$time))
    }
    return(tc)
  }
  
  cdata <- data.frame(Time=numeric(0),Conc=numeric(0),ID=character(0),FCT=character(0))
  pddf  <- data.frame(a=character(0),b=character(0),c=character(0))
  if (case == 1){
    for (d in 1:ndose){
      if (!is.null(doseNm)){ifdf <- indf[indf[,doseNm]==dose[d],]}else{ifdf <- indf}
      if (length(which(is.na(ifdf[,concCol]) | ifdf[,concCol]=="")) != 0){ifdf <- ifdf[-which(is.na(ifdf[,concCol]) | ifdf[,concCol]==""),]}
      if (nrow(ifdf) == 0){next}
      idd        <- unique(ifdf[,idCol])
      DoseNumber <- dose[d]
      if (!is.null(doseAmtNm) & is.null(doseNm)){
        doseAmount <- paste(as.numeric(unique(refdf[refdf[,doseAmtNm]!="." & as.numeric(as.character(refdf[,doseAmtNm])) > 0, doseAmtNm])), collapse=", ")
      }else if(!is.null(doseAmtNm) & !is.null(doseNm)){
        doseAmount <- paste(as.numeric(unique(refdf[refdf[,doseNm]==dose[d] & refdf[,doseAmtNm]!="." & as.numeric(as.character(refdf[,doseAmtNm])) > 0, doseAmtNm])), collapse=", ")
      }else{
        doseAmount <- NA
      }
      # Description
      pddf <- rbind(pddf, data.frame(a=DoseNumber, b=doseAmount, c=length(idd)))
      for (i in 1:length(idd)){
        if (!is.null(doseAmtNm) & is.null(doseNm)){
          idzAmt <- as.numeric(refdf[refdf[,idNmObs]==idd[i] & refdf[,doseAmtNm]!="." & as.numeric(as.character(refdf[,doseAmtNm])) > 0, doseAmtNm][1])
        }else if(!is.null(doseAmtNm) & !is.null(doseNm)){
          idzAmt <- as.numeric(refdf[refdf[,doseNm]==dose[d] & refdf[,idNmObs]==idd[i] & refdf[,doseAmtNm]!="." & as.numeric(as.character(refdf[,doseAmtNm])) > 0, doseAmtNm][1])
        }else{
          idzAmt <- NA
        }
        tc   <- ncaId(ifdf,idd[i])
        if(nrow(tc)==0) next
        time <- as.numeric(tc$time)
        conc <- as.numeric(tc$conc)
        cdata  <- rbind(cdata,cbind(Time=time,Conc=conc,ID=as.character(idd[i]),FCT=paste(oidNm,"-",dose[d],sep="")))
        NCAprm <- est.nca(time=time,conc=conc,backExtrp=backExtrp,negConcExcl=negConcExcl,doseType=doseType,adminType=adminType,doseNm=doseNm,dose=dose,doseNumber=DoseNumber,doseAmt=idzAmt,method=method,AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,LambdaExclude=LambdaExclude,Tau=Tau,TI=TI,simFile=simFile,dset=dset) # calls est.nca function
        outData <- rbind(outData, data.frame(ID=idd[i],DoseNumber=DoseNumber,DoseAmount=idzAmt,t(NCAprm)))
      }
      if(noPlot=="FALSE"){
        plotData    <- subset(outData, DoseNumber=dose[d], select=c(AUClast,AUCINF_obs,Cmax,Tmax))
        if (nrow(plotData)<=5) next
        pltPrm      <- c("AUClast","AUCINF_obs","Cmax","Tmax")
        for (p in 1:length(pltPrm)){if (nrow(plotData[plotData[,p] != "NaN",])<5) pltPrm <- pltPrm[-p]}
        if (length(pltPrm) == 0) next
        figlbl      <- paste(oidNm,"-",dose[d],sep="")
        histobsgrob <- histobs.plot(plotData=plotData,figlbl=figlbl,param=pltPrm,cunit=cunit,tunit=tunit,spread=spread)
        gdr         <- histobsgrob$gdr
        mylegend    <- histobsgrob$legend
        lheight     <- histobsgrob$lheight
        if (printOut=="TRUE"){
          fl <- paste(usrdir,"/HistObs_",figlbl,sep="")
          eval(parse(text=paste(figFormat,"(file=\"",fl,".",figFormat,"\",height=15,width=14,units=\"cm\",res=200)",sep="")))
          suppressMessages(suppressWarnings(grid.arrange(gdr, mylegend, heights = unit.c(unit(1,"npc")-lheight, lheight))))
          dev.off()
        }
        suppressMessages(suppressWarnings(grid.arrange(gdr, mylegend, heights = unit.c(unit(1,"npc")-lheight, lheight))))
        ggr <- grid.grab()
        histobsplot[[length(histobsplot)+1]] <- ggr
      }
    }
    cnm         <- c(paste("OID (",oidNm,")",sep=""),paste("Dose (",dunit,")",sep=""),"No. of individuals")
    names(pddf) <- cnm
    
    if(noPlot=="FALSE"){
      fct <- unique(cdata$FCT); nfct <- length(fct)
      # ggplot for conc vs. time
      for (p in seq(1,nfct,3)){
        str <- ""
        if ((nfct-p)>=2) {
          dht   <- 16
          dwd   <- 15
          str   <- gsub("^_", "", paste(str,fct[p:(p+2)],sep="_",collapse=""))
          pdata <- subset(cdata, FCT == fct[p] | FCT == fct[p+1] | FCT == fct[p+2])
        }else if ((nfct-p)==1){
          dht   <- 11
          dwd   <- 13
          str   <- gsub("^_", "", paste(str,fct[p:(p+1)],sep="_",collapse=""))
          pdata <- subset(cdata, FCT == fct[p] | FCT == fct[p+1])
        }else if ((nfct-p)==0){
          dht   <- 7
          dwd   <- 11
          str   <- gsub("^_", "", paste(str,fct[p],sep="_",collapse=""))
          pdata <- subset(cdata, FCT == fct[p])
        }
        gdr <- dv.plot(pdata=pdata,cunit=cunit,tunit=tunit)
        suppressMessages(suppressWarnings(grid.arrange(gdr)))
        ggr <- grid.grab()
        concplot[[length(concplot)+1]] <- ggr
        if (printOut=="TRUE") ggsave(filename=paste(usrdir,"/TimeConc_",str,".",figFormat,sep=""),plot=gdr,height=dht,width=dwd,units="cm",dpi=200)
      }
    }
  }
  if (case == 2){
    for (g in 1:ngrp){
      for (d in 1:ndose){
        if (!is.null(doseNm)){ifdf <- indf[indf[,grCol]==as.character(grp[g]) & indf[,doseNm]==dose[d],]}else{ifdf <- indf[indf[,grCol]==as.character(grp[g]),]}
        if (length(which(is.na(ifdf[,concCol]) | ifdf[,concCol]=="")) != 0){ifdf <- ifdf[-which(is.na(ifdf[,concCol]) | ifdf[,concCol]==""),]}
        if (nrow(ifdf) == 0){next}
        idd        <- unique(ifdf[,idCol])
        DoseNumber <- dose[d]
        igr        <- as.character(grp[g])
        if (!is.null(doseAmtNm) & is.null(doseNm)){
          doseAmount <- paste(as.numeric(unique(refdf[refdf[,grCol]==as.character(grp[g]) & refdf[,doseAmtNm]!="." & as.numeric(as.character(refdf[,doseAmtNm])) > 0, doseAmtNm])), collapse=", ")
        }else if(!is.null(doseAmtNm) & !is.null(doseNm)){
          doseAmount <- paste(as.numeric(unique(refdf[refdf[,grCol]==as.character(grp[g]) & refdf[,doseNm]==dose[d] & refdf[,doseAmtNm]!="." & as.numeric(as.character(refdf[,doseAmtNm])) > 0, doseAmtNm])), collapse=", ")
        }else{
          doseAmount <- NA
        }
        # Description
        pddf <- rbind(pddf, data.frame(a=as.character(grp[g]), b=DoseNumber, c=doseAmount, d=length(idd)))
        for (i in 1:length(idd)){
          if (!is.null(doseAmtNm) & is.null(doseNm)){
            idzAmt <- as.numeric(refdf[refdf[,grCol]==as.character(grp[g]) & refdf[,idNmObs]==idd[i] & refdf[,doseAmtNm]!="." & as.numeric(as.character(refdf[,doseAmtNm])) > 0, doseAmtNm][1])
          }else if(!is.null(doseAmtNm) & !is.null(doseNm)){
            idzAmt <- as.numeric(refdf[refdf[,grCol]==as.character(grp[g]) & refdf[,doseNm]==dose[d] & refdf[,idNmObs]==idd[i] & refdf[,doseAmtNm]!="." & as.numeric(as.character(refdf[,doseAmtNm])) > 0, doseAmtNm][1])
          }else{
            idzAmt <- NA
          }
          tc   <- ncaId(ifdf,idd[i])
          if(nrow(tc)==0) next
          time <- as.numeric(tc$time)
          conc <- as.numeric(tc$conc)
          cdata  <- rbind(cdata,cbind(Time=time,Conc=conc,ID=as.character(idd[i]),FCT=paste(grNm,"-",as.character(grp[g]),"_",oidNm,"-",dose[d],sep="")))
          NCAprm <- est.nca(time=time,conc=conc,backExtrp=backExtrp,negConcExcl=negConcExcl,doseType=doseType,adminType=adminType,doseNm=doseNm,dose=dose,doseNumber=DoseNumber,doseAmt=idzAmt,method=method,AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,LambdaExclude=LambdaExclude,Tau=Tau,TI=TI,simFile=simFile,dset=dset) # calls est.nca function
          outData <- rbind(outData, cbind(data.frame(GROUP=igr),data.frame(ID=idd[i],DoseNumber=DoseNumber,DoseAmount=idzAmt,t(NCAprm))))
        }
        
        if(noPlot=="FALSE"){
          plotData    <- subset(outData, DoseNumber==dose[d], select=c(AUClast,AUCINF_obs,Cmax,Tmax))
          if (nrow(plotData)<=5) next
          pltPrm      <- c("AUClast","AUCINF_obs","Cmax","Tmax")
          for (p in 1:length(pltPrm)){if (nrow(plotData[plotData[,p] != "NaN",])<5) pltPrm <- pltPrm[-p]}
          if (length(pltPrm) == 0) next
          figlbl      <- paste(grNm,"-",as.character(grp[g]),"_",oidNm,"-",dose[d],sep="")
          histobsgrob <- histobs.plot(plotData=plotData,figlbl=figlbl,param=c("AUClast","AUCINF_obs","Cmax","Tmax"),cunit=cunit,tunit=tunit,spread=spread)
          gdr         <- histobsgrob$gdr
          mylegend    <- histobsgrob$legend
          lheight     <- histobsgrob$lheight
          if (printOut=="TRUE"){
            fl <- paste(usrdir,"/HistObs_",figlbl,sep="")
            eval(parse(text=paste(figFormat,"(file=\"",fl,".",figFormat,"\",height=15,width=14,units=\"cm\",res=200)",sep="")))
            suppressMessages(suppressWarnings(grid.arrange(gdr, mylegend, heights = unit.c(unit(1,"npc")-lheight, lheight))))
            dev.off()
          }
          suppressMessages(suppressWarnings(grid.arrange(gdr, mylegend, heights = unit.c(unit(1,"npc")-lheight, lheight))))
          ggr <- grid.grab()
          histobsplot[[length(histobsplot)+1]] <- ggr
        }
      }
    }
    cnm         <- c(paste("GID (",grNm,")",sep=""),paste("OID (",oidNm,")",sep=""),paste("Dose (",dunit,")",sep=""),"No. of individuals")
    names(pddf) <- cnm
    
    if(noPlot=="FALSE"){
      # ggplot for conc vs. time
      fct <- unique(cdata$FCT); nfct <- length(fct)
      # ggplot for conc vs. time
      for (p in seq(1,nfct,3)){
        str <- ""
        if ((nfct-p)>=2) {
          dht   <- 16
          dwd   <- 15
          str   <- gsub("^_", "", paste(str,fct[p:(p+2)],sep="_",collapse=""))
          pdata <- subset(cdata, FCT == fct[p] | FCT == fct[p+1] | FCT == fct[p+2])
        }else if ((nfct-p)==1){
          dht   <- 11
          dwd   <- 13
          str   <- gsub("^_", "", paste(str,fct[p:(p+1)],sep="_",collapse=""))
          pdata <- subset(cdata, FCT == fct[p] | FCT == fct[p+1])
        }else if ((nfct-p)==0){
          dht   <- 7
          dwd   <- 11
          str   <- gsub("^_", "", paste(str,fct[p],sep="_",collapse=""))
          pdata <- subset(cdata, FCT == fct[p])
        }
        gdr <- dv.plot(pdata=pdata,cunit=cunit,tunit=tunit)
        suppressMessages(suppressWarnings(grid.arrange(gdr)))
        ggr <- grid.grab()
        concplot[[length(concplot)+1]] <- ggr
        if (printOut=="TRUE") ggsave(filename=paste(usrdir,"/TimeConc_",str,".",figFormat,sep=""),plot=gdr,height=dht,width=dwd,units="cm",dpi=200)
      } 
    }
  }
  if (case == 3){
    for (f in 1:nflag){
      for (d in 1:ndose){
        if (!is.null(doseNm)){ifdf <- indf[indf[,flCol]==as.character(flag[f]) & indf[,doseNm]==dose[d],]}else{ifdf <- indf[indf[,flCol]==as.character(flag[f]),]}
        if (length(which(is.na(ifdf[,concCol]) | ifdf[,concCol]=="")) != 0){ifdf <- ifdf[-which(is.na(ifdf[,concCol]) | ifdf[,concCol]==""),]}
        if (nrow(ifdf) == 0){next}
        idd        <- unique(ifdf[,idCol])
        DoseNumber <- dose[d]
        iflag      <- as.character(flag[f])
        if (!is.null(doseAmtNm) & is.null(doseNm)){
          doseAmount <- paste(as.numeric(unique(refdf[refdf[,flCol]==as.character(flag[f]) & refdf[,doseAmtNm]!="." & as.numeric(as.character(refdf[,doseAmtNm])) > 0, doseAmtNm])), collapse=", ")
        }else if(!is.null(doseAmtNm) & !is.null(doseNm)){
          doseAmount <- paste(as.numeric(unique(refdf[refdf[,flCol]==as.character(flag[f]) & refdf[,doseNm]==dose[d] & refdf[,doseAmtNm]!="." & as.numeric(as.character(refdf[,doseAmtNm])) > 0, doseAmtNm])), collapse=", ")
        }else{
          doseAmount <- NA
        }
        # Description
        pddf   <- rbind(pddf, data.frame(a=as.character(flag[f]), b=DoseNumber, c=doseAmount, d=length(idd)))
        for (i in 1:length(idd)){
          if (!is.null(doseAmtNm) & is.null(doseNm)){
            idzAmt <- as.numeric(refdf[refdf[,flCol]==as.character(flag[f]) & refdf[,idNmObs]==idd[i] & refdf[,doseAmtNm]!="." & as.numeric(as.character(refdf[,doseAmtNm])) > 0, doseAmtNm][1])
          }else if(!is.null(doseAmtNm) & !is.null(doseNm)){
            idzAmt <- as.numeric(refdf[refdf[,flCol]==as.character(flag[f]) & refdf[,doseNm]==dose[d] & refdf[,idNmObs]==idd[i] & refdf[,doseAmtNm]!="." & as.numeric(as.character(refdf[,doseAmtNm])) > 0, doseAmtNm][1])
          }else{
            idzAmt <- NA
          }
          tc   <- ncaId(ifdf,idd[i])
          if(nrow(tc)==0) next
          time <- as.numeric(tc$time)
          conc <- as.numeric(tc$conc)
          cdata  <- rbind(cdata,cbind(Time=time,Conc=conc,ID=as.character(idd[i]),FCT=paste(flNm,"-",as.character(flag[f]),"_",oidNm,"-",dose[d],sep="")))
          NCAprm <- est.nca(time=time,conc=conc,backExtrp=backExtrp,negConcExcl=negConcExcl,doseType=doseType,adminType=adminType,doseNm=doseNm,dose=dose,doseNumber=DoseNumber,doseAmt=idzAmt,method=method,AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,LambdaExclude=LambdaExclude,Tau=Tau,TI=TI,simFile=simFile,dset=dset) # calls est.nca function
          outData <- rbind(outData, cbind(data.frame(FLAG=iflag),data.frame(ID=idd[i],DoseNumber=DoseNumber,DoseAmount=idzAmt,t(NCAprm))))
        }
        
        if(noPlot=="FALSE"){
          plotData    <- subset(outData, DoseNumber==dose[d], select=c(AUClast,AUCINF_obs,Cmax,Tmax))
          if (nrow(plotData)<=5) next
          pltPrm      <- c("AUClast","AUCINF_obs","Cmax","Tmax")
          for (p in 1:length(pltPrm)){ if (nrow(plotData[plotData[,p] != "NaN",])<5) pltPrm <- pltPrm[-p]}
          if (length(pltPrm) == 0) next
          figlbl      <- paste(flNm,"-",as.character(flag[f]),"_",oidNm,"-",dose[d],sep="")
          histobsgrob <- histobs.plot(plotData=plotData,figlbl=figlbl,param=c("AUClast","AUCINF_obs","Cmax","Tmax"),cunit=cunit,tunit=tunit,spread=spread)
          gdr         <- histobsgrob$gdr
          mylegend    <- histobsgrob$legend
          lheight     <- histobsgrob$lheight
          if (printOut=="TRUE"){
            fl <- paste(usrdir,"/HistObs_",figlbl,sep="")
            eval(parse(text=paste(figFormat,"(file=\"",fl,".",figFormat,"\",height=15,width=14,units=\"cm\",res=200)",sep="")))
            suppressMessages(suppressWarnings(grid.arrange(gdr, mylegend, heights = unit.c(unit(1,"npc")-lheight, lheight))))
            dev.off()
          }
          suppressMessages(suppressWarnings(grid.arrange(gdr, mylegend, heights = unit.c(unit(1,"npc")-lheight, lheight))))
          ggr <- grid.grab()
          histobsplot[[length(histobsplot)+1]] <- ggr
        }
      }
    }
    cnm         <- c(paste("GID (",flNm,")",sep=""),paste("OID (",oidNm,")",sep=""),paste("Dose (",dunit,")",sep=""),"No. of individuals")
    names(pddf) <- cnm
    
    if(noPlot=="FALSE"){
      # ggplot for conc vs. time
      fct <- unique(cdata$FCT); nfct <- length(fct)
      # ggplot for conc vs. time
      for (p in seq(1,nfct,3)){
        str <- ""
        if ((nfct-p)>=2) {
          dht   <- 16
          dwd   <- 15
          str   <- gsub("^_", "", paste(str,fct[p:(p+2)],sep="_",collapse=""))
          pdata <- subset(cdata, FCT == fct[p] | FCT == fct[p+1] | FCT == fct[p+2])
        }else if ((nfct-p)==1){
          dht   <- 11
          dwd   <- 13
          str   <- gsub("^_", "", paste(str,fct[p:(p+1)],sep="_",collapse=""))
          pdata <- subset(cdata, FCT == fct[p] | FCT == fct[p+1])
        }else if ((nfct-p)==0){
          dht   <- 7
          dwd   <- 11
          str   <- gsub("^_", "", paste(str,fct[p],sep="_",collapse=""))
          pdata <- subset(cdata, FCT == fct[p])
        }
        gdr <- dv.plot(pdata=pdata,cunit=cunit,tunit=tunit)
        suppressMessages(suppressWarnings(grid.arrange(gdr)))
        ggr <- grid.grab()
        concplot[[length(concplot)+1]] <- ggr
        if (printOut=="TRUE") ggsave(filename=paste(usrdir,"/TimeConc_",str,".",figFormat,sep=""),plot=gdr,height=dht,width=dwd,units="cm",dpi=200)
      }
    }
  }
  if (case == 4){
    for (g in 1:ngrp){
      for (f in 1:nflag){
        for (d in 1:ndose){
          if (!is.null(doseNm)){ifdf <- indf[(indf[,flCol]==as.character(flag[f]) & indf[,grCol]==as.character(grp[g]) & indf[,doseNm]==dose[d]),]}else{ifdf <- indf[(indf[,flCol]==as.character(flag[f]) & indf[,grCol]==as.character(grp[g])),]}
          if (length(which(is.na(ifdf[,concCol]) | ifdf[,concCol]=="")) != 0){ifdf <- ifdf[-which(is.na(ifdf[,concCol]) | ifdf[,concCol]==""),]}
          if (nrow(ifdf) == 0){next}
          idd        <- unique(ifdf[,idCol])
          DoseNumber <- dose[d]
          igr        <- as.character(grp[g])
          iflag      <- as.character(flag[f])
          if (!is.null(doseAmtNm) & is.null(doseNm)){
            doseAmount <- paste(as.numeric(unique(refdf[refdf[,grCol]==as.character(grp[g]) & refdf[,flCol]==as.character(flag[f]) & refdf[,doseAmtNm]!="." & as.numeric(as.character(refdf[,doseAmtNm])) > 0, doseAmtNm])), collapse=", ")
          }else if(!is.null(doseAmtNm) & !is.null(doseNm)){
            doseAmount <- paste(as.numeric(unique(refdf[refdf[,grCol]==as.character(grp[g]) & refdf[,flCol]==as.character(flag[f]) & refdf[,doseNm]==dose[d] & refdf[,doseAmtNm]!="." & as.numeric(as.character(refdf[,doseAmtNm])) > 0, doseAmtNm])), collapse=", ")
          }else{
            doseAmount <- NA
          }
          # Description
          pddf  <- rbind(pddf, data.frame(a=as.character(grp[g]), b=as.character(flag[f]), c=DoseNumber, d=doseAmount, e=length(idd)))
          for (i in 1:length(idd)){
            if (!is.null(doseAmtNm) & is.null(doseNm)){
              idzAmt <- as.numeric(refdf[refdf[,grCol]==as.character(grp[g]) & refdf[,flCol]==as.character(flag[f]) & refdf[,idNmObs]==idd[i] & refdf[,doseAmtNm]!="." & as.numeric(as.character(refdf[,doseAmtNm])) > 0, doseAmtNm][1])
            }else if(!is.null(doseAmtNm) & !is.null(doseNm)){
              idzAmt <- as.numeric(refdf[refdf[,grCol]==as.character(grp[g]) & refdf[,flCol]==as.character(flag[f]) & refdf[,doseNm]==dose[d] & refdf[,idNmObs]==idd[i] & refdf[,doseAmtNm]!="." & as.numeric(as.character(refdf[,doseAmtNm])) > 0, doseAmtNm][1])
            }else{
              idzAmt <- NA
            }
            tc     <- ncaId(ifdf,idd[i])
            if(nrow(tc)==0) next
            time   <- as.numeric(tc$time)
            conc   <- as.numeric(tc$conc)
            cdata  <- rbind(cdata,cbind(Time=time,Conc=conc,ID=as.character(idd[i]),FCT=paste(grNm,"-",as.character(grp[g]),"_",flNm,"-",as.character(flag[f]),"_",oidNm,"-",dose[d],sep="")))
            NCAprm <- est.nca(time=time,conc=conc,backExtrp=backExtrp,negConcExcl=negConcExcl,doseType=doseType,adminType=adminType,doseNm=doseNm,dose=dose,doseNumber=DoseNumber,doseAmt=idzAmt,method=method,AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,LambdaExclude=LambdaExclude,Tau=Tau,TI=TI,simFile=simFile,dset=dset) # calls est.nca function
            outData <- rbind(outData, cbind(data.frame(GROUP=igr,FLAG=iflag),data.frame(ID=idd[i],DoseNumber=DoseNumber,DoseAmount=idzAmt,t(NCAprm))))
          }
          
          if(noPlot=="FALSE"){
            plotData    <- subset(outData, DoseNumber==dose[d], select=c(AUClast,AUCINF_obs,Cmax,Tmax))
            if (nrow(plotData)<=5) next
            pltPrm      <- c("AUClast","AUCINF_obs","Cmax","Tmax")
            for (p in 1:length(pltPrm)){ if (nrow(plotData[plotData[,p] != "NaN",])<5) pltPrm <- pltPrm[-p]}
            if (length(pltPrm) == 0) next
            figlbl      <- paste(grNm,"-",as.character(grp[g]),"_",flNm,"-",as.character(flag[f]),"_",oidNm,"-",dose[d],sep="")
            histobsgrob <- histobs.plot(plotData=plotData,figlbl=figlbl,param=c("AUClast","AUCINF_obs","Cmax","Tmax"),cunit=cunit,tunit=tunit,spread=spread)
            gdr         <- histobsgrob$gdr
            mylegend    <- histobsgrob$legend
            lheight     <- histobsgrob$lheight
            if (printOut=="TRUE"){
              fl <- paste(usrdir,"/HistObs_",figlbl,sep="")
              eval(parse(text=paste(figFormat,"(file=\"",fl,".",figFormat,"\",height=15,width=14,units=\"cm\",res=200)",sep="")))
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
    cnm         <- c(paste("GID (",grNm,")",sep=""),paste("SGID (",flNm,")",sep=""),paste("OID (",oidNm,")",sep=""),paste("Dose (",dunit,")",sep=""),"No. of individuals")
    names(pddf) <- cnm
    
    if(noPlot=="FALSE"){
      # ggplot for conc vs. time
      fct <- unique(cdata$FCT); nfct <- length(fct)
      # ggplot for conc vs. time
      for (p in seq(1,nfct,3)){
        str <- ""
        if ((nfct-p)>=2) {
          dht   <- 16
          dwd   <- 15
          str   <- gsub("^_", "", paste(str,fct[p:(p+2)],sep="_",collapse=""))
          pdata <- subset(cdata, FCT == fct[p] | FCT == fct[p+1] | FCT == fct[p+2])
        }else if ((nfct-p)==1){
          dht   <- 11
          dwd   <- 13
          str   <- gsub("^_", "", paste(str,fct[p:(p+1)],sep="_",collapse=""))
          pdata <- subset(cdata, FCT == fct[p] | FCT == fct[p+1])
        }else if ((nfct-p)==0){
          dht   <- 7
          dwd   <- 11
          str   <- gsub("^_", "", paste(str,fct[p],sep="_",collapse=""))
          pdata <- subset(cdata, FCT == fct[p])
        }
        gdr <- dv.plot(pdata=pdata,cunit=cunit,tunit=tunit)
        suppressMessages(suppressWarnings(grid.arrange(gdr)))
        ggr <- grid.grab()
        concplot[[length(concplot)+1]] <- ggr
        if (printOut=="TRUE") ggsave(filename=paste(usrdir,"/TimeConc_",str,".",figFormat,sep=""),plot=gdr,height=dht,width=dwd,units="cm",dpi=200)
      }
    }
  }
  
  if (printOut=="TRUE" && is.null(simFile)){
    # write the output in a file
    tmpdf <- outData
    if (case == 1){
      names(outData)[2]         <- oidNm
      outData[,3:ncol(outData)] <- data.frame(lapply(outData[,3:ncol(outData)], function(x) round(as.numeric(x),digits=4)))
    }
    if (case == 2){
      names(outData)[c(1,3)]    <- c(grNm,oidNm)
      outData[,4:ncol(outData)] <- data.frame(lapply(outData[,4:ncol(outData)], function(x) round(as.numeric(x),digits=4)))
    }
    if (case == 3){
      names(outData)[c(1,3)]    <- c(flNm,oidNm)
      outData[,4:ncol(outData)] <- data.frame(lapply(outData[,4:ncol(outData)], function(x) round(as.numeric(x),digits=4)))
    }
    if (case == 4){
      names(outData)[c(1,2,4)]  <- c(grNm,flNm,oidNm)
      outData[,5:ncol(outData)] <- data.frame(lapply(outData[,5:ncol(outData)], function(x) round(as.numeric(x),digits=4)))
    }
    write.table(outData, file=paste(usrdir,"/ncaOutput.tsv",sep=""), sep="\t", row.names=F, col.names=T, quote=F)
    outData <- tmpdf; rm(tmpdf)
  }
  
  # Statistical analysis for each patient group
  stNm <- c("C0","Tmax","Cmax","Cmax_D","Tlast","Clast","AUClast","AUMClast","MRTlast","No_points_Lambda_z","AUC_pBack_Ext","AUClower_upper","Rsq","Rsq_adjusted","Corr_XY","Lambda_z","Lambda_z_lower","Lambda_z_upper","HL_Lambda_z","AUCINF_obs","AUCINF_D_obs","AUC_pExtrap_obs","Vz_obs","Cl_obs","AUCINF_pred","AUCINF_D_pred","AUC_pExtrap_pred","Vz_pred","Cl_pred","AUMCINF_obs","AUMC_pExtrap_obs","AUMCINF_pred","AUMC_pExtrap_pred","MRTINF_obs","MRTINF_pred","Vss_obs","Vss_pred","Tau","Tmin","Cmin","Cavg","p_Fluctuation","Accumulation_Index","Clss")
  grStat <- data.frame()
  if (case == 1){
    for (d in 1:ndose){
      pm <- data.frame(Ntot=numeric(0),Nunique=numeric(0),Min=numeric(0),Max=numeric(0),Mean=numeric(0),Median=numeric(0),SD=numeric(0),SE=numeric(0),CVp=numeric(0),CI95l=numeric(0),CI95u=numeric(0),gMean=numeric(0),gCVp=numeric(0))
      nm <- data.frame(character(0))
      counter <- 1
      for (i in 1:length(stNm)){
        Nm <- stNm[i]
        tdf <- as.numeric(as.character(outData[(outData$DoseNumber==dose[d] & outData[,Nm]!="NaN"),Nm]))
        if (length(tdf) < 2){
          nm <- rbind(nm, data.frame(Nm))
          pm[counter,] <- rep(NA,13)
        }else{
          nm <- rbind(nm, data.frame(Nm))
          stPrm <- calc.stat(x=tdf)        # calls calc.stat function
          stPrm <- unname(stPrm)
          pm[counter,] <- stPrm
        }
        counter <- counter + 1
      }
      if (nrow(pm) == 0) next
      pm <- data.frame(lapply(pm, function(x){if(is.numeric(x)){round(x,digits=4)}else{x}}))
      tmpStat <- t(cbind(nm, pm))
      rownames(tmpStat)[1] <- "Name"
      tmpStat <- cbind(DOSE=dose[d],Stat=rownames(tmpStat),tmpStat)
      for (cnum in 3:ncol(tmpStat)){colnames(tmpStat)[cnum] <- as.character(tmpStat[1,cnum])}
      tmpStat <- tmpStat[-1,]
      grStat <- rbind(grStat,tmpStat)
    }
    names(grStat)[1] <- oidNm
  }
  if (case == 2){
    for (g in 1:ngrp){
      for (d in 1:ndose){
        pm <- data.frame(Ntot=numeric(0),Nunique=numeric(0),Min=numeric(0),Max=numeric(0),Mean=numeric(0),Median=numeric(0),SD=numeric(0),SE=numeric(0),CVp=numeric(0),CI95l=numeric(0),CI95u=numeric(0),gMean=numeric(0),gCVp=numeric(0))
        nm <- data.frame(character(0))
        counter <- 1
        for (i in 1:length(stNm)){
          Nm <- stNm[i]
          tdf <- as.numeric(as.character(outData[(outData$DoseNumber==dose[d] & outData$GROUP==as.character(grp[g]) & outData[,Nm]!="NaN"),Nm]))
          if (length(tdf) < 2){
            nm <- rbind(nm, data.frame(Nm))
            pm[counter,] <- rep(NA,13)
          }else{
            nm <- rbind(nm, data.frame(Nm))
            stPrm <- calc.stat(x=tdf)        # calls calc.stat function
            stPrm <- unname(stPrm)
            pm[counter,] <- stPrm
          }
          counter <- counter + 1
        }
        if (nrow(pm) == 0) next
        pm <- data.frame(lapply(pm, function(x){if(is.numeric(x)){round(x,digits=4)}else{x}}))
        tmpStat <- t(cbind(nm, pm))
        rownames(tmpStat)[1] <- "Name"
        tmpStat <- cbind(GRP=as.character(grp[g]),DOSE=dose[d],Stat=rownames(tmpStat),tmpStat)
        for (cnum in 4:ncol(tmpStat)){colnames(tmpStat)[cnum] <- as.character(tmpStat[1,cnum])}
        tmpStat <- tmpStat[-1,]
        grStat <- rbind(grStat,tmpStat)
      }
    }
    names(grStat)[c(1,2)] <- c(grNm,oidNm)
  }
  if (case == 3){
    for (f in 1:nflag){
      for (d in 1:ndose){
        pm <- data.frame(Ntot=numeric(0),Nunique=numeric(0),Min=numeric(0),Max=numeric(0),Mean=numeric(0),Median=numeric(0),SD=numeric(0),SE=numeric(0),CVp=numeric(0),CI95l=numeric(0),CI95u=numeric(0),gMean=numeric(0),gCVp=numeric(0))
        nm <- data.frame(character(0))
        counter <- 1
        for (i in 1:length(stNm)){
          Nm <- stNm[i]
          tdf <- as.numeric(as.character(outData[(outData$DoseNumber==dose[d] & outData$FLAG==as.character(flag[f]) & outData[,Nm]!="NaN"),Nm]))
          if (length(tdf) < 2){
            nm <- rbind(nm, data.frame(Nm))
            pm[counter,] <- rep(NA,13)
          }else{
            nm <- rbind(nm, data.frame(Nm))
            stPrm <- calc.stat(x=tdf)        # calls calc.stat function
            stPrm <- unname(stPrm)
            pm[counter,] <- stPrm
          }
          counter <- counter + 1
        }
        if (nrow(pm) == 0) next
        pm <- data.frame(lapply(pm, function(x){if(is.numeric(x)){round(x,digits=4)}else{x}}))
        tmpStat <- t(cbind(nm, pm))
        rownames(tmpStat)[1] <- "Name"
        tmpStat <- cbind(GRP=as.character(flag[f]),DOSE=dose[d],Stat=rownames(tmpStat),tmpStat)
        for (cnum in 4:ncol(tmpStat)){colnames(tmpStat)[cnum] <- as.character(tmpStat[1,cnum])}
        tmpStat <- tmpStat[-1,]
        grStat <- rbind(grStat,tmpStat)
      }
    }
    names(grStat)[c(1,2)] <- c(flNm,oidNm)
  }
  if (case == 4){
    for (g in 1:ngrp){
      for (f in 1:nflag){
        for (d in 1:ndose){
          pm <- data.frame(Ntot=numeric(0),Nunique=numeric(0),Min=numeric(0),Max=numeric(0),Mean=numeric(0),Median=numeric(0),SD=numeric(0),SE=numeric(0),CVp=numeric(0),CI95l=numeric(0),CI95u=numeric(0),gMean=numeric(0),gCVp=numeric(0))
          nm <- data.frame(character(0))
          counter <- 1
          for (i in 1:length(stNm)){
            Nm <- stNm[i]
            tdf <- as.numeric(as.character(outData[(outData$DoseNumber==dose[d] & outData$GROUP==as.character(grp[g]) & outData$FLAG==as.character(flag[f]) & outData[,Nm]!="NaN"),Nm]))
            if (length(tdf) < 2){
              nm <- rbind(nm, data.frame(Nm))
              pm[counter,] <- rep(NA,13)
            }else{
              nm <- rbind(nm, data.frame(Nm))
              stPrm <- calc.stat(x=tdf)        # calls calc.stat function
              stPrm <- unname(stPrm)
              pm[counter,] <- stPrm
            }
            counter <- counter + 1
          }
          if (nrow(pm) == 0) next
          pm <- data.frame(lapply(pm, function(x){if(is.numeric(x)){round(x,digits=4)}else{x}}))
          tmpStat <- t(cbind(nm, pm))
          rownames(tmpStat)[1] <- "Name"
          tmpStat <- cbind(GRP=as.character(grp[g]),FLAG=as.character(flag[f]),DOSE=dose[d],Stat=rownames(tmpStat),tmpStat)
          for (cnum in 5:ncol(tmpStat)){colnames(tmpStat)[cnum] <- as.character(tmpStat[1,cnum])}
          tmpStat <- tmpStat[-1,]
          grStat <- rbind(grStat,tmpStat)
        }
      }
    }
    names(grStat)[c(1:3)] <- c(grNm,flNm,oidNm)
  }
  if(printOut=="TRUE") write.table(grStat, file=paste(usrdir,"/ObsStat.tsv",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
  
  ############################################################################
  # Analyze the simulated data if exists
  if (is.null(simFile)){
    streamsEnv <- parent.frame()
    if(exists("outData")) assign("ncaOutput",  outData,   envir=streamsEnv)
    if(exists("grStat"))  assign("ObsStat",    grStat,    envir=streamsEnv)
    
    if(printOut=="TRUE"){
      if (case == 1){
        prnTab <- head(cbind(outData[,1:3], subset(outData, select = tabCol)), 100)
        names(prnTab)[2] <- oidNm
        prnTab[,3:ncol(prnTab)] <- data.frame(lapply(prnTab[,3:ncol(prnTab)], function(x) round(as.numeric(x),digits=4)))
      }else if (case == 2){
        prnTab <- head(cbind(outData[,1:4], subset(outData, select = tabCol)),100)
        names(prnTab)[c(1,3)] <- c(grNm,oidNm)
        prnTab[,4:ncol(prnTab)] <- data.frame(lapply(prnTab[,4:ncol(prnTab)], function(x) round(as.numeric(x),digits=4)))
      }else if (case == 3){
        prnTab <- head(cbind(outData[,1:4], subset(outData, select = tabCol)),100)
        names(prnTab)[c(1,3)] <- c(flNm,oidNm)
        prnTab[,4:ncol(prnTab)] <- data.frame(lapply(prnTab[,4:ncol(prnTab)], function(x) round(as.numeric(x),digits=4)))
      }else if (case == 4){
        prnTab <- head(cbind(outData[,1:5], subset(outData, select = tabCol)),100)
        names(prnTab)[c(1,2,4)] <- c(grNm,flNm,oidNm)
        prnTab[,5:ncol(prnTab)] <- data.frame(lapply(prnTab[,5:ncol(prnTab)], function(x) round(as.numeric(x),digits=4)))
      }
      prnTab <- data.frame(lapply(prnTab, function(x){if(is.numeric(x)){round(x,digits=4)}else{x}}))
      fnOut <- list(arglist=match.call(), TXT=txt, pddf=pddf, prnTab=prnTab, spread=spread, conc=concplot, histobs=histobsplot)
    }
  }else{
    od <- paste(usrdir,"/SIMDATA",sep="")
    if (file.exists(od)){
      dirTest <- readline("Directory \"SIMDATA\" already exists.\n
                          Overwrite it? (type 1)\n
                          Rename the existing folder and create a new one? (type 2)\n
                          Use data from the existing folder? (type 3)\n")
      if (dirTest == "1"){
        unlink(od, recursive=T)
        dir.create(od)
      }else if (dirTest == "2"){
        print("\nRenaming \"SIMDATA\" to \"SIMDATA_PREVIOUS\"\n")
        file.rename(from="SIMDATA", to="SIMDATA_PREVIOUS")
        dir.create(od)
      }else if (dirTest == "3"){
        file_list <- list.files(path="./SIMDATA/", pattern="sim_[0-9]*.csv", full.names=T)
      }else{setwd(usrdir);stop("Bad choice!!!\n")}
    }else{
      dir.create(od)
      dirTest <- "0"
    }
    
    if (dirTest != 3){
      # read NONMEM output into individual simulation data file
      IPSIM <- function(table.sim,MDV.rm){
        if(missing(MDV.rm)){MDV.rm=T}
        sim <- read.table(table.sim,skip=1,header=T,fill=T,as.is=T)
        sim <- as.data.frame(apply(sim,2,as.numeric))
        Nro <- min(which(is.na(sim[,1])))-1
        sim <- sim[!is.na(sim[,1]),]
        sim$NSUB <- rep(1:(nrow(sim)/Nro),each=Nro)
        if(MDV.rm==T){
          if(any(colnames(sim)=='MDV')){sim <- sim[sim[,'MDV']==0,]
          }else{cat('\nWarning MDV data item not listed in header,
                         Could not remove dose events!')}
        }
        assign("nmdf", sim)
        return(nmdf)
      }
      nmdf <- IPSIM(simFile,MDV.rm=F); simID <- unique(nmdf$NSUB); nsim <- length(unique(nmdf$NSUB))
      if (printOut=="TRUE") write.table(nmdf, file=paste(usrdir,"/ncaSimData.tsv",sep=""), row.names=F, quote=F, sep="\t")
      
      srdf <- nmdf[nmdf$NSUB == 1,]  # copy simulated data before processing

      if (idNmSim%in%colnames(nmdf)==F | timeNmSim%in%colnames(nmdf)==F | concNmSim%in%colnames(nmdf)==F){
        setwd(usrdir);stop("Incorrect column names of ID, TIME and/or DV in simulation output\n")
      }else{
        idCol   <- which(colnames(nmdf) == idNmSim)
        timeCol <- which(colnames(nmdf) == timeNmSim)
        concCol <- which(colnames(nmdf) == concNmSim)
      }
      
      if (!is.null(grNm)){
        if (grNm%in%colnames(nmdf)==F){setwd(usrdir);stop("Incorrect name for the group column in simulation output\n")}else{grCol <- which(colnames(nmdf) == grNm)}
        for (i in 1:ngrp){
          if (length(nmdf[nmdf[,grCol]==grp[i],idCol]) == 0){setwd(usrdir);stop("Group identifier does not match the group column in simulation output\n")}
        }
      }
      if (!is.null(flNm)){
        if (flNm%in%colnames(nmdf)==F){setwd(usrdir);stop("Incorrect name for the flag column\n")}else{flCol <- which(colnames(nmdf) == flNm)}
        if (is.null(flag)){flag <- unique(nmdf[,flCol])}
        nflag <- length(flag)
        for (i in 1:nflag){
          if (length(nmdf[nmdf[,flCol]==flag[i],idCol]) == 0){setwd(usrdir);stop("Flag identifier does not match the flag column\n")}
        }
      }
      
      # Dose identifiers. For missing doseNm argument, data is assumed to have single dose.
      # Dose amount is extracted from doseAmtNm column
      if (!is.null(doseNm)){
        if (doseNm%in%colnames(nmdf)==F){setwd(usrdir);stop("Incorrect name for the dose column\n")}else{doseCol <- which(colnames(nmdf) == doseNm)}
        oidNm <- doseNm
        if (is.null(dose)){dose <- unique(sort(nmdf[,doseCol]))}
        ndose <- length(dose)
        for (i in 1:ndose){
          if (nrow(nmdf[nmdf[,doseCol]==dose[i],]) == 0){setwd(usrdir);stop("Dose identifier does not match the dose column\n")}
        }
      }else{
        dose <- 1; ndose <- 1
      }
      
      if (is.null(doseAmtNm)){
        if ("AMT"%in%colnames(nmdf) == F){setwd(usrdir);stop("Dose amount column is required as AMT column is absent\n")}
        doseAmtNm <- "AMT"
      }
      
      # ignore data with BLQ = 1 or user specified value (optional)
      if (!is.null(blqNm)){
        if (blqNm%in%colnames(nmdf) == T){
          blqCol <- which(colnames(nmdf) == blqNm)
          for (i in 1:length(blqExcl)) {nmdf <- nmdf[nmdf[,blqNm] != blqExcl[i],]}
        }else{setwd(usrdir);stop("Incorrect BLQ column name in simulation output\n")}
      }
      
      # include data based on specific values on EVID column (optional) but keep rows with TIME == 0
      if (evid == "TRUE"){
        if ("EVID"%in%colnames(nmdf) == T){
          # uevid == unique values in EVID column
          # evidIncl == EVID values to be included
          # ievid == EVID values to be ignored
          uevid <- unique(as.numeric(as.character(nmdf$EVID))); ievid <- setdiff(uevid, as.numeric(evidIncl))
          if (length(ievid) != 0){
            for (i in 1:length(ievid)){
              if (ievid[i] != 1){
                nmdf <- nmdf[nmdf$EVID != ievid[i],]
              }else{
                if (length(which(as.numeric(as.character(nmdf[,timeCol])) != 0 & as.numeric(as.character(nmdf$EVID)) == as.numeric(ievid[i]))) == 0) next
                nmdf <- nmdf[-which(as.numeric(as.character(nmdf[,timeCol])) != 0 & as.numeric(as.character(nmdf$EVID)) == as.numeric(ievid[i])),]
              }
            }
          }
        }else{setwd(usrdir);stop("Incorrect EVID column name in simulation output\n")}
      }
      
      # if MDV fiter is present, exclude data for MDV == 1 but keep rows with TIME == 0
      if (mdv == "TRUE"){
        if ("MDV"%in%colnames(nmdf) == T){
          if (length(which(as.numeric(as.character(nmdf[,timeCol])) != 0 & as.numeric(as.character(nmdf$MDV)) == 1)) == 0) next
          nmdf <- nmdf[-which(as.numeric(as.character(nmdf[,timeCol])) != 0 & as.numeric(as.character(nmdf$MDV)) == 1),]
        }else{setwd(usrdir);stop("Incorrect MDV column name in simulation output\n")}
      }
      
      # exclude data based on specific values on filter column (optional)
      if (!is.null(filterNm)){
        if (filterNm%in%colnames(nmdf)==T & !is.null(filterExcl)){
          # filterExcl  == values to be excluded
          filterCol <- which(colnames(nmdf) == filterNm)
          for (i in 1:length(filterExcl)){
            if (grepl("^[-]?[0-9]*[.]?[0-9]*[eE]?[-]?[0-9]*[.]?[0-9]*$", filterExcl[i]) == T) {nmdf <- nmdf[nmdf[,filterCol] != as.numeric(filterExcl[i]),]}
            if (grepl("^[-]?[0-9]*[.]?[0-9]*[eE]?[-]?[0-9]*[.]?[0-9]*$", filterExcl[i]) == F) {nmdf <- eval(parse(text=paste("subset(nmdf, !",filterNm,"%in% nmdf[nmdf[,",filterCol,"]",filterExcl[i],",filterCol])",sep="")))}
          }
        }else{setwd(usrdir);stop("Incorrect filterNm or filterExcl specification in simulation output\n")}
      }
      
      # Calculate AUC parameters for the simulation output
      dset = "sim"
      # Function to extract time and conc data for NCA metrics calculation for simulated data
      simNcaId <- function(ifdf,ID){
        if(adminType == "iv-infusion" & is.null(TI)){
          amt  <- ifdf[ifdf[,idCol]==ID & ifdf$AMT > 0,"AMT"][1]
          rate <- ifdf[ifdf[,idCol]==ID & ifdf$RATE > 0,"RATE"][1]
          if (is.na(amt) | is.na(rate) | rate==0){setwd(usrdir);stop(paste("Incorrect AMT and/or RATE value IN NONMEM output for ",ID,sep=""))}else{TI <- amt/rate}
        }else{TI <- "NaN"}
        conc <- as.numeric(as.character(ifdf[ifdf[,idCol]==ID,concCol]))
        if (simLog == "TRUE"){
          for (c in 1:length(conc)){conc[c] <- ifelse ((conc[c] != 0), exp(conc[c]), 0)}
        }
        time   <- as.numeric(as.character(ifdf[ifdf[,idCol]==ID,timeCol]))
        
        tc <- data.frame(time,conc)
        if(nrow(tc)>0){
          tc <- tc[order(tc$time),]
          tc$time <- tc$time + abs(min(tc$time))
        }
        return(tc)
      }
      
      for (s in 1:nsim){
        simData <- data.frame()
        smdf    <- nmdf[nmdf$NSUB == simID[s],]
        # Calculate NCA parameters
        if (case == 1){
          for (d in 1:ndose){
            if (!is.null(doseNm)){ifdf <- smdf[smdf[,doseNm]==dose[d],]}else{ifdf <- smdf}
            if (length(which(is.na(ifdf[,concCol]) | ifdf[,concCol]=="")) != 0){ifdf <- ifdf[-which(is.na(ifdf[,concCol]) | ifdf[,concCol]==""),]}
            if (nrow(ifdf) == 0){next}
            idd        <- unique(ifdf[,idCol])
            DoseNumber <- dose[d]
            for (i in 1:length(idd)){
              if (!is.null(doseAmtNm) & is.null(doseNm)){
                idzAmt <- as.numeric(srdf[srdf[,idNmSim]==idd[i] & srdf[,doseAmtNm] > 0, doseAmtNm][1])
              }else if(!is.null(doseAmtNm) & !is.null(doseNm)){
                idzAmt <- as.numeric(srdf[srdf[,doseNm]==dose[d] & srdf[,idNmSim]==idd[i] & srdf[,doseAmtNm] > 0, doseAmtNm][1])
              }else{
                idzAmt <- NA
              }
              if(nrow(ifdf[ifdf$ID==idd[i],])==0) next
              stc  <- simNcaId(ifdf,idd[i])
              time <- as.numeric(stc$time)
              conc <- as.numeric(stc$conc)
              NCAprm <- est.nca(time=time,conc=conc,backExtrp=backExtrp,negConcExcl=negConcExcl,doseType=doseType,adminType=adminType,doseNm=doseNm,dose=dose,doseNumber=DoseNumber,doseAmt=idzAmt,method=method,AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,LambdaExclude=LambdaExclude,Tau=Tau,TI=TI,simFile=simFile,dset=dset) # calls est.nca function
              simData <- rbind(simData, data.frame(ID=idd[i],DoseNumber=dose[d],DoseAmount=idzAmt,t(NCAprm),NSIM=s))
            }
          }
        }
        if (case == 2){
          for (g in 1:ngrp){
            for (d in 1:ndose){
              if (!is.null(doseNm)){ifdf <- smdf[smdf[,grCol]==as.character(grp[g]) & smdf[,doseNm]==dose[d],]}else{ifdf <- smdf[smdf[,grCol]==as.character(grp[g]),]}
              if (length(which(is.na(ifdf[,concCol]) | ifdf[,concCol]=="")) != 0){ifdf <- ifdf[-which(is.na(ifdf[,concCol]) | ifdf[,concCol]==""),]}
              if (nrow(ifdf) == 0){next}
              idd        <- unique(ifdf[,idCol])
              DoseNumber <- dose[d]
              igr        <- as.character(grp[g])
              for (i in 1:length(idd)){
                if (!is.null(doseAmtNm) & is.null(doseNm)){
                  idzAmt <- as.numeric(srdf[srdf[,grCol]==as.character(grp[g]) & srdf[,idNmSim]==idd[i] & srdf[,doseAmtNm] > 0, doseAmtNm][1])
                }else if(!is.null(doseAmtNm) & !is.null(doseNm)){
                  idzAmt <- as.numeric(srdf[srdf[,grCol]==as.character(grp[g]) & srdf[,doseNm]==dose[d] & srdf[,idNmSim]==idd[i] & srdf[,doseAmtNm] > 0, doseAmtNm][1])
                }else{
                  idzAmt <- NA
                }
                if(nrow(ifdf[ifdf$ID==idd[i],])==0) next
                stc  <- simNcaId(ifdf,idd[i])
                time <- as.numeric(stc$time)
                conc <- as.numeric(stc$conc)
                NCAprm <- est.nca(time=time,conc=conc,backExtrp=backExtrp,negConcExcl=negConcExcl,doseType=doseType,adminType=adminType,doseNm=doseNm,dose=dose,doseNumber=DoseNumber,doseAmt=idzAmt,method=method,AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,LambdaExclude=LambdaExclude,Tau=Tau,TI=TI,simFile=simFile,dset=dset) # calls est.nca function
                simData <- rbind(simData, cbind(data.frame(GROUP=igr),data.frame(ID=idd[i],DoseNumber=DoseNumber,DoseAmount=idzAmt,t(NCAprm),NSIM=s)))
              }
            }
          }
        }
        if (case == 3){
          for (f in 1:nflag){
            for (d in 1:ndose){
              if (!is.null(doseNm)){ifdf <- smdf[smdf[,flCol]==as.character(flag[f]) & smdf[,doseNm]==dose[d],]}else{ifdf <- smdf[smdf[,flCol]==as.character(flag[f]),]}
              if (length(which(is.na(ifdf[,concCol]) | ifdf[,concCol]=="")) != 0){ifdf <- ifdf[-which(is.na(ifdf[,concCol]) | ifdf[,concCol]==""),]}
              if (nrow(ifdf) == 0){next}
              idd        <- unique(ifdf[,idCol])
              DoseNumber <- dose[d]
              iflag      <- as.character(flag[f])
              if (!is.null(doseAmtNm) & is.null(doseNm)){
                doseAmount <- paste(as.numeric(unique(srdf[srdf[,flCol]==as.character(flag[f]) & srdf[,doseAmtNm] > 0, doseAmtNm])), collapse=", ")
              }else if(!is.null(doseAmtNm) & !is.null(doseNm)){
                doseAmount <- paste(as.numeric(unique(srdf[srdf[,flCol]==as.character(flag[f]) & srdf[,doseNm]==dose[d] & srdf[,doseAmtNm] > 0, doseAmtNm])), collapse=", ")
              }else{
                doseAmount <- NA
              }
              for (i in 1:length(idd)){
                if (!is.null(doseAmtNm) & is.null(doseNm)){
                  idzAmt <- as.numeric(srdf[srdf[,flCol]==as.character(flag[f]) & srdf[,idNmSim]==idd[i] & srdf[,doseAmtNm] > 0, doseAmtNm][1])
                }else if(!is.null(doseAmtNm) & !is.null(doseNm)){
                  idzAmt <- as.numeric(srdf[srdf[,flCol]==as.character(flag[f]) & srdf[,doseNm]==dose[d] & srdf[,idNmSim]==idd[i] & srdf[,doseAmtNm] > 0, doseAmtNm][1])
                }else{
                  idzAmt <- NA
                }
                if(nrow(ifdf[ifdf$ID==idd[i],])==0) next
                stc  <- simNcaId(ifdf,idd[i])
                time <- as.numeric(stc$time)
                conc <- as.numeric(stc$conc)
                NCAprm <- est.nca(time=time,conc=conc,backExtrp=backExtrp,negConcExcl=negConcExcl,doseType=doseType,adminType=adminType,doseNm=doseNm,dose=dose,doseNumber=DoseNumber,doseAmt=idzAmt,method=method,AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,LambdaExclude=LambdaExclude,Tau=Tau,TI=TI,simFile=simFile,dset=dset) # calls est.nca function
                simData <- rbind(simData, cbind(data.frame(FLAG=iflag),data.frame(ID=idd[i],DoseNumber=DoseNumber,DoseAmount=idzAmt,t(NCAprm),NSIM=s)))
              }
            }
          }
        }
        if (case == 4){
          for (f in 1:nflag){
            for (g in 1:ngrp){
              for (d in 1:ndose){
                if (!is.null(doseNm)){ifdf <- smdf[(smdf[,flCol]==as.character(flag[f]) & smdf[,grCol]==as.character(grp[g]) & smdf[,doseNm]==dose[d]),]}else{ifdf <- smdf[smdf[,flCol]==as.character(flag[f]) & smdf[,grCol]==as.character(grp[g]),]}
                if (length(which(is.na(ifdf[,concCol]) | ifdf[,concCol]=="")) != 0){ifdf <- ifdf[-which(is.na(ifdf[,concCol]) | ifdf[,concCol]==""),]}
                if (nrow(ifdf) == 0){next}
                idd        <- unique(ifdf[,idCol])
                DoseNumber <- dose[d]
                igr        <- as.character(grp[g])
                iflag      <- as.character(flag[f])
                for (i in 1:length(idd)){
                  if (!is.null(doseAmtNm) & is.null(doseNm)){
                    idzAmt <- as.numeric(srdf[srdf[,grCol]==as.character(grp[g]) & srdf[,flCol]==as.character(flag[f]) & srdf[,idNmSim]==idd[i] & srdf[,doseAmtNm] > 0, doseAmtNm][1])
                  }else if(!is.null(doseAmtNm) & !is.null(doseNm)){
                    idzAmt <- as.numeric(srdf[srdf[,grCol]==as.character(grp[g]) & srdf[,flCol]==as.character(flag[f]) & srdf[,doseNm]==dose[d] & srdf[,idNmSim]==idd[i] & srdf[,doseAmtNm] > 0, doseAmtNm][1])
                  }else{
                    idzAmt <- NA
                  }
                  if(nrow(ifdf[ifdf$ID==idd[i],])==0) next
                  stc  <- simNcaId(ifdf,idd[i])
                  time <- as.numeric(stc$time)
                  conc <- as.numeric(stc$conc)
                  NCAprm <- est.nca(time=time,conc=conc,backExtrp=backExtrp,negConcExcl=negConcExcl,doseType=doseType,adminType=adminType,doseNm=doseNm,dose=dose,doseNumber=DoseNumber,doseAmt=idzAmt,method=method,AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,LambdaExclude=LambdaExclude,Tau=Tau,TI=TI,simFile=simFile,dset=dset) # calls est.nca function
                  simData <- rbind(simData, cbind(data.frame(GROUP=igr,FLAG=iflag),data.frame(ID=idd[i],DoseNumber=DoseNumber,DoseAmount=idzAmt,t(NCAprm),NSIM=s)))
                }
              }
            }
          }
        }
        if (printOut=="TRUE") write.csv(simData, file=paste(od,"/sim_",s,".csv",sep=""), row.names=F, quote=F)
      }
    }
    
    # Statistical analysis for each individual
    setwd(od)
    
    # read all simulated NCA parameters to a list
    lasdf <- lapply(list.files(pattern="sim_[0-9]*.csv",full.names=T),function(i){read.csv(i, header=T)})
    file_list <- list.files(path=".", pattern="sim_[0-9]*.csv", full.names=T)
    nsim <- length(file_list)
    dasdf <- do.call(rbind, lapply(file_list, read.csv))
    if (printOut=="TRUE") write.table(dasdf, file=paste(usrdir,"/ncaSimEst.tsv",sep=""), row.names=F, quote=F, sep="\t")
    
    # Population histogram
    if (case == 1){
      for (d in 1:ndose){
        if (nrow(dasdf[dasdf$DoseNumber==dose[d],]) == 0) next
        smeanData <- data.frame()
        for (i in 1:length(lasdf)){
          tmdf   <- subset(data.frame(lasdf[[i]]), select=param, DoseNumber==dose[d])
          tmpPrm <- as.data.frame(lapply(tmdf, FUN=function(x) mean(as.numeric(x[!is.na(x)]))))
          smeanData <- rbind(smeanData, tmpPrm)
        }
        if(noPlot=="FALSE"){
          obsdata     <- subset(outData, select=param, ID!="" & DoseNumber==dose[d])
          figlbl      <- paste(oidNm,"-",dose[d],sep="")
          histpopgrob <- histpop.plot(obsdata=obsdata,simdata=smeanData,figlbl=figlbl,param=param,cunit=cunit,tunit=tunit,spread=spread)
          gdr         <- histpopgrob$gdr
          mylegend    <- histpopgrob$legend
          lheight     <- histpopgrob$lheight
          if(printOut=="TRUE"){
            fl <- paste(usrdir,"/popMean_",figlbl,sep="")
            eval(parse(text=paste(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",res=200)",sep="")))
            suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
            dev.off()
          }
          suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
          ggr <- grid.grab()
          popplot[[length(popplot)+1]] <- ggr
        }
      }
    }
    
    # Group based histogram for no multiple flag
    if (case == 2){
      for (g in 1:ngrp){
        for (d in 1:ndose){
          if (nrow(dasdf[dasdf$DoseNumber==dose[d] & dasdf$GROUP==as.character(grp[g]),]) == 0) next
          smeanData <- data.frame()
          for (i in 1:length(lasdf)){
            tmdf   <- subset(data.frame(lasdf[[i]]), select=param, GROUP==as.character(grp[g]) & DoseNumber==dose[d])
            tmpPrm <- as.data.frame(lapply(tmdf, FUN=function(x) mean(as.numeric(x[!is.na(x)]))))
            smeanData <- rbind(smeanData, tmpPrm)
          }
          if(noPlot=="FALSE"){
            obsdata     <- subset(outData, select=param, ID!="" & DoseNumber==dose[d] & GROUP==as.character(grp[g]))
            figlbl      <- paste(grNm,"-",as.character(grp[g]),"_",oidNm,"-",dose[d],sep="")
            histpopgrob <- histpop.plot(obsdata=obsdata,simdata=smeanData,figlbl=figlbl,param=param,cunit=cunit,tunit=tunit,spread=spread)
            gdr         <- histpopgrob$gdr
            mylegend    <- histpopgrob$legend
            lheight     <- histpopgrob$lheight
            if (printOut=="TRUE"){
              fl <- paste(usrdir,"/popMean_",figlbl,sep="")
              eval(parse(text=paste(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",res=200)",sep="")))
              suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
              dev.off()
            }
            suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
            ggr <- grid.grab()
            popplot[[length(popplot)+1]] <- ggr
          }
        }
      }
    }
    
    # Flag based histogram
    if (case == 3){
      for (f in 1:nflag){
        for (d in 1:ndose){
          if (nrow(dasdf[dasdf$DoseNumber==dose[d] & dasdf$FLAG==as.character(flag[f]),]) == 0) next
          smeanData <- data.frame()
          for (i in 1:length(lasdf)){
            tmdf   <- subset(data.frame(lasdf[[i]]), select=param, FLAG==as.character(flag[f]) & DoseNumber==dose[d])
            tmpPrm <- as.data.frame(lapply(tmdf, FUN=function(x) mean(as.numeric(x[!is.na(x)]))))
            smeanData <- rbind(smeanData, tmpPrm)
          }
          if(noPlot=="FALSE"){
            obsdata     <- subset(outData, select=param, ID!="" & DoseNumber==dose[d] & FLAG==as.character(flag[f]))
            figlbl      <- paste(flNm,"-",as.character(flag[f]),"_",oidNm,"-",dose[d],sep="")
            histpopgrob <- histpop.plot(obsdata=obsdata,simdata=smeanData,figlbl=figlbl,param=param,cunit=cunit,tunit=tunit,spread=spread)
            gdr         <- histpopgrob$gdr
            mylegend    <- histpopgrob$legend
            lheight     <- histpopgrob$lheight
            if (printOut=="TRUE"){
              fl <- paste(usrdir,"/popMean_",figlbl,sep="")
              eval(parse(text=paste(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",res=200)",sep="")))
              suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
              dev.off()
            }
            suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
            ggr <- grid.grab()
            popplot[[length(popplot)+1]] <- ggr
          }
        }
      }
    }
    
    # Group and flag based histogram
    if (case == 4){
      for (g in 1:ngrp){
        for (f in 1:nflag){
          for (d in 1:ndose){
            if (nrow(dasdf[dasdf$DoseNumber==dose[d] & dasdf$GROUP==as.character(grp[g]) & dasdf$FLAG==as.character(flag[f]),]) == 0) next
            smeanData <- data.frame()
            for (i in 1:length(lasdf)){
              tmdf   <- subset(data.frame(lasdf[[i]]), select=param, GROUP==as.character(grp[g]) & FLAG==as.character(flag[f]) & DoseNumber==dose[d])
              tmpPrm <- as.data.frame(lapply(tmdf, FUN=function(x) mean(as.numeric(x[!is.na(x)]))))
              smeanData <- rbind(smeanData, tmpPrm)
            }
            if(noPlot=="FALSE"){
              obsdata     <- subset(outData, select=param, ID!="" & DoseNumber==dose[d] & GROUP==as.character(grp[g]) & FLAG==as.character(flag[f]))
              figlbl      <- paste(grNm,"-",as.character(grp[g]),"_",flNm,"-",as.character(flag[f]),"_",oidNm,"-",dose[d],sep="")
              histpopgrob <- histpop.plot(obsdata=obsdata,simdata=smeanData,figlbl=figlbl,param=param,cunit=cunit,tunit=tunit,spread=spread)
              gdr         <- histpopgrob$gdr
              mylegend    <- histpopgrob$legend
              lheight     <- histpopgrob$lheight
              if (printOut=="TRUE"){
                fl <- paste(usrdir,"/popMean_",figlbl,sep="")
                eval(parse(text=paste(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",res=200)",sep="")))
                suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
                dev.off()
              }
              suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
              ggr <- grid.grab()
              popplot[[length(popplot)+1]] <- ggr
            }
          }
        }
      }
    }
    
    devcol  <- paste("d",param,sep="")
    npdecol <- paste("npde",param,sep="")
    # ggplot options for the forest plot
    ggOpt_forest <- list(scale_color_manual(name="",values=c("mean"="red","SD"="darkgreen")),
                         theme(plot.title = element_text(size=10,face="bold"),
                               axis.title.x = element_text(size=9,face="bold"),
                               axis.title.y = element_text(size=9,face="bold"),
                               axis.text.x  = element_text(size=9,face="bold",color="black",angle=45,vjust=1,hjust=1),
                               axis.text.y  = element_text(size=9,face="bold",color="black",hjust=0),
                               legend.text  = element_text(size=9,face="bold"),
                               legend.background = element_rect(),
                               legend.position = "bottom", legend.direction = "horizontal",
                               legend.key.size = unit(0.8, "cm"),
                               panel.margin = unit(0.5, "cm"),
                               plot.margin  = unit(c(0.5,0.5,0.5,0.5), "cm")),
                         facet_wrap(~type, scales="free", ncol=2),
                         theme(strip.text.x = element_text(size=9, face="bold")))
    
    OTL   <- data.frame(No_of_outliers=numeric(0),ID_metric=character(0))
    npde  <- data.frame()
    fpval <- data.frame(type=character(0),mean=numeric(0),mcil=numeric(0),mciu=numeric(0),sdu=numeric(0),sducil=numeric(0),sduciu=numeric(0),str=character(0))
    if (case == 1){
      for (d in 1:ndose){
        if (nrow(dasdf[dasdf$DoseNumber==dose[d],]) == 0) next
        tdasdf <- subset(dasdf, DoseNumber==dose[d])
        id     <- unique(tdasdf$ID)
        pde    <- data.frame()
        metric <- ""
        nout   <- 0
        for (i in 1:length(id)){
          obsdata <- subset(outData, ID==id[i] & DoseNumber==dose[d])
          simdata <- subset(tdasdf, ID==id[i])
          figlbl  <- paste(oidNm,"-",dose[d],sep="")
          pdeout  <- nca.pde.deviation.outlier(obsdata=obsdata,simdata=simdata,idNm="ID",id=id[i],spread=spread,figlbl=figlbl,calcparam=alwprm,diagparam=param,cunit=cunit,tunit=tunit,noPlot=noPlot)
          pde     <- rbind(pde, cbind(data.frame(ID=id[i],DoseNumber=dose[d]), pdeout$pde))
          outData[(outData$ID==id[i] & outData$DoseNumber==dose[d]),] <- pdeout$obsdata
          if (pdeout$metric != ""){
            nout     <- nout + 1
            metric   <- paste(metric,pdeout$metric,sep=", ")
            if(noPlot=="FALSE"){
              gdr      <- pdeout$grob
              mylegend <- pdeout$legend
              lheight  <- pdeout$lheight
              if (printOut=="TRUE"){
                fl <- paste(usrdir,"/Outlier_ID-",id[i],"_",figlbl,sep="")
                eval(parse(text=paste(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",res=200)",sep="")))
                suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
                dev.off()
              }
              suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
              ggr <- grid.grab()
              outlierplot[[length(outlierplot)+1]] <- ggr
            }
          }
        }
        if (metric != "") metric <- gsub("^, ", "", metric)
        OTL  <- rbind(OTL, data.frame(No_of_outlier=nout,ID_metric=metric))
        npde <- rbind(npde,pde)
      }
      npde   <- nca.npde(pdedata=npde,pdecol=alwprm)
      npdeNm <- paste("npde",alwprm,sep="")
      for (r in 1:nrow(outData)){
        if (nrow(npde[(npde$ID==outData$ID[r] & npde$DoseNumber==outData$DoseNumber[r]),])!=1){
          outData[r,npdeNm] <- "NaN"
        }else{
          outData[r,npdeNm] <- npde[(npde$ID==outData$ID[r] & npde$DoseNumber==outData$DoseNumber[r]),npdeNm]
        }
      }
      tmpdf <- outData
      names(outData)[2] <- oidNm
      outData <- as.data.frame(lapply(outData, FUN=function(x) round(as.numeric(x), digits=4)))
      if (printOut=="TRUE") write.table(outData, file=paste(usrdir,"/ncaOutput.tsv",sep=""), sep="\t", row.names=F, col.names=T, quote=F)
      outData <- tmpdf; rm(tmpdf)
      
      if(noPlot=="FALSE"){
        for (d in 1:ndose){
          plotdata <- subset(outData, DoseNumber==dose[d])
          if (nrow(plotdata) == 0) next
          figlbl <- paste(oidNm,"-",dose[d])
          # Deviation plot
          ggdev <- nca.deviation.plot(plotdata=plotdata,xvar="ID",devcol=devcol,figlbl=figlbl,spread=spread,cunit=cunit,tunit=tunit)
          if (!is.null(ggdev)){
            suppressMessages(suppressWarnings(print(ggdev)))
            if (printOut=="TRUE") suppressMessages(suppressWarnings(ggsave(filename=paste(usrdir,"/Deviation_",figlbl,".",figFormat,sep=""),height=hth,width=wth,units="cm",dpi=200)))
            devplot[[length(devplot)+1]] <- ggdev
          }
          # NPDE plot
          npdeout <- nca.npde.plot(plotdata=plotdata,xvar="ID",npdecol=npdecol,figlbl=figlbl,cunit=cunit,tunit=tunit)
          if (is.null(npdeout$forestdata)) next
          forestdata <- npdeout$forestdata
          forestdata$str <- figlbl
          fpval <- rbind(fpval, forestdata)
          
          npdeplot[[length(npdeplot)+1]] <- npdeout$ggnpde
          suppressMessages(suppressWarnings(print(npdeout$ggnpde)))
          if (printOut=="TRUE") suppressMessages(suppressWarnings(ggsave(filename=paste(usrdir,"/NPDE_",figlbl,".",figFormat,sep=""),height=hth,width=wth,units="cm",dpi=200)))
          
          histnpdeplot[[length(histnpdeplot)+1]] <- npdeout$gghnpde
          suppressMessages(suppressWarnings(print(npdeout$gghnpde)))
          if (printOut=="TRUE") suppressMessages(suppressWarnings(ggsave(filename=paste(usrdir,"/histNPDE_",figlbl,".",figFormat,sep=""),height=hth,width=wth,units="cm",dpi=200)))
        }
      }
      
      if(noPlot=="FALSE"){
        # Forest plot for NPDE
        ggplt <- ggplot(fpval) + ggOpt_forest +
          xlab("\nNPDE") + ylab("") +
          labs(title = "Forest plot of NPDE\nErrorbar = 95% confidence interval\n\n") +
          geom_point(aes(mean,str,color="mean"), show_guide=T, size=2) +
          geom_errorbarh(aes(x=mean,y=str,xmin=mcil,xmax=mciu),size=0.6, color="red",height=0.3) +
          geom_point(aes(sdu,str,color="SD"), size=2) +
          geom_errorbarh(aes(x=sdu,y=str,xmin=sducil,xmax=sduciu), size=0.6, color="darkgreen", height=0.4)
        suppressMessages(suppressWarnings(print(ggplt)))
        forestplot[[length(forestplot)+1]] <- ggplt
        if (printOut=="TRUE") suppressMessages(suppressWarnings(ggsave(filename=paste(usrdir,"/forestNPDE.",figFormat,sep=""),height=hth,width=wth,units="cm",dpi=200)))
      }
    }else if (case==2){
      for (g in 1:ngrp){
        for (d in 1:ndose){
          if (nrow(dasdf[dasdf$GROUP==as.character(grp[g]) & dasdf$DoseNumber==dose[d],]) == 0) next
          tdasdf <- subset(dasdf, GROUP==as.character(grp[g]) & DoseNumber==dose[d])
          id     <- unique(tdasdf$ID)
          pde    <- data.frame()
          metric <- ""
          nout   <- 0
          figlbl <- paste(grNm,"-",as.character(grp[g]),"_",oidNm,"-",dose[d],sep="")
          for (i in 1:length(id)){
            obsdata <- subset(outData, ID==id[i] & GROUP==as.character(grp[g]) & DoseNumber==dose[d])
            simdata <- subset(tdasdf, ID==id[i])
            pdeout  <- nca.pde.deviation.outlier(obsdata=obsdata,simdata=simdata,idNm="ID",id=id[i],spread=spread,figlbl=figlbl,calcparam=alwprm,diagparam=param,cunit=cunit,tunit=tunit)
            pde     <- rbind(pde, cbind(data.frame(ID=id[i],GROUP=as.character(grp[g]),DoseNumber=dose[d]), pdeout$pde))
            outData[(outData$ID==id[i] & outData$GROUP==as.character(grp[g]) & outData$DoseNumber==dose[d]),] <- pdeout$obsdata
            if (pdeout$metric != ""){
              nout     <- nout + 1
              metric   <- paste(metric,pdeout$metric,sep=", ")
              if(noPlot=="FALSE"){
                gdr      <- pdeout$grob
                mylegend <- pdeout$legend
                lheight  <- pdeout$lheight
                if (printOut=="TRUE"){
                  fl <- paste(usrdir,"/Outlier_ID-",id[i],"_",figlbl,sep="")
                  eval(parse(text=paste(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",res=200)",sep="")))
                  suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
                  dev.off()
                }
                suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
                ggr <- grid.grab()
                outlierplot[[length(outlierplot)+1]] <- ggr
              }
            }
          }
          if (metric != "") metric <- gsub("^, ", "", metric)
          OTL  <- rbind(OTL, data.frame(No_of_outlier=nout,ID_metric=metric))
          npde <- rbind(npde,pde)
        }
      }
      npde   <- nca.npde(pdedata=npde,pdecol=alwprm)
      npdeNm <- paste("npde",alwprm,sep="")
      for (r in 1:nrow(outData)){
        if (nrow(npde[(npde$ID==outData$ID[r] & npde$GROUP==outData$GROUP[r] & npde$DoseNumber==outData$DoseNumber[r]),])!=1){
          outData[r,npdeNm] <- "NaN"
        }else{
          outData[r,paste("npde",alwprm,sep="")] <- npde[(npde$ID==outData$ID[r] & npde$GROUP==outData$GROUP[r] & npde$DoseNumber==outData$DoseNumber[r]),paste("npde",alwprm,sep="")]
        }
      }
      tmpdf <- outData
      names(outData)[c(1,3)] <- c(grNm,oidNm)
      outData <- as.data.frame(lapply(outData, FUN=function(x) round(as.numeric(x), digits=4)))
      if (printOut=="TRUE") write.table(outData, file=paste(usrdir,"/ncaOutput.tsv",sep=""), sep="\t", row.names=F, col.names=T, quote=F)
      outData <- tmpdf; rm(tmpdf)
      
      if(noPlot=="FALSE"){
        for (g in 1:ngrp){  
          for (d in 1:ndose){
            plotdata <- subset(outData, GROUP==as.character(grp[g]) & DoseNumber==dose[d])
            if (nrow(plotdata) == 0) next
            figlbl <- paste(grNm,"-",as.character(grp[g]),"_",oidNm,"-",dose[d],sep="")
            # Deviation plot
            ggdev <- nca.deviation.plot(plotdata=plotdata,xvar="ID",devcol=devcol,figlbl=figlbl,spread=spread,cunit=cunit,tunit=tunit)
            if (!is.null(ggdev)){
              suppressMessages(suppressWarnings(print(ggdev)))
              if (printOut=="TRUE") suppressMessages(suppressWarnings(ggsave(filename=paste(usrdir,"/Deviation_",figlbl,".",figFormat,sep=""),height=hth,width=wth,units="cm",dpi=200)))
              devplot[[length(devplot)+1]] <- ggdev
            }
            # NPDE plot
            npdeout <- nca.npde.plot(plotdata=plotdata,xvar="ID",npdecol=npdecol,figlbl=figlbl,cunit=cunit,tunit=tunit)
            if (is.null(npdeout$forestdata)) next
            forestdata <- npdeout$forestdata
            forestdata$str <- figlbl
            fpval <- rbind(fpval, forestdata)
            
            npdeplot[[length(npdeplot)+1]] <- npdeout$ggnpde
            suppressMessages(suppressWarnings(print(npdeout$ggnpde)))
            if (printOut=="TRUE") suppressMessages(suppressWarnings(ggsave(filename=paste(usrdir,"/NPDE_",figlbl,".",figFormat,sep=""),height=hth,width=wth,units="cm",dpi=200)))
            
            histnpdeplot[[length(histnpdeplot)+1]] <- npdeout$gghnpde
            suppressMessages(suppressWarnings(print(npdeout$gghnpde)))
            if (printOut=="TRUE") suppressMessages(suppressWarnings(ggsave(filename=paste(usrdir,"/histNPDE_",figlbl,".",figFormat,sep=""),height=hth,width=wth,units="cm",dpi=200)))
          }
        }
      }
      
      if(noPlot=="FALSE"){
        # Forest plot for NPDE
        ggplt <- ggplot(fpval) + ggOpt_forest +
          xlab("\nNPDE") + ylab("") +
          labs(title = "Forest plot of NPDE\nErrorbar = 95% confidence interval\n\n") +
          geom_point(aes(mean,str,color="mean"), show_guide=T, size=2) +
          geom_errorbarh(aes(x=mean,y=str,xmin=mcil,xmax=mciu),size=0.6, color="red",height=0.3) +
          geom_point(aes(sdu,str,color="SD"), size=2) +
          geom_errorbarh(aes(x=sdu,y=str,xmin=sducil,xmax=sduciu), size=0.6, color="darkgreen", height=0.4)
        suppressMessages(suppressWarnings(print(ggplt)))
        forestplot[[length(forestplot)+1]] <- ggplt
        if (printOut=="TRUE") suppressMessages(suppressWarnings(ggsave(filename=paste(usrdir,"/forestNPDE.",figFormat,sep=""),height=hth,width=wth,units="cm",dpi=200)))
      }
    }else if (case == 3){
      for (f in 1:nflag){
        for (d in 1:ndose){
          if (nrow(dasdf[dasdf$FLAG==as.character(flag[f]) & dasdf$DoseNumber==dose[d],]) == 0) next
          tdasdf <- subset(dasdf, FLAG==as.character(flag[f]) & DoseNumber==dose[d])
          id     <- unique(tdasdf$ID)
          pde    <- data.frame()
          metric <- ""
          nout   <- 0
          for (i in 1:length(id)){
            obsdata <- subset(outData, ID==id[i] & FLAG==as.character(flag[f]) & DoseNumber==dose[d])
            simdata <- subset(tdasdf, ID==id[i])
            figlbl  <- paste(flNm,"-",as.character(flag[f]),"_",oidNm,"-",dose[d],sep="")
            pdeout  <- nca.pde.deviation.outlier(obsdata=obsdata,simdata=simdata,idNm="ID",id=id[i],spread=spread,figlbl=figlbl,calcparam=alwprm,diagparam=param,cunit=cunit,tunit=tunit)
            pde     <- rbind(pde, cbind(data.frame(ID=id[i],FLAG=as.character(flag[f]),DoseNumber=dose[d]), pdeout$pde))
            outData[(outData$ID==id[i] & outData$FLAG==as.character(flag[f]) & outData$DoseNumber==dose[d]),] <- pdeout$obsdata
            if (pdeout$metric != ""){
              nout     <- nout + 1
              metric   <- paste(metric,pdeout$metric,sep=", ")
              if(noPlot=="FALSE"){
                gdr      <- pdeout$grob
                mylegend <- pdeout$legend
                lheight  <- pdeout$lheight
                if (printOut=="TRUE"){
                  fl <- paste(usrdir,"/Outlier_ID-",id[i],"_",figlbl,sep="")
                  eval(parse(text=paste(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",res=200)",sep="")))
                  suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
                  dev.off()
                }
                suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
                ggr <- grid.grab()
                outlierplot[[length(outlierplot)+1]]
              }
            }
          }
          if (metric != "") metric <- gsub("^, ", "", metric)
          OTL  <- rbind(OTL, data.frame(No_of_outlier=nout,ID_metric=metric))
          npde <- rbind(npde,pde)
        }
      }
      npde   <- nca.npde(pdedata=npde,pdecol=alwprm)
      npdeNm <- paste("npde",alwprm,sep="")
      for (r in 1:nrow(outData)){
        if (nrow(npde[(npde$ID==outData$ID[r] & npde$FLAG==outData$FLAG[r] & npde$DoseNumber==outData$DoseNumber[r]),])!=1){
          outData[r,npdeNm] <- "NaN"
        }else{
          outData[r,paste("npde",alwprm,sep="")] <- npde[(npde$ID==outData$ID[r] & npde$FLAG==outData$FLAG[r] & npde$DoseNumber==outData$DoseNumber[r]),paste("npde",alwprm,sep="")]
        }
      }
      tmpdf <- outData
      names(outData)[c(1,3)] <- c(flNm,oidNm)
      outData <- as.data.frame(lapply(outData, FUN=function(x) round(as.numeric(x), digits=4)))
      if (printOut=="TRUE") write.table(outData, file=paste(usrdir,"/ncaOutput.tsv",sep=""), sep="\t", row.names=F, col.names=T, quote=F)
      outData <- tmpdf; rm(tmpdf)
      
      if(noPlot=="FALSE"){
        for (f in 1:nflag){  
          for (d in 1:ndose){
            plotdata <- subset(outData, FLAG==as.character(flag[f]) & DoseNumber==dose[d])
            if (nrow(plotdata) == 0) next
            figlbl <- paste(flNm,"-",as.character(flag[f]),"_",oidNm,"-",dose[d],sep="")
            # Deviation plot
            ggdev <- nca.deviation.plot(plotdata=plotdata,xvar="ID",devcol=devcol,figlbl=figlbl,spread=spread,cunit=cunit,tunit=tunit)
            if (!is.null(ggdev)){
              suppressMessages(suppressWarnings(print(ggdev)))
              if (printOut=="TRUE") suppressMessages(suppressWarnings(ggsave(filename=paste(usrdir,"/Deviation_",figlbl,".",figFormat,sep=""),height=hth,width=wth,units="cm",dpi=200)))
              devplot[[length(devplot)+1]] <- ggdev
            }
            # NPDE plot
            npdeout <- nca.npde.plot(plotdata=plotdata,xvar="ID",npdecol=npdecol,figlbl=figlbl,cunit=cunit,tunit=tunit)
            if (is.null(npdeout$forestdata)) next
            forestdata <- npdeout$forestdata
            forestdata$str <- figlbl
            fpval <- rbind(fpval, forestdata)
            
            npdeplot[[length(npdeplot)+1]] <- npdeout$ggnpde
            suppressMessages(suppressWarnings(print(npdeout$ggnpde)))
            if (printOut=="TRUE") suppressMessages(suppressWarnings(ggsave(filename=paste(usrdir,"/NPDE_",figlbl,".",figFormat,sep=""),height=hth,width=wth,units="cm",dpi=200)))
            
            histnpdeplot[[length(histnpdeplot)+1]] <- npdeout$gghnpde
            suppressMessages(suppressWarnings(print(npdeout$gghnpde)))
            if (printOut=="TRUE") suppressMessages(suppressWarnings(ggsave(filename=paste(usrdir,"/histNPDE_",figlbl,".",figFormat,sep=""),height=hth,width=wth,units="cm",dpi=200)))
          }
        }
      }
      
      if(noPlot=="FALSE"){
        # Forest plot for NPDE
        ggplt <- ggplot(fpval) + ggOpt_forest +
          xlab("\nNPDE") + ylab("") +
          labs(title = "Forest plot of NPDE\nErrorbar = 95% confidence interval\n\n") +
          geom_point(aes(mean,str,color="mean"), show_guide=T, size=2) +
          geom_errorbarh(aes(x=mean,y=str,xmin=mcil,xmax=mciu),size=0.6, color="red",height=0.3) +
          geom_point(aes(sdu,str,color="SD"), size=2) +
          geom_errorbarh(aes(x=sdu,y=str,xmin=sducil,xmax=sduciu), size=0.6, color="darkgreen", height=0.4)
        suppressMessages(suppressWarnings(print(ggplt)))
        forestplot[[length(forestplot)+1]] <- ggplt
        if (printOut=="TRUE") suppressMessages(suppressWarnings(ggsave(filename=paste(usrdir,"/forestNPDE.",figFormat,sep=""),height=hth,width=wth,units="cm",dpi=200)))
      }
    }else if (case == 4){
      for (g in 1:ngrp){
        for (f in 1:nflag){
          for (d in 1:ndose){
            if (nrow(dasdf[dasdf$GROUP==as.character(grp[g]) & dasdf$FLAG==as.character(flag[f]) & dasdf$DoseNumber==dose[d],]) == 0) next
            tdasdf <- subset(dasdf, GROUP==as.character(grp[g]) & FLAG==as.character(flag[f]) & DoseNumber==dose[d])
            id     <- unique(tdasdf$ID)
            pde    <- data.frame()
            metric <- ""
            nout   <- 0
            for (i in 1:length(id)){
              obsdata <- subset(outData, ID==id[i] & GROUP==as.character(grp[g]) & FLAG==as.character(flag[f]) & DoseNumber==dose[d])
              simdata <- subset(tdasdf, ID==id[i])
              if(nrow(obsdata)==0 | nrow(simdata)==0) next
              figlbl  <- paste(grNm,"-",as.character(grp[g]),"_",flNm,"-",as.character(flag[f]),"_",oidNm,"-",dose[d],sep="")
              pdeout  <- nca.pde.deviation.outlier(obsdata=obsdata,simdata=simdata,idNm="ID",id=id[i],spread=spread,figlbl=figlbl,calcparam=alwprm,diagparam=param,cunit=cunit,tunit=tunit)
              pde     <- rbind(pde, cbind(data.frame(ID=id[i],GROUP=as.character(grp[g]),FLAG=as.character(flag[f]),DoseNumber=dose[d]), pdeout$pde))
              outData[(outData$ID==id[i] & outData$GROUP==as.character(grp[g]) & outData$FLAG==as.character(flag[f]) & outData$DoseNumber==dose[d]),] <- pdeout$obsdata
              if (pdeout$metric != ""){
                nout     <- nout + 1
                metric   <- paste(metric,pdeout$metric,sep=", ")
                if(noPlot=="FALSE"){
                  gdr      <- pdeout$grob
                  mylegend <- pdeout$legend
                  lheight  <- pdeout$lheight
                  if (printOut=="TRUE"){
                    fl <- paste(usrdir,"/Outlier_ID-",id[i],"_",figlbl,sep="")
                    eval(parse(text=paste(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",res=200)",sep="")))
                    suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
                    dev.off()
                  }
                  suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
                  ggr <- grid.grab()
                  outlierplot[[length(outlierplot)+1]] <- ggr
                }
              }
            }
            if (metric != "") metric <- gsub("^, ", "", metric)
            OTL  <- rbind(OTL, data.frame(No_of_outlier=nout,ID_metric=metric))
            npde <- rbind(npde,pde)
          }
        }
      }
      npde   <- nca.npde(pdedata=npde,pdecol=alwprm)
      npdeNm <- paste("npde",alwprm,sep="")
      for (r in 1:nrow(outData)){
        if (nrow(npde[(npde$ID==outData$ID[r] & npde$GROUP==outData$GROUP[r] & npde$FLAG==outData$FLAG[r] & npde$DoseNumber==outData$DoseNumber[r]),])!=1){
          outData[r,npdeNm] <- "NaN"
        }else{
          outData[r,paste("npde",alwprm,sep="")] <- npde[(npde$ID==outData$ID[r] & npde$GROUP==outData$GROUP[r] & npde$FLAG==outData$FLAG[r] & npde$DoseNumber==outData$DoseNumber[r]),paste("npde",alwprm,sep="")]
        }
      }
      tmpdf <- outData
      names(outData)[c(1,3,4)] <- c(grNm,flNm,oidNm)
      outData <- as.data.frame(lapply(outData, FUN=function(x) round(as.numeric(x), digits=4)))
      if (printOut=="TRUE") write.table(outData, file=paste(usrdir,"/ncaOutput.tsv",sep=""), sep="\t", row.names=F, col.names=T, quote=F)
      outData <- tmpdf; rm(tmpdf)
      
      if(noPlot=="FALSE"){
        for (g in 1:ngrp){
          for (f in 1:nflag){  
            for (d in 1:ndose){
              plotdata <- subset(outData, GROUP==as.character(grp[g]) & FLAG==as.character(flag[f]) & DoseNumber==dose[d])
              if (nrow(plotdata) == 0) next
              figlbl <- paste(grNm,"-",as.character(grp[g]),"_",flNm,"-",as.character(flag[f]),"_",oidNm,"-",dose[d],sep="")
              # Deviation plot
              ggdev <- nca.deviation.plot(plotdata=plotdata,xvar="ID",devcol=devcol,figlbl=figlbl,spread=spread,cunit=cunit,tunit=tunit)
              if (!is.null(ggdev)){
                suppressMessages(suppressWarnings(print(ggdev)))
                if (printOut=="TRUE") suppressMessages(suppressWarnings(ggsave(filename=paste(usrdir,"/Deviation_",figlbl,".",figFormat,sep=""),height=hth,width=wth,units="cm",dpi=200)))
                devplot[[length(devplot)+1]] <- ggdev
              }
              # NPDE plot
              npdeout <- nca.npde.plot(plotdata=plotdata,xvar="ID",npdecol=npdecol,figlbl=figlbl,cunit=cunit,tunit=tunit)
              if (is.null(npdeout$forestdata)) next
              forestdata <- npdeout$forestdata
              forestdata$str <- figlbl
              fpval <- rbind(fpval, forestdata)
              
              npdeplot[[length(npdeplot)+1]] <- npdeout$ggnpde
              suppressMessages(suppressWarnings(print(npdeout$ggnpde)))
              if (printOut=="TRUE") suppressMessages(suppressWarnings(ggsave(filename=paste(usrdir,"/NPDE_",figlbl,".",figFormat,sep=""),height=hth,width=wth,units="cm",dpi=200)))
              
              histnpdeplot[[length(histnpdeplot)+1]] <- npdeout$gghnpde
              suppressMessages(suppressWarnings(print(npdeout$gghnpde)))
              if (printOut=="TRUE") suppressMessages(suppressWarnings(ggsave(filename=paste(usrdir,"/histNPDE_",figlbl,".",figFormat,sep=""),height=hth,width=wth,units="cm",dpi=200)))
            }
          }
        }
      }
      
      if(noPlot=="FALSE"){
        # Forest plot for NPDE
        ggplt <- ggplot(fpval) + ggOpt_forest +
          xlab("\nNPDE") + ylab("") +
          labs(title = "Forest plot of NPDE\nErrorbar = 95% confidence interval\n\n") +
          geom_point(aes(mean,str,color="mean"), show_guide=T, size=2) +
          geom_errorbarh(aes(x=mean,y=str,xmin=mcil,xmax=mciu),size=0.6, color="red",height=0.3) +
          geom_point(aes(sdu,str,color="SD"), size=2) +
          geom_errorbarh(aes(x=sdu,y=str,xmin=sducil,xmax=sduciu), size=0.6, color="darkgreen", height=0.4)
        suppressMessages(suppressWarnings(print(ggplt)))
        forestplot[[length(forestplot)+1]] <- ggplt
        if (printOut=="TRUE") suppressMessages(suppressWarnings(ggsave(filename=paste(usrdir,"/forestNPDE.",figFormat,sep=""),height=hth,width=wth,units="cm",dpi=200)))
      }
    }
    
    # Statistical analysis for each patient group
    stNm <- c("Tmax","Cmax","AUClast","AUClower_upper","AUCINF_obs","AUC_pExtrap_obs","AUCINF_pred","AUC_pExtrap_pred","AUMClast","AUMCINF_obs","AUMC_pExtrap_obs","AUMCINF_pred","AUMC_pExtrap_pred","HL_Lambda_z","Rsq","Rsq_adjusted","No_points_Lambda_z")
    simGrStat <- data.frame()
    if (case == 1){
      for (d in 1:ndose){
        pm <- data.frame(Ntot=numeric(0),Nunique=numeric(0),Min=numeric(0),Max=numeric(0),Mean=numeric(0),Median=numeric(0),SD=numeric(0),SE=numeric(0),CVp=numeric(0),CI95l=numeric(0),CI95u=numeric(0),gMean=numeric(0),gCVp=numeric(0))
        nm <- data.frame(character(0))
        counter <- 1
        for (i in 1:length(stNm)){
          Nm <- stNm[i]
          tdf <- as.numeric(as.character(dasdf[(dasdf$DoseNumber==dose[d] & dasdf[,Nm]!="NaN"),Nm]))
          if (length(tdf) < 2){
            nm <- rbind(nm, data.frame(Nm))
            pm[counter,] <- rep(NA,13)
          }else{
            nm <- rbind(nm, data.frame(Nm))
            stPrm <- calc.stat(x=tdf)        # calls calc.stat function
            stPrm <- unname(stPrm)
            pm[counter,] <- stPrm
          }
          counter <- counter + 1
        }
        if (nrow(pm) == 0) next
        pm <- data.frame(lapply(pm, function(x){if(is.numeric(x)){round(x,digits=4)}else{x}}))
        tmpStat <- t(cbind(nm, pm))
        rownames(tmpStat)[1] <- "Name"
        tmpStat <- cbind(DOSE=dose[d],Stat=rownames(tmpStat),tmpStat)
        for (cnum in 3:ncol(tmpStat)){colnames(tmpStat)[cnum] <- as.character(tmpStat[1,cnum])}
        tmpStat <- tmpStat[-1,]
        simGrStat <- rbind(simGrStat,tmpStat)
      }
      names(simGrStat)[1] <- oidNm
    }
    if (case == 2){
      simGrStat <- data.frame()
      for (g in 1:ngrp){
        for (d in 1:ndose){
          pm <- data.frame(Ntot=numeric(0),Nunique=numeric(0),Min=numeric(0),Max=numeric(0),Mean=numeric(0),Median=numeric(0),SD=numeric(0),SE=numeric(0),CVp=numeric(0),CI95l=numeric(0),CI95u=numeric(0),gMean=numeric(0),gCVp=numeric(0))
          nm <- data.frame(character(0))
          counter <- 1
          for (i in 1:length(stNm)){
            Nm <- stNm[i]
            tdf <- as.numeric(as.character(dasdf[(dasdf$DoseNumber==dose[d] & dasdf$GROUP==as.character(grp[g]) & dasdf[,Nm]!="NaN"),Nm]))
            if (length(tdf) < 2){
              nm <- rbind(nm, data.frame(Nm))
              pm[counter,] <- rep(NA,13)
            }else{
              nm <- rbind(nm, data.frame(Nm))
              stPrm <- calc.stat(x=tdf)        # calls calc.stat function
              stPrm <- unname(stPrm)
              pm[counter,] <- stPrm
            }
            counter <- counter + 1
          }
          if (nrow(pm) == 0) next
          pm <- data.frame(lapply(pm, function(x){if(is.numeric(x)){round(x,digits=4)}else{x}}))
          tmpStat <- t(cbind(nm, pm))
          rownames(tmpStat)[1] <- "Name"
          tmpStat <- cbind(GRP=as.character(grp[g]),DOSE=dose[d],Stat=rownames(tmpStat),tmpStat)
          for (cnum in 4:ncol(tmpStat)){colnames(tmpStat)[cnum] <- as.character(tmpStat[1,cnum])}
          tmpStat <- tmpStat[-1,]
          simGrStat <- rbind(simGrStat,tmpStat)
        }
      }
      names(simGrStat)[c(1,2)] <- c(grNm,oidNm)
    }
    if (case == 3){
      simGrStat <- data.frame()
      for (f in 1:nflag){
        for (d in 1:ndose){
          pm <- data.frame(Ntot=numeric(0),Nunique=numeric(0),Min=numeric(0),Max=numeric(0),Mean=numeric(0),Median=numeric(0),SD=numeric(0),SE=numeric(0),CVp=numeric(0),CI95l=numeric(0),CI95u=numeric(0),gMean=numeric(0),gCVp=numeric(0))
          nm <- data.frame(character(0))
          counter <- 1
          for (i in 1:length(stNm)){
            Nm <- stNm[i]
            tdf <- as.numeric(as.character(dasdf[(dasdf$DoseNumber==dose[d] & dasdf$FLAG==as.character(flag[f]) & dasdf[,Nm]!="NaN"),Nm]))
            if (length(tdf) < 2){
              nm <- rbind(nm, data.frame(Nm))
              pm[counter,] <- rep(NA,13)
            }else{
              nm <- rbind(nm, data.frame(Nm))
              stPrm <- calc.stat(x=tdf)        # calls calc.stat function
              stPrm <- unname(stPrm)
              pm[counter,] <- stPrm
            }
            counter <- counter + 1
          }
          if (nrow(pm) == 0) next
          pm <- data.frame(lapply(pm, function(x){if(is.numeric(x)){round(x,digits=4)}else{x}}))
          tmpStat <- t(cbind(nm, pm))
          rownames(tmpStat)[1] <- "Name"
          tmpStat <- cbind(GRP=as.character(flag[f]),DOSE=dose[d],Stat=rownames(tmpStat),tmpStat)
          for (cnum in 4:ncol(tmpStat)){colnames(tmpStat)[cnum] <- as.character(tmpStat[1,cnum])}
          tmpStat <- tmpStat[-1,]
          simGrStat <- rbind(simGrStat,tmpStat)
        }
      }
      names(simGrStat)[c(1,2)] <- c(flNm,oidNm)
    }
    if (case == 4){
      simGrStat <- data.frame()
      for (g in 1:ngrp){
        for (f in 1:nflag){
          for (d in 1:ndose){
            pm <- data.frame(Ntot=numeric(0),Nunique=numeric(0),Min=numeric(0),Max=numeric(0),Mean=numeric(0),Median=numeric(0),SD=numeric(0),SE=numeric(0),CVp=numeric(0),CI95l=numeric(0),CI95u=numeric(0),gMean=numeric(0),gCVp=numeric(0))
            nm <- data.frame(character(0))
            counter <- 1
            for (i in 1:length(stNm)){
              Nm <- stNm[i]
              tdf <- as.numeric(as.character(dasdf[(dasdf$DoseNumber==dose[d] & dasdf$GROUP==as.character(grp[g]) & dasdf$FLAG==as.character(flag[f]) & dasdf[,Nm]!="NaN"),Nm]))
              if (length(tdf) < 2){
                nm <- rbind(nm, data.frame(Nm))
                pm[counter,] <- rep(NA,13)
              }else{
                nm <- rbind(nm, data.frame(Nm))
                stPrm <- calc.stat(x=tdf)        # calls calc.stat function
                stPrm <- unname(stPrm)
                pm[counter,] <- stPrm
              }
              counter <- counter + 1
            }
            if (nrow(pm) == 0) next
            pm <- data.frame(lapply(pm, function(x){if(is.numeric(x)){round(x,digits=4)}else{x}}))
            tmpStat <- t(cbind(nm, pm))
            rownames(tmpStat)[1] <- "Name"
            tmpStat <- cbind(GRP=as.character(grp[g]),FLAG=as.character(flag[f]),DOSE=dose[d],Stat=rownames(tmpStat),tmpStat)
            for (cnum in 5:ncol(tmpStat)){colnames(tmpStat)[cnum] <- as.character(tmpStat[1,cnum])}
            tmpStat <- tmpStat[-1,]
            simGrStat <- rbind(simGrStat,tmpStat)
          }
        }
      }
      names(simGrStat)[c(1:3)] <- c(grNm,flNm,oidNm)
    }
    
    if (printOut=="TRUE") write.table(simGrStat, file=paste(usrdir,"/SimStat.tsv",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
    
    streamsEnv <- parent.frame()
    if(exists("outData"))   assign("ncaOutput",  outData,   envir=streamsEnv)
    if(exists("grStat"))    assign("ObsStat",    grStat,    envir=streamsEnv)
    if(exists("simGrStat")) assign("SimStat",    simGrStat, envir=streamsEnv)
    if(exists("nmdf"))      assign("ncaSimData", nmdf,      envir=streamsEnv)
    if(exists("dasdf"))     assign("ncaSimEst",  dasdf,     envir=streamsEnv)
    
    if (printOut=="TRUE"){
      # Create HTML output
      simFileNm <- ifelse(is.data.frame(simFile), deparse(substitute(simFile)), simFile)
      txt <- paste(txt,paste0("Name of the NONMEM simulation output file: \"",simFileNm,"\""),sep="\n")
      txt <- paste(txt,paste0("Number of simulations performed: ",nsim),sep="\n")
      pddf <- cbind(pddf,OTL); names(pddf)[c((ncol(pddf)-1),ncol(pddf))] <- c("No. of outlier","Selected outliers ID and NCA metrics")
      #for (i in 1:length(concplot)){print(concplot[i])}
      #for (i in 1:length(histobsplot)){print(histobsplot[i])}
      #for (i in 1:length(popplot)){print(popplot[i])}
      #for (i in 1:length(devplot)){print(devplot[i])}
      #for (i in 1:length(outlierplot)){print(outlierplot[i])}
      #for (i in 1:length(forestplot)){print(forestplot[i])}
      #for (i in 1:length(npdeplot)){print(npdeplot[i])}
      #for (i in 1:length(histnpdeplot)){print(histnpdeplot[i])}
      
      if (case == 1){
        prnTab <- head(cbind(outData[,1:3], subset(outData, select = tabCol)),100)
        prnTab[,4:ncol(prnTab)] <- data.frame(lapply(prnTab[,4:ncol(prnTab)], function(x) round(as.numeric(x),digits=4)))
        names(prnTab)[2] <- oidNm
      }else if (case == 2){
        prnTab <- head(cbind(outData[,1:4], subset(outData, select = tabCol)),100)
        prnTab[,5:ncol(prnTab)] <- data.frame(lapply(prnTab[,5:ncol(prnTab)], function(x) round(as.numeric(x),digits=4)))
        names(prnTab)[c(1,3)] <- c(grNm,oidNm)
      }else if (case == 3){
        prnTab <- head(cbind(outData[,1:4], subset(outData, select = tabCol)),100)
        prnTab[,5:ncol(prnTab)] <- data.frame(lapply(prnTab[,5:ncol(prnTab)], function(x) round(as.numeric(x),digits=4)))
        names(prnTab)[c(1,3)] <- c(flNm,oidNm)
      }else if (case == 4){
        prnTab <- head(cbind(outData[,1:5], subset(outData, select = tabCol)),100)
        prnTab[,6:ncol(prnTab)] <- data.frame(lapply(prnTab[,6:ncol(prnTab)], function(x) round(as.numeric(x),digits=4)))
        names(prnTab)[c(1,2,4)] <- c(grNm,flNm,oidNm)
      }
      prnTab <- data.frame(lapply(prnTab, function(x){if(is.numeric(x)){round(x,digits=4)}else{x}}))
      fnOut <- list(arglist=match.call(),TXT=txt, pddf=pddf, prnTab=prnTab, NSIM=nsim, spread=spread, conc=concplot, histobs=histobsplot, pop=popplot, dev=devplot, outlier=outlierplot, forest=forestplot, npde=npdeplot, histnpde=histnpdeplot, phth=phth, pwth=pwth)
    }
  }
  setwd(usrdir)
  
  if (printOut=="TRUE"){
    misc <- system.file("misc", package = "ncappc")
    if (is.null(simFile)){
      mdFile <- paste(misc,"ncappcReport-NCA.Rmd",sep="/")
      nwFile <- paste(misc,"ncappcReport-NCA.Rnw",sep="/")
    }else{
      mdFile <- paste(misc,"ncappcReport-NCA-PPC.Rmd",sep="/")
      nwFile <- paste(misc,"ncappcReport-NCA-PPC.Rnw",sep="/")
    }
    
    knit2html(mdFile, style=paste(misc,"custom.css",sep="/"))
    knit(nwFile)
    if (.Platform$OS.type == "unix"){
      texcomp <- system('which texi2pdf')
      if (texcomp == 0){
        knit2pdf(nwFile)
      }else{
        print("Please install \"texi2pdf\" to compile the produced tex file into a PDF report")
      }
    }else if (.Platform$OS.type == "windows"){
      texcomp <- system('kpsewhich pdftex --version')
      if (texcomp == 0){
        knit2pdf(nwFile)
      }else{
        print("Please install \"pdftex\" to compile the produced tex file into a PDF report")
      }
    }
  }
  #unlink(list.files(pattern = "ncappcReport.[a,t,l,m,o]"))
  #unlink(list.files(pattern = "sum.tex"))
  #unlink(list.files(pattern = "tab.tex"))
}
