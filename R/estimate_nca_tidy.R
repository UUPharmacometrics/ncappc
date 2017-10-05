if(getRversion() >= '2.15.1')
  utils::globalVariables(c("."))

estimate_nca_tidy <- function(pkData, 
                              all_data,
                              doseAmtNm,
                              dvLog, dataType,
                              idNm, timeNm, concNm,
                              adminType, TI,
                              dateColNm, dateFormat, timeFormat,
                              backExtrp,negConcExcl,doseType,
                              method,AUCTimeRange,LambdaTimeRange,
                              LambdaExclude,doseTime,Tau,simFile,onlyNCA,
                              strat_vars = NULL,
                              #popStrNm1=NULL,popStrNm2=NULL,popStrNm3=NULL,
                              dunit=NULL,
                              extrapolate=FALSE,
                              ...) 
{
  
  # for CRAN checks
  dose <- time <- conc <- NULL
  
  pddf <- data.frame()   # Summary table
  outData   <- data.frame()  # Create empty data frame for output
  cdata <- data.frame()
  
  if(is.null(dunit)){Dcol   <- "Dose"}else{Dcol <- paste0("Dose (",dunit,")")} # Dose
  
  obsData <- tibble::as_tibble(pkData)
  fullData <- tibble::as_tibble(all_data)
  
  #if (nrow(obsData) == 0){next}
  
  #id_name <- rlang::quo(!!rlang::sym(idNm))
  id_name <- rlang::sym(idNm)
  n_label <- "No. of individuals"
  
  strats <- NULL
  if(!is.null(strat_vars)) strats <- rlang::syms(strat_vars)
  
  dose_amt_name <- doseAmtNm
  #if(!is.null(doseAmtNm)) dose_amt_name <- rlang::quo(!!rlang::sym(doseAmtNm))
  if(!is.null(doseAmtNm)) dose_amt_name <- rlang::sym(doseAmtNm)
  
  tmp_data <- fullData %>%  
    dplyr::select(c(!!id_name,!!dose_amt_name,!!!strats))
  if(!is.null(dose_amt_name)) tmp_data <- tmp_data %>% 
    dplyr::filter((!!dose_amt_name) > 0) 
  
  ## create DF of a summary of groups
  pddf <- tmp_data
  if(!is.null(strats)) pddf <- pddf %>% dplyr::group_by(!!!strats,add=TRUE)
  pddf <- pddf %>% dplyr::distinct() 
  if(!is.null(dose_amt_name)){
    tmp_fcn <- function(x){paste(sort(unique(x)),collapse = ", ")}
    pddf <- pddf %>% dplyr::summarize(!!Dcol:=tmp_fcn(!!dose_amt_name),!!n_label:=n())
  }  else {
    pddf <- pddf %>% dplyr::summarize(!!Dcol:=NA,!!n_label:=n())
  }
  pddf <- pddf %>% data.frame(check.names=FALSE)
  
  if(!is.null(dose_amt_name)){
    ind_amt_data <- tmp_data %>% dplyr::group_by(!!id_name) 
    if(!is.null(strats)) ind_amt_data <- ind_amt_data %>% dplyr::group_by(!!!strats,add=TRUE)
    
    ind_amt_data <- ind_amt_data %>% 
      dplyr::distinct() %>% 
      dplyr::summarise("n"=n(),"dose"=first(!!dose_amt_name)) %>% 
      dplyr::mutate(ind_amt=as.double(dplyr::if_else(n==1,as.character(dose),"NA")))
    
    join_vars <- c(id_name,strats) %>% 
      purrr::map(rlang::UQE) %>% paste()

    obsData <- dplyr::left_join(obsData,ind_amt_data,by=join_vars)
    
  } else {
    obsData$ind_amt <- NA
  }
  
  # }
  
  new_data <- obsData %>% dplyr::group_by(!!id_name) 
  if(!is.null(strats)) new_data <- new_data %>% dplyr::group_by(!!!strats,add=TRUE)
  new_data <- new_data %>% 
    dplyr::do(data.frame(nca_ind_data(., dvLog = dvLog, dataType=dataType,
                                      idNm=idNm, timeNm=timeNm, concNm=concNm,
                                      adminType=adminType, TI=TI,
                                      dateColNm=dateColNm, dateFormat=dateFormat, timeFormat=timeFormat),ind_amt=.$ind_amt[1]))
  
  outData <- new_data %>% dplyr::group_by(!!id_name)
  if(!is.null(strats)) outData <- outData %>% dplyr::group_by(!!!strats,add=TRUE)
  outData <- outData %>%
    dplyr::do(data.frame(Dose=.$ind_amt[1],t(est.nca(time=.$time,conc=.$conc,backExtrp=backExtrp,negConcExcl=negConcExcl,doseType=doseType,adminType=adminType,
                                                     doseAmt=.$ind_amt[1],method=method,AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,
                                                     LambdaExclude=LambdaExclude,doseTime=doseTime,Tau=Tau,TI=.$iTI[1],simFile=simFile,dset=dataType,onlyNCA=onlyNCA,extrapolate=extrapolate)))) 
  if(!is.null(strats)) for (i in 1:length(strats)) outData <- outData %>% dplyr::rename(!!paste0("STRAT",i):=!!strats[[i]])

  outData <- outData %>% data.frame(check.names=FALSE)
  
  
  cdata   <- new_data %>% 
    dplyr::select(c(!!id_name,time,conc,!!!strats)) %>% 
    dplyr::rename("ID"=!!id_name,"Time"=time,"Conc"=conc) 
  if(!is.null(strats)) for (i in 1:length(strats)) cdata <- cdata %>% dplyr::rename(!!paste0("STRAT",i):=!!strats[[i]])
  cdata <- cdata %>% data.frame(check.names=FALSE)
  
  return(list(outData=outData,pddf=pddf,cdata=cdata))
}
