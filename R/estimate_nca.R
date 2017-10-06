estimate_nca <- function(case,
                         pkData, 
                         all_data,
                         doseAmtNm,
                         dvLog, dataType,
                         idNm, timeNm, concNm,
                         adminType, TI,
                         dateColNm, dateFormat, timeFormat,
                         backExtrp,negConcExcl,doseType,
                         method,AUCTimeRange,LambdaTimeRange,
                         LambdaExclude,doseTime,Tau,simFile,onlyNCA,
                         npopStr1,npopStr2,npopStr3,
                         popStrNm1,popStrNm2,popStrNm3,
                         popStr1,popStr2,popStr3,
                         dunit=NULL,
                         nca_method=1,
                         tidy=T,
                         extrapolate=FALSE,
                         ...
                         ) {
  
  pddf <- data.frame()   # Summary table
  outData   <- data.frame()  # Create empty data frame for output
  cdata <- data.frame()
  
  if(is.null(dunit)){Dcol   <- "Dose"}else{Dcol <- paste0("Dose (",dunit,")")}                       # Dose

  if(tidy){
    result <- estimate_nca_tidy(pkData, 
                                all_data,
                                doseAmtNm,
                                dvLog, dataType,
                                idNm, timeNm, concNm,
                                adminType, TI,
                                dateColNm, dateFormat, timeFormat,
                                backExtrp,negConcExcl,doseType,
                                method,AUCTimeRange,LambdaTimeRange,
                                LambdaExclude,doseTime,Tau,simFile,onlyNCA,
                                strat_vars = c(popStrNm1,popStrNm2,popStrNm3),
                                dunit,extrapolate=extrapolate,
                                ...) 
    return(result)
  }
  
  
  # Estimate NCA metrics for case = 1
  if (case == 1){
    if(nca_method==1){
      ifdf <- pkData
      if (nrow(ifdf) == 0){next}
      
      idd  <- unique(as.character(ifdf[,idNm]))
      if (is.null(doseAmtNm)){
        doseAmount <- NA
      }else{
        doseData   <- as.numeric(as.character(all_data[,doseAmtNm]))
        doseData   <- doseData[complete.cases(doseData) & doseData>0]
        doseAmount <- ifelse(length(doseData)==0, NA, paste(unique(doseData), collapse=", "))
      }
      
      # Description
      pddf  <- rbind(pddf, data.frame(a=doseAmount, b=length(idd)))
      for (i in 1:length(idd)){
        if (is.null(doseAmtNm)){
          idzAmt <- NA
        }else{
          doseData <- as.numeric(as.character(all_data[all_data[,idNm]==idd[i], doseAmtNm]))
          doseData <- doseData[complete.cases(doseData) & doseData>0]
          idzAmt   <- ifelse(length(doseData)==0, NA, doseData[1])
        }
        tcTI <- nca.ind.data(pkData=ifdf, ID=idd[i], dvLog = dvLog, dataType=dataType,
                             idNm=idNm, timeNm=timeNm, concNm=concNm,
                             adminType=adminType, TI=TI,
                             dateColNm=dateColNm, dateFormat=dateFormat, timeFormat=timeFormat)
        tc      <- tcTI$tc
        iTI     <- tcTI$iTI
        if (nrow(tc)==0) next
        time    <- as.numeric(tc$time)
        conc    <- as.numeric(tc$conc)
        cdata   <- rbind(cdata,cbind(Time=time,Conc=conc,ID=idd[i]))
        NCAprm  <- est.nca(time=time,conc=conc,backExtrp=backExtrp,negConcExcl=negConcExcl,
                           doseType=doseType,adminType=adminType,
                           doseAmt=idzAmt,method=method,AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,
                           LambdaExclude=LambdaExclude,doseTime=doseTime,Tau=Tau,TI=iTI,simFile=simFile,
                           dset=dataType,onlyNCA=onlyNCA,extrapolate=extrapolate)
        outData <- rbind(outData, data.frame(ID=idd[i],Dose=idzAmt,t(NCAprm)))
      }
    }
    if(nca_method==3){
      
      ifdf <- pkData
      if (nrow(ifdf) == 0){next}
      
      idd  <- unique(as.character(ifdf[,idNm]))
      if (is.null(doseAmtNm)){
        doseAmount <- NA
      }else{
        doseData   <- as.numeric(as.character(all_data[,doseAmtNm]))
        doseData   <- doseData[complete.cases(doseData) & doseData>0]
        doseAmount <- ifelse(length(doseData)==0, NA, paste(unique(doseData), collapse=", "))
      }
      
      # Description
      pddf  <- rbind(pddf, data.frame(a=doseAmount, b=length(idd)))
      nca_calc <- function(id,data) {
        if (is.null(doseAmtNm)){
          idzAmt <- NA
        }else{
          doseData <- as.numeric(as.character(data[data[,idNm]==id, doseAmtNm]))
          doseData <- doseData[complete.cases(doseData) & doseData>0]
          idzAmt   <- ifelse(length(doseData)==0, NA, doseData[1])
        }
        tcTI <- nca.ind.data(pkData=data, ID=id, dvLog = dvLog, dataType=dataType,
                             idNm=idNm, timeNm=timeNm, concNm=concNm,
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
                           LambdaExclude=LambdaExclude,doseTime=doseTime,Tau=Tau,TI=iTI,simFile=simFile,dset=dataType,onlyNCA=onlyNCA,extrapolate=extrapolate)
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
    cnm         <- c(Dcol,"No. of individuals")
    names(pddf) <- cnm
  }
  
  # Estimate NCA metrics for case = 2
  if (case == 2){
    for (s1 in 1:npopStr1){
      ifdf  <- pkData[pkData[,popStrNm1]==popStr1[s1],]
      if (nrow(ifdf) == 0){next}
      idd <- unique(as.character(ifdf[,idNm]))
      if (!is.null(doseAmtNm)){
        doseData   <- as.numeric(as.character(all_data[all_data[,popStrNm1]==popStr1[s1], doseAmtNm]))
        doseData   <- doseData[complete.cases(doseData) & doseData>0]
        doseAmount <- ifelse(length(doseData)==0, NA, paste(unique(doseData), collapse=", "))
      }else{
        doseAmount <- NA
      }
      # Description
      pddf <- rbind(pddf, data.frame(a=popStr1[s1], b=doseAmount, c=length(idd)))
      for (i in 1:length(idd)){
        if (!is.null(doseAmtNm)){
          doseData <- as.numeric(as.character(all_data[all_data[,popStrNm1]==popStr1[s1] & all_data[,idNm]==idd[i], doseAmtNm]))
          doseData <- doseData[complete.cases(doseData) & doseData>0]
          idzAmt   <- ifelse(length(doseData)==0, NA, doseData[1])
        }else{
          idzAmt <- NA
        }
        tcTI <- nca.ind.data(pkData=ifdf, ID=idd[i], dvLog = dvLog, dataType=dataType,
                             idNm=idNm, timeNm=timeNm, concNm=concNm,
                             adminType=adminType, TI=TI,
                             dateColNm=dateColNm, dateFormat=dateFormat, timeFormat=timeFormat)
        tc      <- tcTI$tc
        iTI     <- tcTI$iTI
        if (nrow(tc)==0) next
        time    <- as.numeric(tc$time)
        conc    <- as.numeric(tc$conc)
        cdata   <- rbind(cdata,cbind(Time=time,Conc=conc,ID=idd[i],STRAT1=popStr1[s1]))
        NCAprm  <- est.nca(time=time,conc=conc,backExtrp=backExtrp,negConcExcl=negConcExcl,
                           doseType=doseType,adminType=adminType,doseAmt=idzAmt,method=method,
                           AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,LambdaExclude=LambdaExclude,
                           doseTime=doseTime,Tau=Tau,TI=iTI,simFile=simFile,dset=dataType,onlyNCA=onlyNCA,extrapolate=extrapolate)
        outData <- rbind(outData, data.frame(ID=idd[i],STRAT1=popStr1[s1],Dose=idzAmt,t(NCAprm)))
      }
    }
    cnm         <- c(popStrNm1,Dcol,"No. of individuals")
    names(pddf) <- cnm
  }
  
  # Estimate NCA metrics for case = 3
  if (case == 3){
    for (s1 in 1:npopStr1){
      for (s2 in 1:npopStr2){
        ifdf <- pkData[pkData[,popStrNm1]==popStr1[s1] & pkData[,popStrNm2]==popStr2[s2],]
        if (nrow(ifdf) == 0){next}
        idd <- unique(as.character(ifdf[,idNm]))
        if (!is.null(doseAmtNm)){
          doseData   <- as.numeric(as.character(all_data[all_data[,popStrNm1]==popStr1[s1] & all_data[,popStrNm2]==popStr2[s2], doseAmtNm]))
          doseData   <- doseData[complete.cases(doseData) & doseData>0]
          doseAmount <- ifelse(length(doseData)==0, NA, paste(unique(doseData), collapse=", "))
        }else{
          doseAmount <- NA
        }
        # Description
        pddf   <- rbind(pddf, data.frame(a=popStr1[s1], b=popStr2[s2], c=doseAmount, d=length(idd)))
        for (i in 1:length(idd)){
          if (!is.null(doseAmtNm)){
            doseData <- as.numeric(as.character(all_data[all_data[,popStrNm1]==popStr1[s1] & all_data[,popStrNm2]==popStr2[s2] & all_data[,idNm]==idd[i], doseAmtNm]))
            doseData <- doseData[complete.cases(doseData) & doseData>0]
            idzAmt   <- ifelse(length(doseData)==0, NA, doseData[1])
          }else{
            idzAmt <- NA
          }
          tcTI <- nca.ind.data(pkData=ifdf, ID=idd[i], dvLog = dvLog, dataType=dataType,
                               idNm=idNm, timeNm=timeNm, concNm=concNm,
                               adminType=adminType, TI=TI,
                               dateColNm=dateColNm, dateFormat=dateFormat, timeFormat=timeFormat)
          tc      <- tcTI$tc
          iTI     <- tcTI$iTI
          if (nrow(tc)==0) next
          time    <- as.numeric(tc$time)
          conc    <- as.numeric(tc$conc)
          cdata   <- rbind(cdata,cbind(Time=time,Conc=conc,ID=idd[i],STRAT1=popStr1[s1],STRAT2=popStr2[s2]))
          NCAprm  <- est.nca(time=time,conc=conc,backExtrp=backExtrp,negConcExcl=negConcExcl,
                             doseType=doseType,adminType=adminType,doseAmt=idzAmt,method=method,
                             AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,LambdaExclude=LambdaExclude,
                             doseTime=doseTime,Tau=Tau,TI=iTI,simFile=simFile,dset=dataType,onlyNCA=onlyNCA,extrapolate=extrapolate)
          outData <- rbind(outData, data.frame(ID=idd[i],STRAT1=popStr1[s1],STRAT2=popStr2[s2],Dose=idzAmt,t(NCAprm)))
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
          ifdf <- pkData[(pkData[,popStrNm1]==popStr1[s1] & pkData[,popStrNm2]==popStr2[s2] & pkData[,popStrNm3]==popStr3[s3]),]
          if (nrow(ifdf) == 0){next}
          idd <- unique(as.character(ifdf[,idNm]))
          if (!is.null(doseAmtNm)){
            doseData   <- as.numeric(as.character(all_data[all_data[,popStrNm1]==popStr1[s1] & all_data[,popStrNm2]==popStr2[s2] & all_data[,popStrNm3]==popStr3[s3], doseAmtNm]))
            doseData   <- doseData[complete.cases(doseData) & doseData>0]
            doseAmount <- ifelse(length(doseData)==0, NA, paste(unique(doseData), collapse=", "))
          }else{
            doseAmount <- NA
          }
          # Description
          pddf  <- rbind(pddf, data.frame(a=popStr1[s1], b=popStr2[s2], c=popStr3[s3], d=doseAmount, e=length(idd)))
          for (i in 1:length(idd)){
            if (!is.null(doseAmtNm)){
              doseData <- as.numeric(as.character(all_data[all_data[,popStrNm1]==popStr1[s1] & all_data[,popStrNm2]==popStr2[s2] & all_data[,popStrNm3]==popStr3[s3] & all_data[,idNm]==idd[i], doseAmtNm]))
              doseData <- doseData[complete.cases(doseData) & doseData>0]
              idzAmt   <- ifelse(length(doseData)==0, NA, doseData[1])
            }else{
              idzAmt <- NA
            }
            tcTI <- nca.ind.data(pkData=ifdf, ID=idd[i], dvLog = dvLog, dataType=dataType,
                                 idNm=idNm, timeNm=timeNm, concNm=concNm,
                                 adminType=adminType, TI=TI,
                                 dateColNm=dateColNm, dateFormat=dateFormat, timeFormat=timeFormat)
            tc      <- tcTI$tc
            iTI     <- tcTI$iTI
            if (nrow(tc)==0) next
            time    <- as.numeric(tc$time)
            conc    <- as.numeric(tc$conc)
            cdata   <- rbind(cdata,cbind(Time=time,Conc=conc,ID=idd[i],STRAT1=popStr1[s1],STRAT2=popStr2[s2],STRAT3=popStr3[s3]))
            NCAprm  <- est.nca(time=time,conc=conc,backExtrp=backExtrp,negConcExcl=negConcExcl,
                               doseType=doseType,adminType=adminType,doseAmt=idzAmt,method=method,
                               AUCTimeRange=AUCTimeRange,LambdaTimeRange=LambdaTimeRange,LambdaExclude=LambdaExclude,
                               doseTime=doseTime,Tau=Tau,TI=iTI,simFile=simFile,dset=dataType,onlyNCA=onlyNCA,extrapolate=extrapolate)
            outData <- rbind(outData, data.frame(ID=idd[i],STRAT1=popStr1[s1],STRAT2=popStr2[s2],STRAT3=popStr3[s3],Dose=idzAmt,t(NCAprm)))
          }
        }
      }
    }
    cnm         <- c(popStrNm1,popStrNm2,popStrNm3,Dcol,"No. of individuals")
    names(pddf) <- cnm
  }
  
  # tests
  # browser()
  # result$pddf
  # pddf
  # dplyr::all_equal(result$pddf,pddf,convert = T)
  # 
  # tibble::as_tibble(result$outData)
  # tibble::as_tibble(outData)
  # 
  # result$outData$ID <- as.factor(result$outData$ID)
  # result$outData$STRAT1 <- as.factor(result$outData$STRAT1)
  # result$outData$STRAT2 <- as.factor(result$outData$STRAT2)
  # result$outData$STRAT3 <- as.factor(result$outData$STRAT3)
  # 
  # dplyr::all_equal(result$outData,outData,convert = T)
  # 
  # new_cdata <- result$cdata %>% purrr::map_df(as.factor)
  # new_cdata_2 <- cdata %>% purrr::map_df(as.factor)
  # 
  # dplyr::all_equal(new_cdata,new_cdata_2,convert = T)

  return(list(outData=outData,pddf=pddf,cdata=cdata))
}
