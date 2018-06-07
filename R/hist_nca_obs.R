

hist_nca_obs <- function(case, outData, AUClast, AUCINF_obs, Cmax, Tmax, cunit, tunit, spread, 
                         printOut, usrdir, figFormat, histobsplot, 
                         npopStr1, STRAT1, popStr1, popStrNm1, 
                         npopStr2, STRAT2, popStr2, popStrNm2, 
                         npopStr3, STRAT3, popStr3, popStrNm3) {
  
  
  if (case == 1){
    # Obs hist plot
    if (nrow(outData)>=5){
      plotData <- subset(outData, select=c(AUClast,AUCINF_obs,Cmax,Tmax))
      numPrm   <- sapply(plotData, FUN=function(x){x <- as.numeric(as.character(x)); length(x[complete.cases(x)])})
      if (length(numPrm[numPrm>=5]) == 0) return(NULL)
      
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
  
  if (case == 2){
    for (s1 in 1:npopStr1){
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
  if (case == 3){
    for (s1 in 1:npopStr1){
      for (s2 in 1:npopStr2){
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
  
  if (case == 4){
    for (s1 in 1:npopStr1){
      for (s2 in 1:npopStr2){
        for (s3 in 1:npopStr3){
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
  return(histobsplot)


}