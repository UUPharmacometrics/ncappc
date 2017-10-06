

dv_vs_idv <- function(case, cdata, concplot, printOut, usrdir, figFormat,
                      cunit, tunit,
                      npopStr1, STRAT1, popStr1, popStrNm1, 
                      npopStr2, STRAT2, popStr2, popStrNm2, 
                      npopStr3, STRAT3, popStr3, popStrNm3) {
  
  if(is.null(cunit)){DVleg  <- "Concentration"}else{DVleg <- paste0("Concentration (",cunit,")")}    # DV legend
  if(is.null(tunit)){IDVleg <- "Time"}else{IDVleg <- paste0("Time (",tunit,")")}                     # IDV legend
  
  
  if (case == 1){
    figlbl <- "All-data"
    gdr    <- dv.plot(df=cdata,xvar="Time",yvar="Conc",myXlab=IDVleg,myYlab=DVleg,color="ID",title=figlbl)
    #suppressMessages(suppressWarnings(grid.arrange(gdr)))
    #ggr <- grid.grab()
    #concplot[[length(concplot)+1]] <- ggr
    concplot[[length(concplot)+1]] <- gdr
    if (printOut){
      fl <- paste0(usrdir,"/TimeConc_",figlbl,".",figFormat)
      if (figFormat=="tiff"){
        eval(parse(text=paste0(figFormat,"(file=\"",fl,"\",height=10,width=15,units=\"cm\",compression=\"lzw\",res=120)")))
      }else{
        eval(parse(text=paste0(figFormat,"(file=\"",fl,"\",height=10,width=15,units=\"cm\",res=120)")))
      }
      suppressMessages(suppressWarnings(grid.arrange(gdr)))
      grDevices::dev.off()
    }
  }
  
  if (case == 2){
    for (s1 in 1:npopStr1){
      plotData <- subset(cdata, STRAT1==popStr1[s1])
      figlbl <- paste0(popStrNm1,"-",popStr1[s1])
      gdr    <- dv.plot(df=plotData,xvar="Time",yvar="Conc",myXlab=IDVleg,myYlab=DVleg,color="ID",title=figlbl)
      #suppressMessages(suppressWarnings(grid.arrange(gdr)))
      #ggr <- grid.grab()
      #concplot[[length(concplot)+1]] <- ggr
      concplot[[length(concplot)+1]] <- gdr
      if (printOut){
        fl <- paste0(usrdir,"/TimeConc_",figlbl,".",figFormat)
        if (figFormat=="tiff"){
          eval(parse(text=paste0(figFormat,"(file=\"",fl,"\",height=10,width=15,units=\"cm\",compression=\"lzw\",res=120)")))
        }else{
          eval(parse(text=paste0(figFormat,"(file=\"",fl,"\",height=10,width=15,units=\"cm\",res=120)")))
        }
        suppressMessages(suppressWarnings(grid.arrange(gdr)))
        grDevices::dev.off()
      }
    }
  }  
  if (case == 3){
    for (s1 in 1:npopStr1){
      for (s2 in 1:npopStr2){
        plotData <- subset(cdata, STRAT1==popStr1[s1] & STRAT2==popStr2[s2])
        figlbl <- paste0(popStrNm1,"-",popStr1[s1],"_",popStrNm2,"-",popStr2[s2])
        gdr    <- dv.plot(df=plotData,xvar="Time",yvar="Conc",myXlab=IDVleg,myYlab=DVleg,color="ID",title=figlbl)
        #suppressMessages(suppressWarnings(grid.arrange(gdr)))
        #ggr <- grid.grab()
        #concplot[[length(concplot)+1]] <- ggr
        concplot[[length(concplot)+1]] <- gdr
        
        if (printOut){
          fl <- paste0(usrdir,"/TimeConc_",figlbl,".",figFormat)
          if (figFormat=="tiff"){
            eval(parse(text=paste0(figFormat,"(file=\"",fl,"\",height=10,width=15,units=\"cm\",compression=\"lzw\",res=120)")))
          }else{
            eval(parse(text=paste0(figFormat,"(file=\"",fl,"\",height=10,width=15,units=\"cm\",res=120)")))
          }
          suppressMessages(suppressWarnings(grid.arrange(gdr)))
          grDevices::dev.off()
        }
      }
    }
  }
  
  if (case == 4){
    for (s1 in 1:npopStr1){
      for (s2 in 1:npopStr2){
        for (s3 in 1:npopStr3){
          plotData <- subset(cdata, STRAT1==popStr1[s1] & STRAT2==popStr2[s2] & STRAT3==popStr3[s3])
          if(nrow(plotData)==0) next
          figlbl <- paste0(popStrNm1,"-",popStr1[s1],"_",popStrNm2,"-",popStr2[s2],"_",popStrNm3,"-",popStr3[s3])
          gdr    <- dv.plot(df=plotData,xvar="Time",yvar="Conc",myXlab=IDVleg,myYlab=DVleg,color="ID",title=figlbl)
          #suppressMessages(suppressWarnings(grid.arrange(gdr)))
          #ggr <- grid.grab()
          #concplot[[length(concplot)+1]] <- ggr
          concplot[[length(concplot)+1]] <- gdr
          if (printOut){
            fl <- paste0(usrdir,"/TimeConc_",figlbl,".",figFormat)
            if (figFormat=="tiff"){
              eval(parse(text=paste0(figFormat,"(file=\"",fl,"\",height=10,width=15,units=\"cm\",compression=\"lzw\",res=120)")))
            }else{
              eval(parse(text=paste0(figFormat,"(file=\"",fl,"\",height=10,width=15,units=\"cm\",res=120)")))
            }
            suppressMessages(suppressWarnings(grid.arrange(gdr)))
            grDevices::dev.off()
          }
        }
      }
    }
  }
  
  return(concplot)
  
}