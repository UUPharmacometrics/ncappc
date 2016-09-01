# Population histogram

# roxygen comments
#' Plots population histogram of the NCA metrics selected for model diagnosis.
#'
#' \pkg{histpop.plot} plots population histogram of the NCA metrics selected 
#' for model diagnosis (e.g. AUClast, AUCINF_obs, Cmax and Tmax).
#'
#' \pkg{histpop.plot} plots histogram of the NCA metrics selected for the model 
#' diagnosis and compares with the corresponding metrics estimated from the 
#' observed data. The allowed NCA metrics for this histograms are "AUClast", 
#' "AUClower_upper", "AUCINF_obs", "AUCINF_pred", "AUMClast", "Cmax", "Tmax" and
#' "HL_Lambda_z". By default, this function produces histogram of AUClast and 
#' Cmax.
#' 
#' @param obsdata Data frame with the values of the NCA metrics estimated from
#'   the observed data
#' @param simdata Data frame with the values of the NCA metrics estimated from
#'   the simulated data
#' @param figlbl Figure label based on dose identifier and/or population
#'   stratifier (\strong{NULL})
#' @param param A character array of the NCA metrics. The allowed NCA metrics 
#'   for this histograms are "AUClast", "AUClower_upper", "AUCINF_obs", 
#'   "AUCINF_pred", "AUMClast", "Cmax", "Tmax" and "HL_Lambda_z". 
#'   (\strong{c("AUClast", "Cmax")})
#' @param cunit Unit for concentration (\strong{"M.L^-3"})
#' @param tunit Unit for time (\strong{"T"})
#' @param spread Measure of the spread of simulated data (ppi (95\% parametric
#'   prediction interval) or npi (95\% nonparametric prediction interval))
#'   (\strong{"npi"})
#'
#' @return returns a graphical object created by arrangeGrob function
#' @export
#'

histpop.plot <- function(obsdata=outData,
                         simdata=smedianData,
                         figlbl=NULL,
                         param=c("AUClast","Cmax"),
                         cunit="M.L^-3",
                         tunit="T",
                         spread="npi"){
  
  "..density.." <- "TYPE" <- "obs" <- "sim" <- "arrangeGrob" <- "scale_linetype_manual" <- "scale_color_manual" <- "xlab" <- "ylab" <- "guides" <- "guide_legend" <- "theme" <- "element_text" <- "unit" <- "element_rect" <- "geom_histogram" <- "aes" <- "geom_vline" <- "melt" <- "ggplot" <- "labs" <- "coord_cartesian" <- "facet_wrap" <- "gtable_filter" <- "ggplot_gtable" <- "ggplot_build" <- "textGrob" <- "gpar" <- "..count.." <- "..PANEL.." <- "scale_y_continuous" <- "percent" <- "sd" <- "quantile" <- "packageVersion" <- NULL
  rm(list=c("..density..","TYPE","obs","sim","arrangeGrob","scale_linetype_manual","scale_color_manual","xlab","ylab","guides","guide_legend","theme","element_text","unit","element_rect","geom_histogram","aes","geom_vline","melt","ggplot","labs","coord_cartesian","facet_wrap","gtable_filter","ggplot_gtable","ggplot_build","textGrob","gpar","..count..","..PANEL..","scale_y_continuous","percent","sd","quantile","packageVersion"))
  
  outData <- obsdata; smedianData <- simdata
  
  alwprm <- c("AUClast","AUClower_upper","AUCINF_obs","AUCINF_pred","AUMClast","Cmax","Tmax","HL_Lambda_z")
  npr    <- length(param)
  fctNm  <- data.frame()
  nc <- ifelse(npr<2, 1, ifelse(npr>=2 & npr<=6, 2, 3))
  if (!all(param%in%alwprm)){setwd("..");stop("Incorrect NCA metrics. Please select NCA metrics from \"AUClast\", \"AUClower_upper\", \"AUCINF_obs\", \"AUCINF_pred\", \"AUMClast\", \"Cmax\", \"Tmax\", \"HL_Lambda_z\".")}
  
  # ggplot variables
  ggOpt_pop <- list(scale_linetype_manual(name="",values=c("median(obs)"="solid","median(medianSim)"="solid","+/-spread"="dashed")),
                    scale_color_manual(name = "", values=c("median(obs)"="red","median(medianSim)"="blue","+/-spread"="blue")),
                    xlab(""), ylab(""),
                    guides(fill = guide_legend(override.aes = list(linetype = 0 )), shape = guide_legend(override.aes = list(linetype = 0))),
                    theme(axis.text.x  = element_text(angle=45,vjust=1,hjust=1),
                          axis.text.y  = element_text(hjust=0),
                          legend.position = "bottom", legend.direction = "horizontal",
                          legend.background = element_rect()),
                    geom_vline(aes(xintercept=as.numeric(obs), color="median(obs)", linetype="median(obs)"), size=1, show.legend=T),
                    geom_vline(aes(xintercept=as.numeric(median), color="median(medianSim)", linetype="median(medianSim)"), size=1),
                    geom_vline(aes(xintercept=as.numeric(sprlow), color="+/-spread", linetype="+/-spread"), size=1),
                    geom_vline(aes(xintercept=as.numeric(sprhgh), color="+/-spread", linetype="+/-spread"), size=1),
                    scale_y_continuous(labels = percent))
  
  obsVal       <- sapply(obsdata, FUN=function(x) median(as.numeric(x), na.rm=T))
  medianMedian <- sapply(simdata, FUN=function(x) median(as.numeric(x), na.rm=T))
  sdMedian     <- sapply(simdata, FUN=function(x) sd(as.numeric(x), na.rm=T))
  #obsVal   <- sapply(obsdata, FUN=function(x) mean(as.numeric(x), na.rm=T))
  meanMean <- sapply(simdata, FUN=function(x) mean(as.numeric(x), na.rm=T))
  sdMean   <- sapply(simdata, FUN=function(x) sd(as.numeric(x), na.rm=T))
  xlow     <- sapply(simdata, FUN=function(x) unname(quantile(as.numeric(x),0.01,na.rm=T)))
  xhgh     <- sapply(simdata, FUN=function(x) unname(quantile(as.numeric(x),0.99,na.rm=T)))
  if (spread=="ppi"){
    sprlow <- meanMean-1.96*sdMean
    sprhgh <- meanMean+1.96*sdMean
  }else if (spread=="npi"){
    sprlow <- sapply(simdata, FUN=function(x) unname(quantile(as.numeric(x),0.025,na.rm=T)))
    sprhgh <- sapply(simdata, FUN=function(x) unname(quantile(as.numeric(x),0.975,na.rm=T)))
  }
  
  longData <- melt(simdata,measure=param)
  names(longData) <- c("TYPE","sim")
  longData <- cbind(longData,median=0,mean=0,sd=0,sprlow=0,sprhgh=0,obs=0,xlow=0,xhgh=0)
  
  for (p in 1:npr){
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
    longData[longData$TYPE==param[p],"median"] <- medianMedian[param[p]]
    longData[longData$TYPE==param[p],"mean"]   <- meanMean[param[p]]
    longData[longData$TYPE==param[p],"sd"]     <- sdMean[param[p]]
    longData[longData$TYPE==param[p],"sprlow"] <- sprlow[param[p]]
    longData[longData$TYPE==param[p],"sprhgh"] <- sprhgh[param[p]]
    longData[longData$TYPE==param[p],"obs"]    <- obsVal[param[p]]
    longData[longData$TYPE==param[p],"xlow"]   <- min(xlow[param[p]],sprlow[param[p]],obsVal[param[p]])
    longData[longData$TYPE==param[p],"xhgh"]   <- max(xhgh[param[p]],sprhgh[param[p]],obsVal[param[p]])
  }
  
  devtag <- ifelse (spread=="ppi","95% parametric prediction interval","95% nonparametric prediction interval")
  gplt <- list()
  for (p in 1:npr){
    df <- subset(longData, TYPE==param[p])
    df$TYPE <- factor(df$TYPE, levels=param[p], labels=fctNm[fctNm$prmNm==param[p],"prmUnit"])
    df$FCT  <- paste0(df$TYPE,"\nmedian(obs)=",out.digits(df$obs[1],dig=4),"\nmedian(medianSim)=",out.digits(df$median[1],dig=4),"\n+/-spread=(",out.digits(df$sprlow[1],dig=4),",",out.digits(df$sprhgh[1],dig=4),")")
    xl <- df$xlow[1]; xu <- df$xhgh[1]
    bw <- diff(unname(quantile(as.numeric(df$sim),c(0.005,0.985))))/(2*IQR(as.numeric(df$sim)))/length(as.numeric(df$sim))^(1/3)
    gplt[[p]] <- ggplot(df,aes(x=as.numeric(sim))) +
      geom_histogram(aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]), size=0.6, color="black", fill="white", binwidth = bw) +
      ggOpt_pop + coord_cartesian(xlim=c(xl,xu)) + facet_wrap(~FCT, scales="free")
  }
  mylegend <- suppressMessages(suppressWarnings(gtable_filter(ggplot_gtable(ggplot_build(gplt[[1]])), "guide-box", trim=T)))
  lheight  <- sum(mylegend$heights)
  for (p in 1:npr){gplt[[p]] <- gplt[[p]] + theme(legend.position="none")}
  
  if(is.null(figlbl)){
    Label <- paste("Histogram of simulated population medians\n(spread = ",devtag,")\n\n",sep="")
  }else{
    Label <- paste("Histogram of simulated population medians (",figlbl,")\n(spread = ",devtag,")\n\n",sep="")
  }
  
  plot_args <- list(top = textGrob(Label,vjust=1,gp=gpar(cex = 1.5)),
                    bottom = textGrob("Value\n\n",vjust=1,gp=gpar(cex = 1.5)),
                    ncol=nc)
  if(packageVersion("gridExtra") < "0.9.2"){
    arg_names <- names(plot_args)
    arg_names <- sub("top","main",arg_names)
    arg_names <- sub("bottom","sub",arg_names)
    names(plot_args) <- arg_names
  }  
  gdr <- suppressMessages(suppressWarnings(do.call(arrangeGrob,c(gplt,plot_args))))
  histpopgrob <- list(gdr=gdr,legend=mylegend,lheight=lheight)
  return(histpopgrob)
}

