# Plot DV vs Time data

# roxygen comments
#' Plots drug plasma concentration vs time data
#'
#' \pkg{dv.plot} plots DV vs Time data.
#' 
#' \pkg{dv.plot} plots DV vs Time data
#' 
#' @param df A data frame to be used for the plot
#' @param xvar is the independent variable, default is \strong{"TIME"}
#' @param yvar is the dependent variable, default is \strong{"DV"}
#' @param obsLog is a logical variable (\code{TRUE}, \code{FALSE}). If
#'   \code{TRUE}, concentration in observed data is assumed to be in logarithmic
#'   scale. Default is \strong{\code{FALSE}}
#' @param myXlab is the x-axis label, default is \strong{"Time"}
#' @param myYlab is the y-axis label, defaults is \strong{"Concentration"}
#' @param color is the column name of the color stratification variable, e.g. 
#'   "DOSEF". Default is \strong{NULL}
#' @param group is the column name of the variable used to group data, default
#'   is \strong{"ID}
#' @param guide if \strong{TRUE}, show guide, default is \strong{TRUE}
#' @param onlyLin if \strong{TRUE}, presents only the linear version of the
#'   plot, default is \strong{FALSE}
#' @param onlyLog if \strong{TRUE}, presents only the log version of the plot,
#'   default is \strong{FALSE}
#' @param XYlog if \strong{TRUE}, both X and Y axes of the log version of the
#'   plot is shown on the logarithmic scale; if \strong{FALSE}, only the Y-axis
#'   is shown on the logarithmic scale. Default is \strong{FALSE}.
#' @param STRATY is the row stratification variable, default is \strong{"."}
#' @param STRATX is the column stratification variable, default is \strong{"."}
#' @param myYBr are the breaks for the Y-axis for the linear plot
#' @param myXBr are the breaks for the X-axis for the linear plot
#' @param myYBrLog are the breaks for the Y-axis for the log plot
#' @param myXBrLog are the breaks for the X-axis for the log plot
#' @param myYlim sets Y-axis limits for the linear plot
#' @param myXlim sets X-axis limits for the linear plot
#' @param myYlimLog sets the Y-axis limit for the log plot
#' @param myXlimLog sets the X-axis limit for the log plot
#' 
#' @return returns a graphical object created by arrangeGrob function
#' @export
#'

dv.plot <- function(df,
                    xvar = "Time",
                    yvar = "Conc",
                    obsLog = FALSE,
                    myXlab = "Time",
                    myYlab = "Concentration",
                    color = NULL,
                    group = NULL,
                    guide = TRUE,
                    onlyLin = FALSE,
                    onlyLog = FALSE,
                    XYlog = FALSE,
                    STRATY = ".",
                    STRATX = ".",
                    myYBr = waiver(),
                    myXBr = waiver(),
                    myYBrLog = waiver(),
                    myXBrLog = waiver(),
                    myYlim = NULL,
                    myXlim = NULL,
                    myYlimLog = NULL,
                    myXlimLog = NULL){
  
  "ID" <- "Time" <- "Conc" <- "theme" <- "unit" <- "element_text" <- "xlab" <- "ylab" <- "geom_line" <- "aes_string" <- "geom_point" <- "ggplot" <- "facet_wrap" <- "scale_y_log10" <- "arrangeGrob" <- "textGrob" <- "gpar" <- "packageVersion" <- NULL
  rm(list=c("ID","Time","Conc","theme","unit","element_text","xlab","ylab","geom_line","aes_string","geom_point","ggplot","facet_wrap","scale_y_log10","arrangeGrob","textGrob","gpar","packageVersion"))
  
  df[,xvar] <- as.numeric(as.character(df[,xvar]))
  df[,yvar] <- as.numeric(as.character(df[,yvar]))
  if(obsLog) df[,yvar] <- exp(df[,yvar])
  
  if(is.null(color)){
    if(is.null(group)){
      p01 <- ggplot(df,aes_string(x=xvar,y=yvar,group=1))
    }else{
      p01 <- ggplot(df,aes_string(x=xvar,y=yvar,group=group))
    }
  }else{
    color <- paste0("factor(",color,")")
    if(is.null(group)){
      p01 <- ggplot(df,aes_string(x=xvar,y=yvar,color=color))
    }else{
      p01 <- ggplot(df,aes_string(x=xvar,y=yvar,group=group,color=color))
    }
  }
  
  p01 <- p01 + geom_line(alpha = 0.5)
  p01 <- p01 + geom_point()
  p01 <- p01 + xlab(myXlab)
  p01 <- p01 + ylab(myYlab)
  p01 <- p01 + theme(legend.position = "none",
                     axis.text.x = element_text(size=10),
                     axis.text.y = element_text(size=10),
                     strip.text.x = element_text(size=10),
                     title = element_text(size=10))
  
  facets <- paste(STRATY, '~',  STRATX)
  if (facets != '. ~ .') p01 <- p01 + facet_grid(facets, scales = "free")
  
  p02 <- p01 + scale_y_log10(breaks=myYBrLog,labels=myYBrLog)
  if(XYlog) p02 <- p02 + scale_x_log10(breaks=myXBrLog,labels=myXBrLog)
  
  p01 <- p01 + scale_x_continuous(breaks=myXBr, labels=myXBr)
  p01 <- p01 + scale_y_continuous(breaks=myYBr, labels=myYBr)
  
  if (!is.null(myYlim) & is.null(myXlim)){
    p01 <- p01 + coord_cartesian(ylim=myYlim)
  }else if (is.null(myYlim) & !is.null(myXlim)){
    p01 <- p01 + coord_cartesian(xlim=myXlim)
  }else if (!is.null(myYlim) & !is.null(myXlim)){
    p01 <- p01 + coord_cartesian(xlim=myXlim, ylim=myYlim)
  }
  
  if (!is.null(myYlimLog) & is.null(myXlimLog)){
    p02 <- p02 + coord_cartesian(ylim=myYlimLog)
  }else if (is.null(myYlimLog) & !is.null(myXlimLog)){
    p02 <- p02 + coord_cartesian(xlim=myXlimLog)
  }else if (!is.null(myYlimLog) & !is.null(myXlimLog)){
    p02 <- p02 + coord_cartesian(xlim=myXlimLog, ylim=myYlimLog)
  }
  
  if(onlyLin){
    p01 <- p01 + ggtitle("Concentration vs. Time profile\n")
    return(p01)
  }
  
  if(onlyLog){
    p02 <- p02 + ggtitle("Concentration vs. Time profile\n")
    return(p02)
  }
  
  if (!onlyLin & !onlyLog){
    p01 <- p01 %+% xlab("") %+% ylab("")
    p02 <- p02 %+% xlab("") %+% ylab("")
    
    plot_args <- list(p01,p02,ncol=2,
                      top=textGrob("Concentration vs. Time profile\n",vjust=1,gp=gpar(cex=0.8,fontface="bold")),
                      left=textGrob(myYlab,gp=gpar(cex=1,fontface="bold"),rot=90),
                      bottom=textGrob(myXlab,gp=gpar(cex=1,fontface="bold")))
    
    if(packageVersion("gridExtra") < "0.9.2"){
      arg_names <- names(plot_args)
      arg_names <- sub("top","main",arg_names)
      arg_names <- sub("bottom","sub",arg_names)
      names(plot_args) <- arg_names
    }  
    
    gdr <- suppressMessages(suppressWarnings(do.call(arrangeGrob,plot_args)))
    return(gdr)
  }
}


