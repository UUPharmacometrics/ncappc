# Plot DV vs Time data

# roxygen comments
#' Plots drug plasma concentration vs time data
#'
#' \pkg{dv.plot} plots DV vs Time data.
#'
#' \pkg{dv.plot} plots DV vs Time data
#' 
#' @param pdata A data frame with four mandatory columns:
#' \enumerate{
#'  \item Time: Time column
#'  \item Conc: DV column
#'  \item ID: Individual ID
#'  \item FCT: Dose identifier and/or population stratifier
#' }
#' @param cunit Unit for concentration (\strong{"[M].[L]^-3"})
#' @param tunit Unit for time (\strong{"[T]"})
#'
#' @return returns a graphical object created by arrangeGrob function
#' @export
#'

dv.plot <- function(pdata,cunit="[M].[L]^-3",tunit="[T]"){
  ID <- Time <- Conc <- NULL
  rm(list=c("ID","Time","Conc"))
  
  ggOpt_conc <- list(theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
                           panel.margin = unit(0.2, "cm"),
                           axis.text.x  = element_text(size=9,face="bold",color="black",vjust=1,hjust=0.5),
                           axis.text.y  = element_text(size=9,face="bold",color="black",hjust=0),
                           legend.position = "none",
                           strip.text.x = element_text(size=9, face="bold")),
                     xlab(""),ylab(""),geom_line(aes(colour=factor(ID)), size=0.5),geom_point(aes(colour=factor(ID)), size=1))
  p01 <- ggplot(pdata, aes(x=as.numeric(as.character(Time)), y=as.numeric(as.character(Conc)), color=factor(ID))) +
    ggOpt_conc + #geom_smooth(method="loess",size=1,color=rgb(139, 0, 0, maxColorValue=255),se=F) +
    facet_wrap(~FCT, scales="free", ncol=1)
  p02 <- p01 + scale_y_log10()
  gdr <- suppressMessages(suppressWarnings(
    arrangeGrob(p01,p02,ncol=2,
                main=textGrob("Concentration vs. Time profile\n",vjust=1,gp=gpar(fontface="bold",cex = 0.7)),
                left=textGrob(paste("\nPlasma concentration (",cunit,")",sep=""),gp=gpar(fontface="bold",cex = 0.7),rot=90),
                sub=textGrob(paste("Time (",tunit,")\n"),gp=gpar(fontface="bold",cex = 0.7)))))
  return(gdr)
}