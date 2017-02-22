# Read NONMEM simulation output file
#
# roxygen comments
#' Check observed data
#'
#' \pkg{nca.read.sim} Reads NONMEM simulation output file.
#' 
#' \pkg{nca.read.sim} Reads NONMEM simulation output file.
#' 
#' @param simFile NONMEM simulation output with the simulated
#'   concentration-time data from an internal data frame or an external table.
#'   Default is \strong{"nca_simulation.1.npctab.dta"}.
#' @param MDV.rm If \code{TRUE} MDV column is used to remove non-observation
#'   rows. Default is \strong{\code{FALSE}}
#' 
#' @return A list of objects
#' @export
#'

# read NONMEM output into individual simulation data file
nca.read.sim <- function(simFile="nca_simulation.1.npctab.dta",
                         MDV.rm=FALSE){
  sim <- read.table(simFile, skip=1, header=T, fill=T, as.is=T)
  sim <- as.data.frame(apply(sim,2,function(x) suppressWarnings(as.numeric(x))))
  Nro <- min(which(is.na(sim[,1])))-1
  sim <- sim[!is.na(sim[,1]),]
  sim$NSUB <- rep(1:(nrow(sim)/Nro),each=Nro)
  if(MDV.rm){
    if(any(colnames(sim)=='MDV')){sim <- sim[sim[,'MDV']==0,]
    }else{cat('\nWarning MDV data item not listed in header,
                         Could not remove dose events!')}
  }
  assign("nmdf", sim)
  return(nmdf)
}