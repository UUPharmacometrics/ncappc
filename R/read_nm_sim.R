#' Read nonmem table files produced from simulation.
#' 
#' The function reads in nonmem table files produced from the \code{$SIM} line 
#' in a NONMEM model file.
#' 
#' Currently the function expects the \code{$TABLE} to have a header for each 
#' new simulation. This means that the \code{NOHEADER} option or 
#' \code{ONEHEADER} option in the table file is not allowed.
#' 
#' @param table_sim The simulated table file to read. A text string.
#' @param only_obs Should the non-observation lines in the data set be removed? 
#'   Currently filtered uisng the expected \code{MDV} column. \code{TRUE} or 
#'   \code{FALSE}.
#'   
#' @return Returns a data frame of the simulated table with an added column for 
#'   the simulation number. The data frame is given class \code{c("tbl_df", 
#'   "tbl", "data.frame")} for easy use with \code{\link[dplyr]{dplyr}}.
#'
#' 

read_nm_sim <- function (table_sim, only_obs=TRUE) {
  
  ## change header names 
  sim <- readr::read_table(table_sim, col_names = TRUE, skip = 1)
  names(sim) <- gsub("\\n.*","",names(sim))
  
  ## create simulation number
  args <- lazyeval::interp(~ cumsum(is.na(var))+1, var = as.name(names(sim)[1]))
  sim <- dplyr::mutate_(sim,NSUB=args)

  ## filter out NA columns
  args <- lazyeval::interp(~ !is.na(var), var = as.name(names(sim)[1]))
  sim <- dplyr::filter_(sim,args)
  
  ## remove non-observation rows
  if(only_obs){
    if(any("MDV"==names(sim))){
      sim <- dplyr::filter(sim,MDV==0)   
    } else {
      warning('\nWarning MDV data item not listed in header, 
              Could not remove dose events!\n')
    }
  }
  #summarise_each(sim,funs(mean,sd))
  return(sim)
}
