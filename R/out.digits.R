#' @title output value with correct digits and trailing zero
#'
#' @description Function to present a value with correct digits and trailing zero
#'
#' @details This is a function to present a value with correct digits and trailing zero. Numbers >= 10000, or <= 0.0001 will be presented in scientific format
#' @name out.digits
#' @param x is the value
#' @param dig is the number of significant digits
#' @export
#' @examples \dontrun{
#'out.digits(1234)
#' }
#'


## Format parameter values in output table
out.digits <- function(x, dig=3){
  if(length(x)==0){
    tmp <- ""
  }else if(length(x)>1){
    tmp <- NULL
    for(i in 1:length(x)){
      if(is.na(x[i])){
        tmp <- c(tmp,NA)
      }else if(!grepl("^[-]?[0-9]*[.]?[0-9]*[eE]?[-]?[0-9]*[.]?[0-9]*$", x[i])){
        tmp <- c(tmp,x[i])
      }else{
        y <- as.numeric(x[i])
        t <- ifelse(y==0, "0", ifelse(abs(round(y,1)) < 1*10^(dig-1), as.character(formatC(signif(y,dig), digits=dig,format="fg", flag="#")),
                                      ifelse(abs(round(y,1)) < 1*10^(dig+1), as.character(signif(y, dig)),
                                             as.character(formatC(y, digits=dig-1, format="E")))))
        tmp <- c(tmp,t)
      }
    }
  }else{
    if(is.na(x)){
      tmp <- NA
    }else if(!grepl("^[-]?[0-9]*[.]?[0-9]*[eE]?[-]?[0-9]*[.]?[0-9]*$", x)){
      tmp <- x
    }else{
      x <- as.numeric(x)
      tmp <- ifelse(x==0, "0", ifelse(abs(round(x,1)) < 1*10^(dig-1), as.character(formatC(signif(x,dig), digits=dig,format="fg", flag="#")),
                                      ifelse(abs(round(x,1)) < 1*10^(dig+1), as.character(signif(x, dig)),
                                             as.character(formatC(x, digits=dig-1, format="E")))))
    }
  }
  return(tmp)
}

