## ----setup, include=FALSE------------------------------------------------
library(knitr)
library(ncappc)
library(tibble)
library(dplyr)
opts_chunk$set(echo = FALSE, comment=NA, tidy = TRUE, 
               fig.align = 'center',
               fig.show='asis', out.width='90%', fig.pos='H')

## ----echo=TRUE,eval=FALSE------------------------------------------------
#  ncappc(obsFile = "file_with_observations.txt")

## ----out.height='90%',out.width='90%',fig.align='center',fig.cap="**Figure 1:** Schematic work-flow of the ```ncappc()``` function in the ***ncappc*** package.",echo = FALSE----
knitr::include_graphics('./ncappc-Schema.png')

## ----echo=TRUE,tidy=FALSE------------------------------------------------
data_1 <- tibble(
  ID=1,
  TIME = c(0,0.25,0.5,1,1.5,2,3,4,6,8,12,16,24),
  DV=c(0, 0.07, 0.14, 0.21, 0.24, 0.27, 0.26, 0.25, 0.22, 0.19, 0.13, 0.081, 0.033)
)
  
  
out <- ncappc(obsFile=data_1,
              onlyNCA = T,
              extrapolate = T,
              printOut = F,
              evid = FALSE,
              noPlot = T)

out$ncaOutput %>% select(c(AUClast,AUCINF_pred,Cmax,Tmax))

## ----echo=TRUE,eval=FALSE------------------------------------------------
#  ncappc(concUnit="ng/ml", doseAmtNm = "DOSE")

