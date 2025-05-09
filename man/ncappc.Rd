% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ncappc.R
\name{ncappc}
\alias{ncappc}
\title{Performs NCA calculations and population PK model diagnosis.}
\usage{
ncappc(
  obsFile = "nca_original.npctab.dta",
  simFile = "nca_simulation.1.npctab.dta.zip",
  str1Nm = NULL,
  str1 = NULL,
  str2Nm = NULL,
  str2 = NULL,
  str3Nm = NULL,
  str3 = NULL,
  concUnit = NULL,
  timeUnit = NULL,
  doseUnit = NULL,
  obsLog = FALSE,
  simLog = obsLog,
  psnOut = TRUE,
  idNmObs = "ID",
  timeNmObs = "TIME",
  concNmObs = "DV",
  idNmSim = idNmObs,
  timeNmSim = timeNmObs,
  concNmSim = concNmObs,
  onlyNCA = FALSE,
  AUCTimeRange = NULL,
  backExtrp = FALSE,
  LambdaTimeRange = NULL,
  LambdaExclude = NULL,
  doseAmtNm = NULL,
  adminType = "extravascular",
  doseType = "ns",
  doseTime = NULL,
  Tau = NULL,
  TI = NULL,
  method = "linearup-logdown",
  blqNm = NULL,
  blqExcl = 1,
  evid = TRUE,
  evidIncl = 0,
  mdv = FALSE,
  filterNm = NULL,
  filterExcl = NULL,
  negConcExcl = FALSE,
  param = c("AUClast", "Cmax"),
  timeFormat = "number",
  dateColNm = NULL,
  dateFormat = NULL,
  spread = "npi",
  tabCol = c("AUClast", "Cmax", "Tmax", "AUCINF_obs", "Vz_obs", "Cl_obs", "HL_Lambda_z"),
  figFormat = "tiff",
  noPlot = FALSE,
  printOut = TRUE,
  studyName = NULL,
  new_data_method = TRUE,
  overwrite_SIMDATA = NULL,
  overwrite_sim_est_file = NULL,
  outFileNm = NULL,
  out_format = "html",
  gg_theme = theme_bw(),
  parallel = FALSE,
  extrapolate = FALSE,
  timing = FALSE,
  ...
)
}
\arguments{
\item{obsFile}{Observed concentration-time data from an internal data frame
or an external table with comma, tab or space as separators.}

\item{simFile}{NONMEM simulation output with the simulated concentration-time
data from an internal data frame or an external table. \code{NULL} produces
just the NCA output, a filename or data frame produces the NCA output as
well as the PopPK diagnosis. If \code{new_data_method=TRUE} then this can
be a compressed file as well.}

\item{str1Nm}{Column name for 1st level population stratifier. Default is
\strong{\code{NULL}}}

\item{str1}{Stratification ID of the members within 1st level stratification
(e.g c(1,2)). Default is \strong{\code{NULL}}}

\item{str2Nm}{Column name for 2nd level population stratifier. Default is
\strong{\code{NULL}}}

\item{str2}{Stratification ID of the members within 2nd level stratification
(e.g c(1,2)). Default is \strong{\code{NULL}}}

\item{str3Nm}{Column name for 3rd level population stratifier. Default is
\strong{\code{NULL}}}

\item{str3}{Stratification ID of the members within 3rd level stratification
(e.g c(1,2)). Default is \strong{\code{NULL}}}

\item{concUnit}{Unit of concentration (e.g. "ng/mL"). Default is 
\strong{\code{NULL}}}

\item{timeUnit}{Unit of time (e.g. "h"). Default is \strong{\code{NULL}}}

\item{doseUnit}{Unit of dose amount (e.g. "ng"). Default is
\strong{\code{NULL}}}

\item{obsLog}{If \code{TRUE} concentration in observed data is in logarithmic
scale. Default is \strong{\code{FALSE}}}

\item{simLog}{If \code{TRUE} concentration in simulated data is in 
logarithmic scale. Default is \strong{\code{FALSE}}}

\item{psnOut}{If \code{TRUE} observed data is an output from PsN or in NONMEM
output format. Default is \strong{\code{TRUE}}}

\item{idNmObs}{Column name for ID in observed data. Default is \strong{"ID"}}

\item{timeNmObs}{Column name for time in observed data. Default is
\strong{"TIME"}}

\item{concNmObs}{Column name for concentration in observed data. Default is
\strong{"DV"}}

\item{idNmSim}{Column name for ID in simulated data. Default is \strong{"ID"}}

\item{timeNmSim}{Column name for time in simulated data. Default is
\strong{"TIME"}}

\item{concNmSim}{Column name for concentration in simulated data. Default is
\strong{"DV"}}

\item{onlyNCA}{If \code{TRUE} only NCA is performed and ppc part is ignored
although simFile is not \code{NULL}. Default is \strong{\code{FALSE}}}

\item{AUCTimeRange}{User-defined window of time used to estimate AUC. Default
is \strong{\code{NULL}}}

\item{backExtrp}{If \code{TRUE} back-extrapolation is performed while
estimating AUC. Default is \strong{\code{FALSE}}}

\item{LambdaTimeRange}{User-defined window of time to estimate elimination
rate-constant. This argument lets the user to choose a specific window of
time to be used to estimate the elimination rate constant (Lambda) in the
elimination phase. The accepted format for the input to this argument is a
numeric array of two elements; \code{c(14,24)} will estimate the Lambda
using the data within the time units 14 to 24. Default is
\strong{\code{NULL}}}

\item{LambdaExclude}{User-defined excluded observation time points for
estimation of Lambda. This can be numeric value or logical condition (e.g.
c(1, 2, "<20", ">=100", "!=100")). Default is \strong{\code{NULL}}}

\item{doseAmtNm}{Column name to specify dose amount. Default is
\strong{\code{NULL}}}

\item{adminType}{Route of administration. Allowed options are iv-bolus,
iv-infusion or extravascular. Default is \strong{"extravascular"}}

\item{doseType}{Steady-state (ss) or non-steady-state (ns) dose. Default is
\strong{"ns"}}

\item{doseTime}{Dose time prior to the first observation for steady-state
data. Default is \strong{\code{NULL}}}

\item{Tau}{Dosing interval for steady-state data. Default is
\strong{\code{NULL}}}

\item{TI}{Infusion duration. If TI is a single numeric value, TI is the same
for all individuals. If TI is the name of a column with numeric data
present in the data set, TI is set to the unique value of the column for a
given individual. Default is \strong{\code{NULL}}}

\item{method}{Method to estimate AUC. \code{linear} method applies the linear
trapezoidal rule to estimate the area under the curve. \code{log} method 
applies the logarithmic trapezoidal rule to estimate the area under the 
curve. \code{linearup-logdown} method applies the linear trapezoidal rule
to estimate the area under the curve for the ascending part of the curve
and the logarithmic trapezoidal rule to estimate the area under the curve
for the descending part of the curve. Default is
\strong{"linearup-logdown"}}

\item{blqNm}{Name of BLQ column if used to exclude data. Default is
\strong{\code{NULL}}}

\item{blqExcl}{Excluded BLQ value; either a numeric value or a logical
condition (e.g. 1 or ">=1" or c(1,">3")). Used only if the \code{blqNm} is
not \code{NULL}. Default is \strong{"1"}}

\item{evid}{If \code{TRUE} EVID is used to filter data. Default is
\strong{\code{TRUE}}}

\item{evidIncl}{Included values in EVID. Default is \strong{"0"}}

\item{mdv}{If \code{TRUE} MDV is used to include data when MDV=0. Default is
\strong{\code{FALSE}}}

\item{filterNm}{Column name to filter data. Default is \strong{\code{NULL}}}

\item{filterExcl}{Row exclusion criteria based on the column defined by
\code{filterNm}. This can be numeric value or logical condition (e.g. c(1,
2, "<20", ">=100", "!=100")). Default is \strong{\code{NULL}}}

\item{negConcExcl}{If \code{TRUE} negative concentrations are excluded.
Default is \strong{\code{FALSE}}}

\item{param}{NCA parameters (AUClast, AUClower_upper, AUCINF_obs,
AUCINF_pred, AUMClast, Cmax, Tmax, HL_Lambda_z). Default is
\strong{(c"AUClast", "Cmax")}}

\item{timeFormat}{time format (number, H:M, H:M:S). Default is
\strong{"number"}}

\item{dateColNm}{column name for date if used (e.g. "Date", "DATE"). Default
is \strong{\code{NULL}}}

\item{dateFormat}{date format (D-M-Y, D/M/Y or any other combination of
D,M,Y). Default is \strong{\code{NULL}}}

\item{spread}{Measure of the spread of simulated data (\code{"ppi"} (95\%
parametric prediction interval) or \code{"npi"} (95\% nonparametric
prediction interval)). Default is \strong{"npi"}}

\item{tabCol}{Output columns to be printed in the report in addition to ID,
dose and population strata information (list of NCA metrics in a string
array). Default is \strong{c("AUClast", "Cmax", "Tmax", "AUCINF_obs",
"Vz_obs", "Cl_obs", "HL_Lambda_z")}}

\item{figFormat}{format of the produced figures (bmp, jpeg, tiff, png).
Default is \strong{"tiff"}}

\item{noPlot}{If \code{TRUE} only NCA calculations are performed without any
plot generation. Default is \strong{\code{FALSE}}}

\item{printOut}{If \code{TRUE} tabular and graphical outputs are saved on the
disk. Default is \strong{\code{TRUE}}}

\item{studyName}{Name of the study to be added as a description in the
report. Default is \strong{\code{NULL}}}

\item{new_data_method}{If \code{TRUE} a faster method of reading data is
tested. Default is \strong{\code{TRUE}}}

\item{overwrite_SIMDATA}{If \code{TRUE} new information is created in the
SIMDATA directory. If \code{FALSE} the information in the SIMDATA directory
is used. If \code{NULL} a dialog will come up to ask the user what to do.
Default is \strong{\code{NULL}}}

\item{overwrite_sim_est_file}{If \code{TRUE} The NCA metrics are created again based
on the simulation data.  If \code{FALSE} the information in the ncaSimEst file
  is used. If \code{NULL} a dialog will come up to ask the user what to do.
  Default is \strong{\code{NULL}}}

\item{outFileNm}{Additional tag to the name of the output html and pdf output
file hyphenated to the standard ncappc report file name standard ncappc
report file name. Default is \strong{\code{NULL}}}

\item{out_format}{What type of output format should the NCA report have? 
Pass "all" to render all formats defined within the rmarkdown file.
Pass "first" to render the first format defined within the rmarkdown file.
Pass "html" to render in HTML.
Pass "pdf" to render in PDF.}

\item{gg_theme}{Which ggplot theme should be used for the plots?}

\item{parallel}{Should the nca computations for the simulated data be run in parallel? See
\code{\link[PopED]{start_parallel}} for a description and additional arguments that can be 
added to this function and passed to \code{\link[PopED]{start_parallel}}.}

\item{extrapolate}{Should the NCA calculations extrapolate from the last observation to infinity?}

\item{timing}{Should timings of calculations be reported to the screen?}

\item{...}{Additional arguments passed to other functions, including \code{\link[PopED]{start_parallel}}.}
}
\value{
NCA results and diagnostic test results
}
\description{
\pkg{ncappc} is a flexible tool, to
\enumerate{
 \item perform a traditional NCA
 \item perform simulation-based posterior predictive checks for a
population PK model using NCA metrics.
}
}
\details{
Non-compartmental analysis (NCA) calculates pharmacokinetic (PK) metrics 
related to the systemic exposure to a drug following administration, e.g. 
area under the concentration-time curve and peak concentration. \pkg{ncappc} 
performs a traditional NCA using the observed plasma concentration-time data.
In the presence of simulated plasma concentration-time data, \pkg{ncappc} 
also performs simulation-based posterior predictive checks (ppc) using NCA 
metrics for the corresponding population PK (PopPK) model used to generate 
the simulated data. The diagnostic analysis is performed at the population as
well as the individual level. The distribution of the simulated population 
means of each NCA metric is compared with the corresponding observed 
population mean. The individual level comparison is performed based on the 
deviation of the mean of any NCA metric based on simulations for an 
individual from the corresponding NCA metric obtained from the observed data.
Additionally, \pkg{ncappc} reports the normalized prediction distribution 
error (NPDE) of the simulated NCA metrics for each individual and their 
distribution within a population. \pkg{ncappc} produces two default outputs 
depending on the type of analysis performed, i.e., traditional NCA and PopPK 
diagnosis. The PopPK diagnosis feature of \pkg{ncappc} produces 7 sets of 
graphical outputs to assess the ability of a population model to simulate the
concentration-time profile of a drug and thereby identify model 
misspecification. In addition, tabular outputs are generated showing the 
values of the NCA metrics estimated from the observed and the simulated data,
along with the deviation, NPDE, regression parameters used to estimate the 
elimination rate constant and the related population statistics. The default 
values of the arguments used in \pkg{ncappc} are shown in the \strong{Usage}
section of this document and/or in \strong{bold} in the \strong{Arguments} 
section.
}
\examples{
out <- ncappc(obsFile=system.file("extdata","pkdata.csv",package="ncappc"), 
  onlyNCA = TRUE,
  extrapolate = TRUE,
  printOut = FALSE,
  evid = FALSE,
  psnOut=FALSE)
  

data_1 <- data.frame(
  ID=1,
  TIME = c(0,0.25,0.5,1,1.5,2,3,4,6,8,12,16,24),
  DV=c(0, 0.07, 0.14, 0.21, 0.24, 0.27, 0.26, 0.25, 0.22, 0.19, 0.13, 0.081, 0.033)
)
out_1 <- ncappc(obsFile=data_1,
                onlyNCA = TRUE,
                extrapolate = TRUE,
                printOut = FALSE,
                evid = FALSE,
                timing=TRUE)


data_2 <- dplyr::filter(data_1,TIME>17|TIME<3)
out_2 <- ncappc(obsFile=data_2,
                onlyNCA = TRUE,
                extrapolate = TRUE,
                printOut = FALSE,
                evid = FALSE,
                force_extrapolate=TRUE)

}
