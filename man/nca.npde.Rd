\name{nca.npde}
\alias{nca.npde}
\title{Calculates individual normalized prediction distribution errors (NPDE) from
PDE data.}
\usage{
nca.npde(pdedata, pdecol)
}
\arguments{
  \item{pdedata}{A data frame containing the prediction
  distribution errors (PDE) of NCA metrics}

  \item{pdecol}{The range of column numbers in the data
  frame containing the PDE values, which will be used to
  calculate the corresponding NPDE}
}
\value{
returns the data frame with the NPDE values based on the
input data.
}
\description{
\pkg{nca.npde} calculates individual normalized prediction
distribution errors (NPDE) of selected NCA metrics from the
PDE data.
}
\details{
\pkg{nca.npde} calculates individual normalized prediction
distribution errors (NPDE) of selected NCA metrics from PDE
data. The The deviation of each estimated NCA metrics is
scaled by the "spread" of the simulated values. By default,
this function calculates the NPDE values of all columns of
the input data frame.
}

