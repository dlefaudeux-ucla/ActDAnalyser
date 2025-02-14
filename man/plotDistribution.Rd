% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ploting.R
\name{plotDistribution}
\alias{plotDistribution}
\title{Plot distribution of half-life across experiments}
\usage{
plotDistribution(HL, metadata, r2_th = NULL, control = NULL)
}
\arguments{
\item{HL}{list. half-life list derived from \code{\link{convertSlopeToHL}}
function.}

\item{metadata}{data.frame. Metadata of the data. Must have specific column
names see \code{\link{checkMetadata}} function.}

\item{r2_th}{vector of 2 elements. Interval of adjusted R² to select genes
for ploting. Default to NULL: will not do any selection on adjusted R² values.}

\item{control}{list. Further control ggplot2 output. This must be a named
list with the names corresponding to ggplot2 function and the associated
value is also a named list of the ggplot2 function argument name and its
associated value.}
}
\description{
This function plots the overall distribution of half-lives infered.
}
\examples{
attach(MEF_dataset)
plotDistribution(HL = halflife, metadata = metadata)
detach(MEF_dataset)
}
