% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flowplot.R
\name{flowplot}
\alias{flowplot}
\title{Draw expression dotplots}
\usage{
flowplot(
  fcs.SCE,
  assay.i = "normalized",
  x.axis,
  y.axis,
  densities = T,
  color.by,
  select.values.color = "all",
  size = 0.5,
  colors = NULL
)
}
\arguments{
\item{fcs.SCE}{A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{FlowCT::fcs.SCE()}}. By default, the matrix used is the \code{arcsinh} transformed.}

\item{assay.i}{Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"normalized"}.}

\item{x.axis}{Vector of markers to draw x-axis on the dotplots. Length must be the same of \code{y.axis} because each marker has its parter in \code{y.axis}.}

\item{y.axis}{Vector of markers to draw y-axis on the dotplots. Length must be the same of \code{x.axis} because each marker has its parter in \code{x.axis}.}

\item{densities}{Logical indicating if densities must be drawn in each plot. Default = \code{TRUE}.}

\item{color.by}{Variable from (from \code{colData(fcs.SCE)}) for events coloring.}

\item{select.values.color}{Vector of values taken from \code{color.by} option to include in the coloring. Default = \code{"all"}.}

\item{size}{Point (event) size. Default = \code{0.5}.}

\item{colors}{Vector with colors for plotting. Default = \code{NULL} (i.e., it will choose automatically a vector of colors according to \code{\link[FlowCT:div.colors]{FlowCT::div.colors()}}).}
}
\description{
It draws the classical flow cytometry dotplot faceting two markers. Multiple dot plots can be plotted if \code{x.axis} and \code{y.axis} are specified with multiple markers.
}
\examples{
\dontrun{
flowplot(fcs, x.axis = c("FSC_A", "CD4", "CD62L"), y.axis = c("SSC_A", "CD8", "CD45"),
    color.by = "SOM", select.values.color = 1:10)
}
}
\keyword{dotplot}
\keyword{marker}
\keyword{scatterplot}
