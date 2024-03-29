% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unify.FCSheaders.R
\name{unify.FCSheaders}
\alias{unify.FCSheaders}
\title{Unify all FCS header for multiple FCS}
\usage{
unify.FCSheaders(
  filelist = NULL,
  directory = getwd(),
  pattern = ".fcs$",
  dataset = 1,
  fix = T,
  select.freq = 1,
  events = "all",
  view = F
)
}
\arguments{
\item{filelist}{A vector with full path of FCS files to be read, commonly generated through \code{\link[base:list.files]{base::list.files()}}. If \code{NULL}, this file list will be generated as indicated below.}

\item{directory}{If \code{filelist = NULL}, those files stored in this location will be read. Default = \code{getwd()} (current directory).}

\item{pattern}{Pattern for reading files within \code{directory}. Default = \code{"fcs"}.}

\item{dataset}{Specify data segment to be read in the FCS file (same behaviour as \code{\link[flowCore:read.FCS]{flowCcre::read.FCS()}}). Default = \code{1}.}

\item{fix}{Logical indicating if user wants to perform the unification of FCS's headers or only display header frecuencies. If \code{TRUE} (default value), you should to specify the chosen frequency.}

\item{select.freq}{Of diplayed frequencies (row-numered), indicate which is your option for unifiying FCS files. Default = \code{1} (i.e., the more frequent).}

\item{events}{Numeric vector indicating how many events to read in each FCS file (same behaviour as \code{\link[flowCore:read.FCS]{flowCcre::read.FCS()}}), those are the number of events for new and unified generated FCS files. Default = \code{NULL} (i.e., all events will be read).}

\item{view}{Logical indicating whether table with frequencies should be displayed in the terminal or in a external table. Default = \code{FALSE}.}
}
\value{
The final output if \code{fix = F} is a table with three columns: \code{names} for channels, \code{freq} for frequency of apparition of these channel names and \code{length} for number of channels within that frequency.

If \code{fix = TRUE}, those FCS files with a distinct header from selected frequency will be renamed/reordered, added with the suffix \code{fixed.fcs} and the original ones (unchanged) will be stored in a new \code{original_files} folder. In the case a FCS would have a different number of channels, it will moved to a new folder called \code{discarded_files_because_diffs} and dicarded from downstream analysis (this will be changed in the future).
}
\description{
It checks if FCS files within a specific folder, or indicated in a vector, have the same header (i.e., same channels:markers nomenclature and order) and offer the possibility to unify them creating new ones.
}
\examples{
\dontrun{
## detect header's frequencies
unify.FCSheaders(directory = "../data/", pattern = "fcs", fix = F)
 #                                                                            names freq length
 # 1 FSC-H, SSC-H, FITC-A, PE-A, PerCP-Cy5-5-A, APC-A, APC-H7-A, V500-A, V450-A       12     10
 # 2 FSC-H, SSC-H, FITC-A, PE-A, PerCP-C5-A, APC-A, APC-H7-A, Pac-Orange, Pac-Blue    10     10
## unify FCS's headers
unify.FCSheaders(directory = "../data/", pattern = "fcs", fix = T, select.freq = 2)
}
}
\keyword{FCS}
\keyword{differing}
\keyword{files}
\keyword{headers}
\keyword{unify}
