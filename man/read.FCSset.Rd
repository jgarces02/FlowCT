% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.FCSset.R
\name{read.FCSset}
\alias{read.FCSset}
\title{Load multiple FCS files}
\usage{
read.FCSset(
  filelist = NULL,
  directory = getwd(),
  pattern = ".fcs$",
  events = "all",
  dataset = 1,
  num.threads = NULL
)
}
\arguments{
\item{filelist}{A vector with full path of FCS files to be read, commonly generated through \code{\link[base:list.files]{base::list.files()}}. If \code{NULL}, this file list will be generated as indicated below.}

\item{directory}{If \code{filelist = NULL}, those files stored in this location will be read. Default = \code{getwd()} (current directory).}

\item{pattern}{Pattern for reading files within \code{directory}. Default = \code{"fcs"}.}

\item{events}{Numeric vector indicating how many events to read in each FCS file (same behaviour as \code{\link[flowCore:read.FCS]{flowCcre::read.FCS()}}). Default = \code{NULL}, i.e., all events will be read.}

\item{dataset}{Specify data segment to be read in the FCS file (same behaviour as \code{\link[flowCore:read.FCS]{flowCcre::read.FCS()}}). Default = \code{1}.}

\item{num.threads}{Number of threads for reading the FCS files. For Windows users, you should to install \code{parallelsugar} package from (\href{https:https://github.com/nathanvan/parallelsugar}{Nathanvan's GitHub}). Defult = \code{NULL} (i.e., not parallelization).}
}
\description{
It reads multiple FCS files (in a computationaly-distributed way) contained in a specific folder or given through a vector. Important note: all files should have an identicar header (i.e., same name and markers order), if not only those common will be read (you can use \code{\link[FlowCT:unify.FCSheaders]{FlowCT::unify.FCSheaders()}} for doing that).
}
\examples{
\dontrun{
# option 1: trough a vector with filenames (full path)
filelist <- list.files(pattern = "fcs", path = "../data/", full.names = T)
fcs1 <- read.FCSset(filelist = filelist, which.lines = 1000)

# option 2: specifiying a directory
fcs2 <- read.FCSset(directory = "../data", pattern = ".LMD", num.threads = 4)
}
}
\keyword{FCS}
\keyword{parallel}
\keyword{reading}
