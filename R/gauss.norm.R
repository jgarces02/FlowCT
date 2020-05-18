#'gauss.norm
#'
#' It performs a Gaussian normalization based on \code{\link[flowStats:gaussNorm]{flowStats::gaussNorm()}} function.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param assay.i Name of matrix stored in the \code{FCS.SE} object from which calculate correlation. Default = \code{"transformed"}.
#' @param marker.to.norm Marker name(s) to normalize (based on visual inspection with \code{\link[FlowCT:multidensity]{FlowCT::multidensity()}}).
#' @param norm.matrix.name New normalized matrix name (it will stored within the \code{FCS.SE} object). Default = \code{"normalized"}.
#' @keywords Gaussian normalization
#' @keywords marker alignment
#' @keywords marker normalization
#' @export
#' @examples
#' \dontrun{
#'  fcs_se <- gauss.norm(fcs.SE = fcs_se,marker.to.norm = c("CD62L", "CCR4", "SSC_A"))
#' }

gauss.norm <- function(fcs.SE, assay.i = "transformed", marker.to.norm, norm.matrix.name = "normalized"){
  data <- as.flowSet.SE(fcs.SE, assay.i = assay.i)

  for(i in marker.to.norm) data <- flowStats::gaussNorm(as.flowSet.SE(data, assay.i = "transformed"), i)$flowset
  SummarizedExperiment::assay(fcs.SE, i = norm.matrix.name) <- t(flowCore::fsApply(data, exprs))
  return(fcs.SE)
}
