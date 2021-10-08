#' Combine multiple subclustering with initial one
#'
#' It combines initial \code{fcs.SCE} object (without subclustering) with other \code{fcs.SCE} objects with subclustering analysis coming from downstream steps and generates a new \code{fcs.SCE} object. This final \code{fcs.SCE} object has an additional column combining all information from initial and subclustering analysis.
#' @param global.fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{FlowCT::fcs.SCE()}}. The initial one, without extracting any cell population.
#' @param subclustering.fcs.SCE A list with all \code{fcs.SCE} object(s) generated in the subclustering analysis (they have to come from the original \code{initial.fcs.SCE}).
#' @param clusters.named Column names from the \code{global.fcs.SCE} and \code{global.fcs.SCE} objects which contains renamed clusters (through \code{\link[FlowCT:clusters.rename]{FlowCT::clusters.rename()}}).
#' @keywords final fcs.SCE object
#' @keywords combine subclustering
#' @export
#' @importFrom SingleCellExperiment colData
#' @import dplyr
#' @examples
#' \dontrun{
#'  fcs_final <- combine.subclusterings(global.fcs.SCE = fcs, 
#'    clusters.named = c("SOM_named", "SOM_named_lymhpos", "SOM_named_monos"), 
#'    subclustering.fcs.SCE = list(fcs_lymphos, fcs_monos))
#' }

combine.subclusterings <- function(global.fcs.SCE, subclustering.fcs.SCE, clusters.named){
  subclusterings <- unique(unlist(lapply(subclustering.fcs.SCE, function(x) colData(x)[[clusters.named[1]]])))

  ## combine subclusterings
  subs <- unlist(lapply(1:length(subclustering.fcs.SCE), function(x){
    aux <- colData(subclustering.fcs.SCE[[x]])[clusters.named[1+x]]
    setNames(aux[[1]], rownames(aux)) #pop/id labelling
  }))

  # ## extract removed populations and substract them from global, v1 (metadata, deprecated)
  # rems <- lapply(subclustering.fcs.SCE, function(x) x@metadata$removed_populations)
  # names(rems) <- subclusterings
  # rems2 <- unlist(rems, use.names = F)
  
  # global.fcs.SCE <- global.fcs.SCE[,!(colnames(global.fcs.SCE) %in% rems2)] #time-consuming

  ## extract removed populations and substract them from global, v2 (no metadata)
  rems3 <- mapply(function(x, y) { #diffs between global and subclustered ones
    aux <- global.fcs.SCE[,global.fcs.SCE[[clusters.named[1]]] == x]
    setdiff(colnames(aux), colnames(y))
    }, subclusterings, subclustering.fcs.SCE)

  global.fcs.SCE <- global.fcs.SCE[,!(colnames(global.fcs.SCE) %in% unlist(rems3))] #time-consuming

  ## final clustering
  global.fcs.SCE$final_clustering <- as.character(global.fcs.SCE[[clusters.named[1]]])
  global.fcs.SCE$final_clustering[match(names(subs), global.fcs.SCE$cell_id)] <- as.character(subs)

  global.fcs.SCE$final_clustering <- factor(global.fcs.SCE$final_clustering)

  ## add subclustering information into metadata
  # global.fcs.SCE@metadata$subclusterings$populations <- paste(subclusterings, collapse = " + ")
  # global.fcs.SCE@metadata$subclusterings$removed_populations <- rems
  global.fcs.SCE@metadata$combined_subclusterings <- paste(subclusterings, collapse = " + ")

  return(global.fcs.SCE)
}
