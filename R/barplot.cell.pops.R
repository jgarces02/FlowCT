#' barplot.cell.pops
#'
#' This function calculates cluster proportions (or raw counts) for each identified cluster and plot them on a stacked barplot.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param assay.i Name of matrix stored in the \code{FCS.SE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param cell.clusters A vector with clusters identified through \code{\link[FlowCT:fsom.clustering]{FlowCT::fsom.clustering()}} (and, normaly, later renamed).
#' @param plot Logical indicating whether plotting stacked barplot. Default = \code{TRUE}.
#' @param count.by Variable name (from \code{colData(fcs.SE)}) for calculating proportions (or counts) and drawing the x-axis in the stacked bar plotting.
#' @param facet.by Variable name (from \code{colData(fcs.SE)}) for splitting the stacked bar plotting. Default = \code{NULL}.
#' @param return.mode String for specifying if final resuls should be proportions ("percentage", default) or raw counts ("counts").
#' @keywords proportions
#' @keywords barplot
#' @export barplot.cell.pops
#' @method barplot cell.pops
#' @import ggplot2
#' @importFrom SummarizedExperiment colData assay
#' @examples
#' \dontrun{
#' prop_table <- barplot.cell.pops(fcs.SE = fcs_se, cell.clusters = fcs_se$SOM_named, 
#'     count.by = "sample_id", facet.by = "condition", 
#'     return.mode = "percentage")
#' counts_table <- barplot.cell.pops(fcs.SE = fcs_se, cell.clusters = fcs_se$SOM_named, 
#'     count.by = "condition", return.mode = "counts")
#' }

barplot.cell.pops <- function(fcs.SE, assay.i = "normalized", cell.clusters, plot = T, count.by, facet.by = NULL, return.mode = "percentage"){
  data <- t(assay(fcs.SE, i = assay.i))
  metadata <- colData(fcs.SE)
  colors_palette <- div.colors(length(unique(cell.clusters)))
  
  counts_table <- table(cell.clusters, metadata[,count.by])
  prop_table <- prop.table(counts_table, margin = 2)*100
  
  if(plot){
    ggdf <- data.table::melt(prop_table, value.name = "proportion")
    colnames(ggdf)[2] <- count.by
    
    mm <- match(ggdf[,count.by], metadata[,count.by]) #add other infos
    ggdf <- data.frame(metadata[mm,], ggdf)
    
    g <- ggplot(ggdf, aes_string(x = count.by, y = "proportion", fill = "cell.clusters")) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_fill_manual(values = colors_palette)
    if(is.null(facet.by)) print(g) else print(g + facet_wrap(~ eval(parse(text = facet.by)), scales = "free_x"))
  }
  if(return.mode == "percentage"){
    return(prop_table)
  }else if(return.mode == "counts"){
    return(counts_table)
  }else{
    cat("Please, specify a valid 'return.mode' value, i.e.: counts or percentage")
  }
}