#' median.heatmap
#'
#' This function draws a heatmap with median values for each FCS file or for identified cluster with \code{\link[FlowCT:fsom.clustering]{FlowCT::fsom.clustering()}}
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param assay.i Name of matrix stored in the \code{FCS.SE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param cell.clusters A vector with clusters identified through \code{\link[FlowCT:fsom.clustering]{FlowCT::fsom.clustering()}} (and, normaly, later renamed). Default = \code{NULL} (i.e., median values will be calculated for each FCS file).
#' @param markers.to.use Vector with markers to use. Default = \code{"all"}.
#' @param not.metadata Vector with variable names (from \code{colData(fcs.SE)}) for not including in the heatmap annotation. Default = \code{"filename"}.
#' @keywords heatmap
#' @keywords cell cluster percentages
#' @keywords median expression values
#' @export median.heatmap
#' @import dplyr
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @examples
#' \dontrun{
#' # option 1: general heatmap (by FCS file)
#' median.heatmap(fcs_se, not.metadata = c("sample_id", "file_name"))
#' 
#' # option 2: heatmap by SOM-identified clusters
#' median.heatmap(fcs.SE = fcs_se, assay.i = "normalized", cell.clusters = fcs_se$SOM)
#' }

median.heatmap <- function(fcs.SE, assay.i = "normalized", cell.clusters = NULL, markers.to.use = "all", not.metadata = "filename"){
  data <- t(assay(fcs.SE, i = assay.i))
  metadata <- fcs.SE@metadata$reduced_metadata
  
  if(markers.to.use == "all") markers.to.use <- colnames(data)
  
  ## prepare median tables
  if(is.null(cell.clusters)){
    med <- median.values(fcs.SE, assay.i = assay.i)
  }else{
    expr_median <- data.frame(cell_clusters = cell.clusters, data[,markers.to.use]) %>%
      group_by(.data$cell_clusters) %>% summarize_all(list(stats::median)) %>% as.data.frame(.data)
    
    expr_saturated_median <- data.frame(cell_clusters = cell.clusters, exprs.saturate(data[,markers.to.use])) %>%
      group_by(.data$cell_clusters) %>% summarize_all(list(stats::median)) %>% as.data.frame(.data)
  }
  
  ## heatmap
  if(is.null(cell.clusters)){
    annotation_colors <- col.annot.pheatmap(metadata[,!(colnames(metadata) %in% not.metadata), drop = F])
    color <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)
    
    print(pheatmap::pheatmap(t(med[,markers.to.use]), color = color, display_numbers = FALSE,
                   number_color = "black", fontsize_number = 5, clustering_method = "average",
                   annotation = metadata[,!(colnames(metadata) %in% not.metadata), drop = F], 
                   annotation_colors = annotation_colors, 
                   show_colnames = F))
  }else{
    ## calculate cluster frequencies
    clustering_table <- table(cell.clusters)
    clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
    labels_row <- paste0(expr_saturated_median$cell_clusters, " (", clustering_prop ,"%)")
    
    d <- stats::dist(expr_median[,markers.to.use], method = "euclidean")
    cluster_rows <- stats::hclust(d, method = "average")
    expr_heat <- as.matrix(expr_saturated_median[,markers.to.use])
    # rownames(expr_heat) <- paste0("c.", expr_saturated_median$cell_clusters) #force rownames to not crash heatmap (??)
    rownames(expr_heat) <- rownames(expr_saturated_median) #force rownames to not crash heatmap (??)
    
    ## annot colors
    annot_row <- expr_saturated_median[,"cell_clusters", drop = F]
    # rownames(annot_row) <- paste0("c.", rownames(annot_row)) #force rownames to not crash heatmap (??)
    annotation_colors <- col.annot.pheatmap(expr_saturated_median[,"cell_clusters", drop = F])
    color <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)
    legend_breaks <- seq(from = 0, to = 1, by = 0.1)
    
    print(pheatmap::pheatmap(expr_heat, color = color, annotation_legend = F,
                   cluster_cols = FALSE, cluster_rows = cluster_rows, labels_row = labels_row,
                   display_numbers = FALSE, number_color = "black",
                   fontsize = 8, fontsize_number = 6, legend_breaks = legend_breaks,
                   annotation_row = annot_row, annotation_colors = annotation_colors))
  }
}