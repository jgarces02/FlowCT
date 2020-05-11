# 'diffdots.cell.clustering
#'
#' It draws a differential dot plot (longitudinally) according condition for each cell cluster identified.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param assay.i Name of matrix stored in the \code{FCS.SE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param cell.clusters A vector with clusters identified through \code{\link[FlowCT:fsom.clustering]{FlowCT::fsom.clustering()}} (and, normaly, later renamed).
#' @param condition.column Column name from the \code{colData(fcs.SE)} object which contains condition information. De
#' @param psig.cutoff P-value cutoff. Default = \code{0.05}.
#' @param return.stats Logical indicating if calculated statistics should be returned in a new variable. Default = \code{TRUE}.
#' @keywords differential dotplot
#' @keywords Dumbbell plot
#' @keywords longitudinal dotplot
#' @export
#' @import ggplot2
#' @examples
#' \dontrun{
#' diffDots.cell.clustering(fcs.SE = fcs_se, cell.clusters = fcs_se$SOM_named, return.stats = F)
#' }

diffdots.cell.clustering <- function(fcs.SE, assay.i = "normalized", cell.clusters, condition.column, psig.cutoff = 0.05, return.stats = T){
  ## prepare tables
  prop_table <- as.data.frame.matrix(t(barplot.cell.pops(fcs.SE, cell.clusters, count.by = "filename", plot = F, assay.i = "normalized")))
  
  prop_table_md <- merge(fcs.SE@metadata$reduced_metada, prop_table, by.x = "filename", by.y = "row.names")
  
  dfm <- data.table::melt(prop_table_md, measure.vars = unique(cell.clusters))
  dfma <- stats::aggregate(dfm$value ~ dfm$variable + dfm[,condition.column], FUN = stats::median)
  colnames(dfma) <- c("variable", condition.column, "value")
  dfma <- transform(dfma, pct = log(stats::ave(dfma$value, dfma[,condition.column], FUN = function(x) x/sum(x)*100))) #transform to percentaje
  
  conditions <- unique(dfma[,condition.column])
  
  ## statistics table
  resultskw <- matrixTests::col_kruskalwallis(prop_table_md[,as.character(unique(cell.clusters))], prop_table_md[,condition.column])
  KWsig <- rownames(resultskw[resultskw$pvalue < psig.cutoff,])
  dfma$sig <- ifelse(dfma$variable %in% KWsig, "1", "0")
  
  if(length(unique(prop_table_md[,condition.column])) > 2){
    kw_posthoc <- list()
    for(i in KWsig) kw_posthoc[[i]] <- stats::pairwise.wilcox.test(prop_table_md[,i], prop_table_md[,condition.column])
  }
  
  
  ## plotting
  print(ggplot(dfma, aes_string(x = "condition", y = "pct", fill = "variable")) + 
          geom_point(size = 3, shape = 21, color = "gray63") + 
          # scale_color_manual(values = div.colors(30)) + 
          geom_line(aes_string(group = "variable", color = "sig"), size = 1) + 
          scale_colour_manual(values = c("gray63", "brown1"), labels = c("no sig.", "sig.")) + 
          # geom_label(data = subset(dfma, dfma$sig == 1 & dfma$condition == conditions[1]), 
          #            aes(x = 1, y = pct, label = variable, fill = variable), nudge_x = -0.1) +
          ggrepel::geom_label_repel(data = subset(dfma, dfma$sig == 1 & dfma$condition == conditions[1]),
                     aes_string(x = 1, y = "pct", label = "variable", fill = "variable"), 
                     nudge_x = -0.1, show.legend = F) +
          theme(panel.background = element_blank(), axis.line = element_line(color = "black")) + 
          labs(x = "\nCondition", y = "% of cells (log-transf.)\n", color = "", fill = "Cell clusters"))
  
  if(return.stats & length(unique(prop_table_md[,condition.column])) > 2){
    return(list(kruskal_results = resultskw, kw_posthoc_results = kw_posthoc))
  }else if(return.stats){
    return(resultskw)
  }
}