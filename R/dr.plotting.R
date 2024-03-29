#' Plot DR data
#'
#' This function plots the indicated dimensional reduction (DR) from a previously calculated \code{\link[FlowCT:dim.reduction]{FlowCT::dim.reduction()}} object.
#' @param data A object with DR generated with \code{\link[FlowCT:dim.reduction]{FlowCT::dim.reduction()}} or a \code{data.frame} with DR, expression and metadata information (like the first element list of the object generated with \code{\link[FlowCT:dim.reduction]{FlowCT::dim.reduction()}}).
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param plot.dr String indicating the desired DR to plot (this indicated DR should be prevoulsy calculated to being plotted). If `pca.loadings`, the plotting result will show the weight of each marker for each PC (the first two components).
#' @param n.dims Vector indicating the two DR components to plot. Default = \code{c(1,2)} (by now, these are the only dims allowed).
#' @param color.by Variable from (from \code{colData(fcs.SCE)}) for dots coloring. If \code{color.by = "expression"} (default), plot will be automatically splitted for each marker.
#' @param shape.by Variable from (from \code{colData(fcs.SCE)}) for dots shaping. Default = \code{NULL}.
#' @param facet.by Variable from (from \code{colData(fcs.SCE)}) for plot spliting. Default = \code{NULL}.
#' @param omit.markers Vector with markers to omit when plotting with \code{color.by = "expression"}. Default = \code{NULL}.
#' @param title Title to add to the plot.
#' @param label.by Variable from (from \code{colData(fcs.SCE)}) for dots labeling. Default = \code{NULL}.
#' @param size Point size. Default = \code{1}.
#' @param raster Number of pixels to consider for rastering image, if none indicated, not rasterization will performed. It is based on \href{https://github.com/exaexa/scattermore}{\code{scattermore} package}.
#' @param return.df Logical indicating if built \code{data.frame} with DR information and metadata must be returned. Default = \code{FALSE}.
#' @param colors Vector with colors for plotting. Default = \code{NULL} (i.e., it will choose automatically a vector of colors according to \code{\link[FlowCT:div.colors]{FlowCT::div.colors()}}).
#' @keywords dimensional reduction plotting
#' @keywords tSNE
#' @keywords PCA
#' @keywords UMAP
#' @export
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#' @importFrom data.table melt as.data.table
#' @importFrom SingleCellExperiment colData reducedDims
#' @importFrom SummarizedExperiment assay
#' @importFrom RColorBrewer brewer.pal
#' @examples
#' \dontrun{
#' dr.plotting(fcs, plot.dr = "tSNE", color.by = "condition")
#' dr.plotting(fcs, plot.dr = "UMAP", color.by = "patient_id")
#' dr <- dr.plotting(fcs, plot.dr = "PCA", color.by = "SOM", facet.by = "condition", return.df = T)
#' }

dr.plotting <- function(data, assay.i = "normalized", plot.dr, color.by = "expression", shape.by = NULL, facet.by = NULL, omit.markers = NULL, title = "", label.by = NULL, size = 1, raster = F, return.df = F, colors, num.loadings = 5){
  if(class(data)[1] == "SingleCellExperiment"){
    if(plot.dr == "pca.loadings"){
      x <- data@metadata$PCA.rotation
    }else{
      pos <- match(tolower(plot.dr), tolower(names(data@int_colData@listData$reducedDims)))
    if(is.na(pos)) stop('The DR indicated has not been calculated yet or is differently named (please, check the output of reducedDimNames(data) to see the correct DR name).\n', call. = F)
    
    dr_calculated <- names(data@int_colData@listData$reducedDims)[pos]
    dr <- data@int_colData@listData$reducedDims@listData[[dr_calculated]]
    
    no.omit.markers <- rownames(data)[!(rownames(data) %in% omit.markers)]
    drmd <- as.data.frame(cbind(colData(data), dr, t(assay(data, i = assay.i))[,no.omit.markers]))
    }
  }else{
    if(tolower(plot.dr) != "pca.loadings"){
      drmd <- data
      dr <- data[,grep(tolower(plot.dr), names(data))]
    }else{
      x <- data$PCA.rotation
    }
  }
  
  if(tolower(plot.dr) == "pca.loadings"){
    dims = c(1,2)

    top_markers <- data.table::melt(x) %>% filter(Var2 %in% paste0("PC", dims)) %>%
      group_by(Var2) %>% arrange(desc(abs(value))) %>% slice(1:num.loadings)

    var_loadings <- x %>% as_tibble(rownames = "marker") %>% filter(marker %in% top_markers$Var1)

    g <- ggplot(var_loadings) + 
      geom_segment(aes_string(x = 0, y = 0, xend = paste0("PC", dims[1]), yend = paste0("PC", dims[2])), 
                     arrow = arrow(length = unit(0.1, "in")), colour = "brown") +
        ggrepel::geom_text_repel(data = var_loadings, aes_string(x = paste0("PC", dims[1]), y = paste0("PC", dims[2]), label = "marker")) +
        labs(x = "PCA-1", y = "PCA-2") + ggtitle("PCA (loadings)") +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom",
                panel.background = element_blank(), panel.border = element_rect(color = "black", fill = NA))

    return(g)
  }else{
    if(color.by != "expression"){
    if(missing(colors)) colors <- div.colors(length(unique(drmd[,color.by])))
    g <- ggplot(drmd, aes_string(x = colnames(dr)[1], y = colnames(dr)[2], color = color.by)) +
      labs(x = paste0(toupper(substr(colnames(dr)[1], 1, nchar(colnames(dr)[1])-1)), "-1"), 
           y = paste0(toupper(substr(colnames(dr)[2], 1, nchar(colnames(dr)[2])-1)), "-2")) + ggtitle(title) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom",
            panel.background = element_blank(), panel.border = element_rect(color = "black", fill = NA))
    
    if(!raster){
      g <- g + geom_point(aes_string(color = color.by), size = size)
    }else{
      if(!requireNamespace("scattermore", quietly = TRUE)) stop("Package \"scattermore\" (https://github.com/exaexa/scattermore) needed for this function to work. Please install it.", call. = FALSE)
      
      g <- g + scattermore::geom_scattermore(pointsize = size, pixels = rep(raster, 2)) #devtools::install_github('exaexa/scattermore')      
    }
    
    if(is.numeric(drmd[,color.by])){
      g <- g + scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50), name = color.by)
    }else{
      g <- g + scale_color_manual(values = colors, name = color.by) +
        guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
    }
    
    if(!is.null(facet.by)) g <- g + facet_wrap(~ eval(parse(text = facet.by)))
    
    if(!is.null(label.by)){
      g <- g + ggrepel::geom_label_repel(aes_string(label = label.by), nudge_y = 0.05)
    }
  }else{
    if (!requireNamespace("cowplot", quietly = TRUE)) stop("Package \"cowplot\" needed for this function to work. Please install it.", call. = FALSE)
    
    drmd <- as.data.frame(melt(as.data.table(drmd), measure.vars = no.omit.markers, value.name = "expression", variable.name = "antigen"))
    glist <- lapply(unique(drmd$antigen), function(x){
      g <- ggplot(drmd[drmd$antigen == x,], aes_string(x = colnames(dr)[1], y = colnames(dr)[2], color = "expression")) +
        labs(x = NULL, y = NULL) + ggtitle(x) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom",
              panel.background = element_blank(), panel.border = element_rect(color = "black", fill = NA))
      
      if(!raster){
        g + geom_point(aes_string(color = color.by), size = size) +
          scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50), name = NULL)
      }else{
        g + scattermore::geom_scattermore(pointsize = size, pixels = rep(raster, 2)) +
          scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50), name = NULL)
      }
      
    })
    g <- do.call(cowplot::plot_grid, glist)
    
    ## add general title
    title <- cowplot::ggdraw() + cowplot::draw_label(plot.dr, fontface = 'bold', x = 0, hjust = 0) +
      theme(plot.margin = margin(0, 0, 0, 7))
    g <- cowplot::plot_grid(title, g, ncol = 1, rel_heights = c(0.1, 1.5))
  }
  
  if(return.df) return(list(data = drmd, plot = g)) else return(g)
  }
}
