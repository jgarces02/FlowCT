## INFORMATION ##########################################################################
# FlowCT version ---> 2.0
# backbone script version ---> 1.0
#### Changes:
# -> Add maker.names for renaming marker names
##########################################################################################

### Environment setting ##################################################################
## Load packages and functions
# Sys.setenv(http_proxy  = "http://proxy.unav.es:8080")
# Sys.setenv(https_proxy = "http://proxy.unav.es:8080")
# devtools::install_github("jgarces02/FlowCT@v2.0", auth_token = "21ea9880f944d42755479e54a5b19ddd00fe17f6")
library(FlowCT.v2)
library(SummarizedExperiment)

## Working directory
setwd("G:/Mi unidad/Proyectos/FlowCT_Ciro/FlowCT.v2/results/")

### unify all FCS nomenclatures ##########################################################
unify.FCSheaders(directory = "../data/", pattern = "fcs", fix = F)
unify.FCSheaders(directory = "../data/", pattern = "fcs", fix = T, select.freq = 2)


### Prepare metadata and fcs.SE object ###################################################
## metadata matrix construction
(filenames <- list.files(pattern = "fcs", path = "../data/"))

md <- data.frame(filename = filenames, #mandatory
                 sample_id = 1:length(filenames),
                 condition = sapply(filenames,function(x) strsplit(x, split = "_|\\.")[[1]][1]),
                 patient_id = sapply(filenames,function(x) strsplit(x, split = "_|\\.")[[1]][2]))
(md)

## FCS reading and transforming
fcs_se <- fcs.SE(directory = "../data/", events = 1000, metadata = md, transf.cofactor = 500)

## adjust maker names
marker.names(fcs_se)
new_names <- c("FSC_A", "FSC_H", "SSC_A", "SSC_H", "CD62L", "CXCR3", "CD8", "CCR4", "CCR6", "CD4", "CD45", "CD27")
marker.names(fcs_se, new.names = new_names)

## select markers for downstrean analysis
surface_markers <- c("CD62L", "CXCR3", "CD8", "CCR4", "CCR6", "CD4", "CD45", "CD27")
physical_markers <- c("FSC_A", "SSC_A")


### QC and doublets removal ##############################################################
fcs_se <- qc.and.removeDoublets(fcs_se, physical.markers = c("FSC_A", "FSC_H", "SSC_A", "SSC_H"))


### Normalization and data alignment #####################################################
multidensity(fcs.SE = fcs_se, assay.i = "transformed", subsampling = 1000)

fcs_se <- gauss.norm(fcs.SE = fcs_se, marker.to.norm = c("CD62L", "CCR4", "SSC_A"))

multidensity(fcs.SE = fcs_se, assay.i = "normalized", subsampling = 1000)
multidensity(fcs_se, assay.i = 2, ridgeline.lim = 10, color.by = "filename",
             show.markers = c("CD62L", "CD4"), interactive = T)


### Descriptive and exploratory analysis #################################################
## boxplot with cell numbers by condition
cell.count.bx(fcs.SE = fcs_se, x.axis = "condition") 

## dimensinal reduction and heatmap with median values
median.dr(fcs_se, color.by = "filename")
median.heatmap(fcs_se, not.metadata = c("sample_id", "filename"))

## PCA and heatmap on single cell expression
fcs_se100 <- sub.samples(fcs.SE = fcs_se, subsampling = 100)

dr_se <- dim.reduction(fcs_se100, dr.method = "tSNE")
dr.plotting(dr_se, dr.calculated = "tSNE", color.by = "condition")

sc.heatmap(fcs_se100, subsampling = NULL, markers.to.use = surface_markers)


### Clustering ###########################################################################
## FlowSOM
fsom <- fsom.clustering(fcs.SE = fcs_se, markers.to.use = surface_markers, k.metaclustering = 40)
fcs_se$SOM <- fsom$metaclusters

median.heatmap(fcs.SE = fcs_se, assay.i = "normalized", cell.clusters = fcs_se$SOM)

## dimensional reduction
fcs_se1000 <- sub.samples(fcs.SE = fcs_se, subsampling = 1000)

dr <- dim.reduction(fcs_se1000, dr.method = "all")

dr.plotting(dr, dr.calculated = "tSNE", color.by = "condition")
dr.plotting(dr, dr.calculated = "UMAP", color.by = "patient_id")
dr.plotting(dr, dr.calculated = "PCA", color.by = "SOM", facet.by = "condition")


## FCS exporting and cluster analysis
## Generate FCS files with dr data
export.metaFCS(fcs.SE = fcs_se1000, dr.object = dr, output.name = "clust1.fcs")

# ++++++++++++++++++++++++++++
# external step >>> manual analysis in a flow cytometry software to identify clusters >>> Excel file
# ++++++++++++++++++++++++++++

## rename and merge clusters according external analysis
replacedata <- readxl::read_excel("../data/Th_new_clustering.xlsx", col_types = "text")

dr$dr_melted$SOM_named <- clusters.rename(dr$dr_melted$SOM, cluster = replacedata$original_cluster, name = replacedata$new_cluster)
dr$dr$SOM_named <- clusters.rename(dr$dr$SOM, cluster = replacedata$original_cluster, name = replacedata$new_cluster)
fcs_se1000$SOM_named <- clusters.rename(fcs_se$SOM, cluster = replacedata$original_cluster, name = replacedata$new_cluster)
fsom$plotStars_value_named <- clusters.rename(fsom$plotStars_value, cluster = replacedata$original_cluster, name = replacedata$new_cluster)

## draw PCA, tSNE, MST and heatmap colored by new merged clusters
dr.plotting(dr, dr.calculated = "UMAP", color.by = "SOM_named")

PlotStars(fsom$fsom, backgroundValues = fsom$plotStars_value_named,  
  backgroundColor = alpha(div.colors(40), alpha = 0.4))

median.heatmap(fcs.SE = fcs_se1000, assay.i = "normalized", cell.clusters = fcs_se1000$SOM_named)


### Subclustering ########################################################################
## prepare fcs.SE object
fcs_seL <- fcs_se1000[,fcs_se1000$SOM_named == "lymphocytes"]
metadata(fcs_seL)$subclustering <- "lymphocytes"

## exploratory analysis
cell.count.bx(fcs_seL, assay.i = "normalized", x.axis = "condition")

median.heatmap(fcs_seL, not.metadata = c("sample_id", "filename"))

## FlowSOM
fsomL <- fsom.clustering(fcs.SE = fcs_seL, markers.to.use = surface_markers, k.metaclustering = 40)
fcs_seL$SOM_L <- fsomL$metaclusters

median.heatmap(fcs.SE = fcs_seL, assay.i = "normalized", cell.clusters = fcs_seL$SOM_L)

## DR
drL <- dim.reduction(fcs_seL, dr.method = "all", markers.to.use = surface_markers)

dr.plotting(drL, dr.calculated = "tSNE")
dr.plotting(drL, dr.calculated = "tSNE", color.by = "SOM_L")

## Generate FCS files with subclustering data
export.metaFCS(fcs.SE = fcs_seL, dr.object = drL, output.name = "subclust1.fcs")

# ++++++++++++++++++++++++++++
# external step >>> manual analysis in a flow cytometry software to identify clusters >>> Excel file
# ++++++++++++++++++++++++++++

## rename and merge clusters according external analysis
replacedataL <- readxl::read_excel("../data/Th_new_clusteringL.xlsx", col_types = "text")

drL$dr_melted$SOM_L_named <- clusters.rename(drL$dr_melted$SOM, cluster = replacedataL$original_cluster, name = replacedataL$new_cluster)
drL$dr$SOM_L_named <- clusters.rename(drL$dr$SOM, cluster = replacedataL$original_cluster, name = replacedataL$new_cluster)
fcs_seL$SOM_L_named <- clusters.rename(fcs_seL$SOM, cluster = replacedataL$original_cluster, name = replacedataL$new_cluster)
fsomL$plotStars_value_named <- clusters.rename(fsomL$plotStars_value, cluster = replacedataL$original_cluster, name = replacedataL$new_cluster)

## draw PCA, tSNE, MST and heatmap colored by new merged clusters
dr.plotting(drL, dr.calculated = "tSNE", color.by = "SOM_L_named")

PlotStars(fsomL$fsom, backgroundValues = fsomL$plotStars_value_named,  
          backgroundColor = alpha(div.colors(40), alpha = 0.4))

median.heatmap(fcs.SE = fcs_seL, assay.i = "normalized", cell.clusters = fcs_seL$SOM_L_named)

## Circular hyerarchical clustering tree
circ.tree(fcs.SE = fcs_seL, cell.clusters = fcs_seL$SOM_L_named, nodes = "display")
circ.tree(fcs.SE = fcs_seL, cell.clusters = fcs_seL$SOM_L_named, nodes = c(17, 25), scale.size = 11)

## heatmap for single cells
sc.heatmap(fcs_seL, assay.i = "normalized", 
           not.metadata = c("cell_ID", "filename", "SOM", "SOM_named", "SOM_L"))


### cluster identif (misc.) ##############################################################
## differential dotplots
dotplot.DE(fcs.SE = fcs_se1000, markers.to.use = surface_markers)


### data exporting #######################################################################
## delete non-interesting populations
delete_pops <- "unclassified"
fcs_se1000_rm <- remove.pop(fcs_se1000, clusters.named = "SOM_named", population = delete_pops)
fcs_seL_rm <- remove.pop(fcs_seL, clusters.named = "SOM_L_named", population = delete_pops)

## combine clustering and subclustering
fcs_se_final <- combine.subclusterings(initial.fcs.SE = fcs_se1000_rm, clusters.named = "SOM_named",
                                       subclustering.fcs.SE = list(fcs_seL_rm))

## final FCS files with named clusters and DR info
dr_final <- dr$dr[dr$dr$cell_id %in% fcs_se_final$cell_id,]

export.metaFCS(fcs.SE = fcs_se_final, dr.object = dr_final, output.name = "prueba_final4C.fcs")

## cellular proportions (or raw counts)
prop_tableL <- barplot.cell.pops(fcs.SE = fcs_se_final, cell.clusters = fcs_se_final$SOM_named_final, 
                                 count.by = "sample_id", facet.by = "condition", 
                                 return.mode = "percentage")

dataset <- data.frame(md, as.data.frame.matrix(t(prop_tableL)))
write.table(dataset, file = "results.txt", sep = "\t")

## median values (from non-transformed data)
med <- median.values(fcs.SE = fcs_se_final)
# med <- median.values(fcs.SE = fcs_se_final, var = "SOM_named_final")
dataset2 <- data.frame(md, med)
write.table(dataset2, file = "results_median.txt", sep = "\t")

### statistics ###########################################################################
## Correlation between (two) conditions
corplot.conditions(fcs.SE = fcs_seL, cell.clusters = fcs_seL$SOM_L_named, 
                   condition.column = "condition")

corplot.conditions(fcs.SE = fcs_se_final, cell.clusters = fcs_se_final$SOM_named_final, 
                   condition.column = "condition")

## boxplot
boxplot.cell.clustering(fcs.SE = fcs_seL, cell.clusters = fcs_seL$SOM_L_named, facet = T,
                        return.stats = F, plot.only.sig = c(T, 0.1), facet.free.scale = "free")
sig <- boxplot.cell.clustering(fcs.SE = fcs_se_final, cell.clusters = fcs_se_final$SOM_named_final,
                               return.stats = T)

## Dumbbell plot
dumbPlot.cell.clustering(fcs.SE = fcs_se_final, cell.clusters = fcs_se_final$SOM_named_final, 
                         return.stats = F, condition.column = "condition")

## diffdots plot
diffdots.cell.clustering(fcs.SE = fcs_se_final, cell.clusters = fcs_se_final$SOM_named_final, 
                         return.stats = F, condition.column = "condition")
